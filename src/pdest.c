/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Copyright (C) 2011 Individual developers
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include <pdest_p4est.h>
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <pdest_p8est.h>
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_communication.h>

int
pdest_object_is_valid (const pdest_object_t * ob)
{
  return
    (ob != NULL) &&
    sc_refcount_is_active (&ob->rc) &&
    (ob->dim == 2 || ob->dim == 3) &&
    (ob->typ != PDEST_TYPE_INVALID) && (ob->ject != NULL);
}

int
pdest_object_is_last (const pdest_object_t * ob)
{
  P4EST_ASSERT (pdest_object_is_valid (ob));

  return sc_refcount_is_last (&ob->rc);
}

void
pdest_object_ref (pdest_object_t * ob)
{
  P4EST_ASSERT (pdest_object_is_valid (ob));
  sc_refcount_ref (&ob->rc);
}

int
pdest_object_unref (pdest_object_t * ob)
{
  P4EST_ASSERT (pdest_object_is_valid (ob));
  if (sc_refcount_unref (&ob->rc)) {
    if (ob->jown) {
      /* only destroy payload if we own it */
      if (ob->dim == 2) {
#ifdef P4EST_BUILD_2D
        pdest_p4est_ject_destroy (ob->typ, ob->ject);
#else
        SC_ABORT_NOT_REACHED ();
#endif
      }
      else {
        P4EST_ASSERT (ob->dim == 3);
#ifdef P4EST_BUILD_3D
        pdest_p8est_ject_destroy (ob->typ, ob->ject);
#else
        SC_ABORT_NOT_REACHED ();
#endif
      }
    }

    if (ob->refed != NULL) {
      (void) pdest_object_unref (ob->refed);
    }
    P4EST_FREE (ob);
    return 1;
  }
  return 0;
}

pdest_object_t     *
pdest_object_new (int dim, pdest_type_t typ,
                  void *ject, int jown, pdest_object_t * refed)
{
  pdest_object_t     *ob;

  P4EST_ASSERT (dim == 2 || dim == 3);
#ifndef P4EST_BUILD_2D
  P4EST_ASSERT (dim != 2);
#endif
#ifndef P4EST_BUILD_3D
  P4EST_ASSERT (dim != 3);
#endif
  P4EST_ASSERT (typ != PDEST_TYPE_INVALID);
  P4EST_ASSERT (ject != NULL);

  ob = P4EST_ALLOC (pdest_object_t, 1);

  sc_refcount_init (&ob->rc, p4est_package_id);
  ob->dim = dim;
  ob->typ = typ;
  ob->ject = ject;
  ob->jown = jown;

  if (refed != NULL) {
    pdest_object_ref (ob->refed = refed);
  }
  else {
    ob->refed = NULL;
  }

  return ob;
}

void
pdest_object_destroy (pdest_object_t * ob)
{
  P4EST_ASSERT (pdest_object_is_valid (ob));
  P4EST_ASSERT (sc_refcount_is_last (&ob->rc));
  P4EST_EXECUTE_ASSERT_TRUE (pdest_object_unref (ob));
}

int
pdest_object_dim (const pdest_object_t * ob)
{
  P4EST_ASSERT (pdest_object_is_valid (ob));

  return ob->dim;
}

int
pdest_quadrant_is_valid (const pdest_quadrant_t * quadrant)
{
  if (quadrant == NULL) {
    return 0;
  }
  if (quadrant->z == -1) {
    return pdest_p4est_quadrant_is_valid (quadrant);
  }
  else {
    return pdest_p8est_quadrant_is_valid (quadrant);
  }
}

int
pdest_quadrant_is_equal (const pdest_quadrant_t * q1,
                         const pdest_quadrant_t * q2)
{
  P4EST_ASSERT (pdest_quadrant_is_valid (q1));
  P4EST_ASSERT (pdest_quadrant_is_valid (q2));

  /* we compare the z coordinate in any case to check 2D against 3D */
  return q1->x == q2->x && q1->y == q2->y && q1->z == q2->z &&
    q1->level == q2->level;
}

int
pdest_quadrant_is_ancestor (const pdest_quadrant_t * q,
                            const pdest_quadrant_t * r)
{
  const int           dq = pdest_quadrant_dim (q);
  const int           dr = pdest_quadrant_dim (r);

  if (dq != dr) {
    return 0;
  }
  if (dq == 2) {
    p4est_quadrant_t    pq, pr;
    pdest_quadrant_to_p4est (&pq, q);
    pdest_quadrant_to_p4est (&pr, r);
    return p4est_quadrant_is_ancestor (&pq, &pr);
  }
  else {
    P4EST_ASSERT (dq == 3);
    p8est_quadrant_t    pq, pr;
    pdest_quadrant_to_p8est (&pq, q);
    pdest_quadrant_to_p8est (&pr, r);
    return p8est_quadrant_is_ancestor (&pq, &pr);
  }
}

int
pdest_quadrant_dim (const pdest_quadrant_t * quadrant)
{
  P4EST_ASSERT (pdest_quadrant_is_valid (quadrant));

  return quadrant->z == -1 ? 2 : 3;
}

void
pdest_quadrant_length (const pdest_quadrant_t * quadrant, int *dim,
                       p4est_qcoord_t * len, p4est_qcoord_t * root_len)
{
  P4EST_ASSERT (pdest_quadrant_is_valid (quadrant));

  if (quadrant->z == -1) {
    if (dim != NULL) {
      *dim = P4EST_DIM;
    }
    if (len != NULL) {
      *len = P4EST_QUADRANT_LEN (quadrant->level);
    }
    if (root_len != NULL) {
      *root_len = P4EST_ROOT_LEN;
    }
  }
  else {
    if (dim != NULL) {
      *dim = P8EST_DIM;
    }
    if (len != NULL) {
      *len = P8EST_QUADRANT_LEN (quadrant->level);
    }
    if (root_len != NULL) {
      *root_len = P8EST_ROOT_LEN;
    }
  }
}

int
pdest_qmaxlevel (const pdest_object_t * ob)
{
  P4EST_ASSERT (pdest_object_is_valid (ob));
  P4EST_ASSERT (ob->typ == PDEST_TYPE_FOREST);

  return ob->dim == 2 ? P4EST_QMAXLEVEL : P8EST_QMAXLEVEL;
}

int
pdest_is_valid (pdest_object_t * ob_p4est)
{
  P4EST_ASSERT (pdest_object_is_valid (ob_p4est));
  if (ob_p4est->dim == 2) {
    return p4est_is_valid (pdest_p4est_access (ob_p4est));
  }
  else {
    return p8est_is_valid (pdest_p8est_access (ob_p4est));
  }
}

long
pdest_revision (pdest_object_t * ob_p4est)
{
  P4EST_ASSERT (pdest_object_is_valid (ob_p4est));
  if (ob_p4est->dim == 2) {
    return p4est_revision (pdest_p4est_access (ob_p4est));
  }
  else {
    return p8est_revision (pdest_p8est_access (ob_p4est));
  }
}

sc_MPI_Comm
pdest_mpicomm (pdest_object_t * ob_p4est)
{
  P4EST_ASSERT (pdest_object_is_valid (ob_p4est));
  if (ob_p4est->dim == 2) {
    return pdest_p4est_access (ob_p4est)->mpicomm;
  }
  else {
    return pdest_p8est_access (ob_p4est)->mpicomm;
  }
}

int
pdest_info_is_equal (const pdest_info_t * i1, const pdest_info_t * i2)
{
  P4EST_ASSERT (i1 != NULL);
  P4EST_ASSERT (i2 != NULL);

  return i1->dim == i2->dim && i1->revision == i2->revision &&
    i1->first_local_tree == i2->first_local_tree &&
    i1->last_local_tree == i2->last_local_tree &&
    i1->local_num_quadrants == i2->local_num_quadrants &&
    i1->global_num_quadrants == i2->global_num_quadrants &&
    i1->global_quadrants_offset == i2->global_quadrants_offset &&
    i1->data_size == i2->data_size;
}

void
pdest_info (pdest_object_t * ob_p4est, pdest_info_t * info)
{
  P4EST_ASSERT (pdest_object_is_valid (ob_p4est));
  P4EST_ASSERT (info != NULL);

  if ((info->dim = ob_p4est->dim) == 2) {
    p4est_t            *p4est = pdest_p4est_access (ob_p4est);
    info->revision = p4est->revision;
    info->first_local_tree = p4est->first_local_tree;
    info->last_local_tree = p4est->last_local_tree;
    info->local_num_quadrants = p4est->local_num_quadrants;
    info->global_num_quadrants = p4est->global_num_quadrants;
    info->global_quadrants_offset =
      p4est->global_first_quadrant[p4est->mpirank];
    info->data_size = p4est->data_size;
  }
  else {
    p8est_t            *p8est = pdest_p8est_access (ob_p4est);
    info->revision = p8est->revision;
    info->first_local_tree = p8est->first_local_tree;
    info->last_local_tree = p8est->last_local_tree;
    info->local_num_quadrants = p8est->local_num_quadrants;
    info->global_num_quadrants = p8est->global_num_quadrants;
    info->global_quadrants_offset =
      p8est->global_first_quadrant[p8est->mpirank];
    info->data_size = p8est->data_size;
  }
}

void
pdest_tree_info (pdest_object_t * ob_p4est, p4est_topidx_t which_tree,
                 pdest_tree_info_t * tree_info)
{
  P4EST_ASSERT (pdest_object_is_valid (ob_p4est));
  P4EST_ASSERT (tree_info != NULL);

  if ((tree_info->dim = ob_p4est->dim) == 2) {
    p4est_t            *p4est = pdest_p4est_access (ob_p4est);
    p4est_tree_t       *tree;

    P4EST_ASSERT (p4est->connectivity != NULL);
    P4EST_ASSERT (0 <= which_tree);
    P4EST_ASSERT (which_tree < p4est->connectivity->num_trees);

    /* this information is always available */
    tree = p4est_tree_array_index (p4est->trees, which_tree);
    tree_info->tree_quadrants_offset = tree->quadrants_offset;
    tree_info->tree_num_quadrants =
      (p4est_locidx_t) tree->quadrants.elem_count;

    /* we allow calling this function on empty trees */
    if (tree_info->tree_num_quadrants == 0) {
      tree_info->quadrants = NULL;
    }
    else {
      P4EST_ASSERT (p4est->first_local_tree <= which_tree);
      P4EST_ASSERT (which_tree <= p4est->last_local_tree);
      tree_info->quadrants = &tree->quadrants;
    }
  }
  else {
    p8est_t            *p8est = pdest_p8est_access (ob_p4est);
    p8est_tree_t       *tree;

    P4EST_ASSERT (p8est->connectivity != NULL);
    P4EST_ASSERT (0 <= which_tree);
    P4EST_ASSERT (which_tree < p8est->connectivity->num_trees);

    /* this information is always available */
    tree = p8est_tree_array_index (p8est->trees, which_tree);
    tree_info->tree_quadrants_offset = tree->quadrants_offset;
    tree_info->tree_num_quadrants =
      (p4est_locidx_t) tree->quadrants.elem_count;

    /* we allow calling this function on empty trees */
    if (tree_info->tree_num_quadrants == 0) {
      tree_info->quadrants = NULL;
    }
    else {
      P4EST_ASSERT (p8est->first_local_tree <= which_tree);
      P4EST_ASSERT (which_tree <= p8est->last_local_tree);
      tree_info->quadrants = &tree->quadrants;
    }
  }
}

void
pdest_quadrant_index (pdest_quadrant_t * quadrant, const
                      pdest_tree_info_t * tree_info, p4est_locidx_t n)
{
  P4EST_ASSERT (quadrant != NULL);
  P4EST_ASSERT (tree_info != NULL);
  P4EST_ASSERT (tree_info->quadrants != NULL);
  P4EST_ASSERT (tree_info->quadrants->elem_count ==
                (size_t) tree_info->tree_num_quadrants);

  P4EST_ASSERT (0 <= n && n < tree_info->tree_num_quadrants);
  if (tree_info->dim == 2) {
    P4EST_ASSERT (tree_info->quadrants->elem_size ==
                  sizeof (p4est_quadrant_t));
    pdest_quadrant_from_p4est
      (quadrant, p4est_quadrant_array_index (tree_info->quadrants, n));
  }
  else {
    P4EST_ASSERT (tree_info->dim == 3);
    P4EST_ASSERT (tree_info->quadrants->elem_size ==
                  sizeof (p8est_quadrant_t));
    pdest_quadrant_from_p8est
      (quadrant, p8est_quadrant_array_index (tree_info->quadrants, n));
  }
}

int
pdest_quadrant_is_contained (pdest_object_t * ob_p4est,
                             p4est_topidx_t which_tree,
                             const pdest_quadrant_t * quadrant, int rank)
{
  P4EST_ASSERT (pdest_object_is_valid (ob_p4est));
  P4EST_ASSERT (pdest_quadrant_is_valid (quadrant));
  if (ob_p4est->dim == 2) {
    p4est_quadrant_t    p4est_quadrant;
    pdest_quadrant_to_p4est (&p4est_quadrant, quadrant);
    return p4est_comm_is_contained (pdest_p4est_access (ob_p4est),
                                    which_tree, &p4est_quadrant, rank);
  }
  else {
    p8est_quadrant_t    p8est_quadrant;
    pdest_quadrant_to_p8est (&p8est_quadrant, quadrant);
    return p8est_comm_is_contained (pdest_p8est_access (ob_p4est),
                                    which_tree, &p8est_quadrant, rank);
  }
}

pdest_object_t     *
pdest_new_conn (sc_MPI_Comm mpicomm, pdest_object_t * ob_conn,
                int min_level, size_t data_size, void *p4est_user_pointer,
                pdest_init_t init_fn, void *user_pointer)
{
  P4EST_ASSERT (pdest_object_is_valid (ob_conn));
  if (ob_conn->dim == 2) {
    return pdest_p4est_new_conn (mpicomm, ob_conn, min_level, data_size,
                                 p4est_user_pointer, init_fn, user_pointer);
  }
  else {
    return pdest_p8est_new_conn (mpicomm, ob_conn, min_level, data_size,
                                 p4est_user_pointer, init_fn, user_pointer);
  }
}

pdest_object_t     *
pdest_new_copy (pdest_object_t * ob_p4est, int copy_data)
{
  P4EST_ASSERT (pdest_object_is_valid (ob_p4est));
  if (ob_p4est->dim == 2) {
    return pdest_p4est_new_copy (ob_p4est, copy_data);
  }
  else {
    return pdest_p8est_new_copy (ob_p4est, copy_data);
  }
}

void
pdest_refine (pdest_object_t * ob_p4est,
              int refine_recursive, int maxlevel,
              pdest_refine_t refine_fn, pdest_init_t init_fn,
              pdest_replace_t replace_fn, void *user_pointer)
{
  P4EST_ASSERT (pdest_object_is_valid (ob_p4est));
  if (ob_p4est->dim == 2) {
    pdest_p4est_refine (ob_p4est, refine_recursive, maxlevel,
                        refine_fn, init_fn, replace_fn, user_pointer);
  }
  else {
    pdest_p8est_refine (ob_p4est, refine_recursive, maxlevel,
                        refine_fn, init_fn, replace_fn, user_pointer);
  }
}

void
pdest_coarsen (pdest_object_t * ob_p4est,
               int coarsen_recursive, int callback_orphans,
               pdest_coarsen_t coarsen_fn, pdest_init_t init_fn,
               pdest_replace_t replace_fn, void *user_pointer)
{
  P4EST_ASSERT (pdest_object_is_valid (ob_p4est));
  if (ob_p4est->dim == 2) {
    pdest_p4est_coarsen (ob_p4est, coarsen_recursive, callback_orphans,
                         coarsen_fn, init_fn, replace_fn, user_pointer);
  }
  else {
    pdest_p8est_coarsen (ob_p4est, coarsen_recursive, callback_orphans,
                         coarsen_fn, init_fn, replace_fn, user_pointer);
  }
}

void
pdest_partition (pdest_object_t * ob_p4est, int partition_for_coarsening,
                 pdest_weight_t weight_fn, void *user_pointer)
{
  P4EST_ASSERT (pdest_object_is_valid (ob_p4est));
  if (ob_p4est->dim == 2) {
    pdest_p4est_partition (ob_p4est, partition_for_coarsening,
                           weight_fn, user_pointer);
  }
  else {
    pdest_p8est_partition (ob_p4est, partition_for_coarsening,
                           weight_fn, user_pointer);
  }
}

void               *
pdest_user_pointer (pdest_object_t * ob_p4est)
{
  P4EST_ASSERT (pdest_object_is_valid (ob_p4est));
  if (ob_p4est->dim == 2) {
    return pdest_p4est_user_pointer (ob_p4est);
  }
  else {
    return pdest_p8est_user_pointer (ob_p4est);
  }
}
