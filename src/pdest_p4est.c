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

#ifndef P4_TO_P8
#include <pdest_p4est.h>
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_extended.h>
#else
#include <pdest_p8est.h>
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#endif

pdest_object_t     *
pdest_p4est_connectivity_new (p4est_connectivity_t * conn, int absorb_conn)
{
  P4EST_ASSERT (p4est_connectivity_is_valid (conn));

  return pdest_object_new (P4EST_DIM, PDEST_TYPE_CONNECTIVITY,
                           (void *) conn, absorb_conn, NULL);
}

p4est_connectivity_t *
pdest_p4est_connectivity_access (pdest_object_t * ob_conn)
{
  P4EST_ASSERT (pdest_object_is_valid (ob_conn));
  P4EST_ASSERT (ob_conn->typ == PDEST_TYPE_CONNECTIVITY);
  P4EST_ASSERT (ob_conn->dim == P4EST_DIM);

  return (p4est_connectivity_t *) ob_conn->ject;
}

void
pdest_quadrant_from_p4est (pdest_quadrant_t * quadrant,
                           const p4est_quadrant_t * p4est_quadrant)
{
  P4EST_ASSERT (quadrant != NULL);
  P4EST_ASSERT (p4est_quadrant != NULL);
  P4EST_ASSERT (p4est_quadrant_is_valid (p4est_quadrant));

  quadrant->x = p4est_quadrant->x;
  quadrant->y = p4est_quadrant->y;
#ifdef P4_TO_P8
  quadrant->z = p4est_quadrant->z;
#else
  quadrant->z = -1;
#endif
  quadrant->level = (int) p4est_quadrant->level;
  quadrant->user_data = p4est_quadrant->p.user_data;
}

void
pdest_quadrant_to_p4est (p4est_quadrant_t * p4est_quadrant,
                         const pdest_quadrant_t * quadrant)
{
  P4EST_ASSERT (p4est_quadrant != NULL);
  P4EST_ASSERT (quadrant != NULL);

  p4est_quadrant->x = quadrant->x;
  p4est_quadrant->y = quadrant->y;
#ifdef P4_TO_P8
  p4est_quadrant->z = quadrant->z;
#else
  P4EST_ASSERT (quadrant->z == -1);
#endif
  p4est_quadrant->level = (int8_t) quadrant->level;
  P4EST_ASSERT (p4est_quadrant_is_valid (p4est_quadrant));
}

int
pdest_p4est_quadrant_is_valid (const pdest_quadrant_t * quadrant)
{
  p4est_quadrant_t    p4est_squadrant, *p4est_quadrant = &p4est_squadrant;

  if (quadrant == NULL
#ifndef P4_TO_P8
      || quadrant->z != -1
#endif
    ) {
    return 0;
  }

  /* duplicate assignment to avoid assertions in pdest_quadrant_to_p4est */
  p4est_quadrant->x = quadrant->x;
  p4est_quadrant->y = quadrant->y;
#ifdef P4_TO_P8
  p4est_quadrant->z = quadrant->z;
#endif
  p4est_quadrant->level = (int8_t) quadrant->level;

  /* defer result to p4est function */
  return p4est_quadrant_is_valid (p4est_quadrant);
}

#ifdef P4EST_ENABLE_DEBUG

static int
pdest_p4est_is_clan (pdest_quadrant_t * parent, pdest_quadrant_t * children[])
{
  int                 i;
  p4est_quadrant_t    p4est_parent;
  p4est_quadrant_t    p4est_children[P4EST_CHILDREN];

  pdest_quadrant_to_p4est (&p4est_parent, parent);
  for (i = 0; i < P4EST_CHILDREN; ++i) {
    pdest_quadrant_to_p4est (&p4est_children[i], children[i]);
  }
  if (!p4est_quadrant_is_parent (&p4est_parent, &p4est_children[0])) {
    return 0;
  }
  if (!p4est_quadrant_is_familyv (p4est_children)) {
    return 0;
  }
  return 1;
}

#endif /* P4EST_ENABLE_DEBUG */

typedef struct pdest_p4est_context
{
  int                 dim;
  void               *user_pointer;
  int                 ifamily;
  pdest_init_t        init_fn;
  pdest_refine_t      refine_fn;
  pdest_coarsen_t     coarsen_fn;
  pdest_replace_t     replace_fn;
  pdest_weight_t      weight_fn;
  pdest_quadrant_t    squadrant, *quadrant;
  pdest_quadrant_t    squadrants[P4EST_CHILDREN], *quadrants[P4EST_CHILDREN];
  pdest_quadrant_t    squadrant2, *quadrants2[P4EST_CHILDREN];
}
pdest_p4est_context_t;

static void
pdest_p4est_init_cb (p4est_t * p4est,
                     p4est_topidx_t which_tree,
                     p4est_quadrant_t * p4est_quadrant)
{
  pdest_quadrant_t   *quadrant;
  pdest_p4est_context_t *p4est_cont;

  p4est_cont = (pdest_p4est_context_t *) p4est->user_pointer;
  P4EST_ASSERT (p4est_cont != NULL);
  P4EST_ASSERT (p4est_cont->dim == P4EST_DIM);
  P4EST_ASSERT (p4est_cont->quadrant != NULL);
  P4EST_ASSERT (p4est_cont->init_fn != NULL);
  P4EST_ASSERT (p4est_cont->weight_fn == NULL);

  if (p4est_cont->refine_fn != NULL) {
    /* refine */
    P4EST_ASSERT (p4est_cont->coarsen_fn == NULL);
    P4EST_ASSERT (0 <= p4est_cont->ifamily &&
                  p4est_cont->ifamily < P4EST_CHILDREN);

    /* we are filling the quadrants of the new family in succession */
    quadrant = p4est_cont->quadrants[p4est_cont->ifamily++];
    pdest_quadrant_from_p4est (quadrant, p4est_quadrant);
    p4est_cont->init_fn (p4est_cont->user_pointer, which_tree, quadrant);
  }
  else {
    /* new or coarsen */
    P4EST_ASSERT (p4est_cont->ifamily == 0);

    /* there is exactly one new quadrant to initialize */
    pdest_quadrant_from_p4est (p4est_cont->quadrant, p4est_quadrant);
    p4est_cont->init_fn (p4est_cont->user_pointer,
                         which_tree, p4est_cont->quadrant);
  }
}

static int
pdest_p4est_refine_cb (p4est_t * p4est, p4est_topidx_t which_tree,
                       p4est_quadrant_t * p4est_quadrant)
{
  pdest_p4est_context_t *p4est_cont;

  p4est_cont = (pdest_p4est_context_t *) p4est->user_pointer;
  P4EST_ASSERT (p4est_cont != NULL);
  P4EST_ASSERT (p4est_cont->dim == P4EST_DIM);
  P4EST_ASSERT (p4est_cont->quadrant != NULL);
  P4EST_ASSERT (p4est_cont->refine_fn != NULL);

  p4est_cont->ifamily = 0;
  pdest_quadrant_from_p4est (p4est_cont->quadrant, p4est_quadrant);
  return p4est_cont->refine_fn (p4est_cont->user_pointer,
                                which_tree, p4est_cont->quadrant);
}

static int
pdest_p4est_coarsen_cb (p4est_t * p4est, p4est_topidx_t which_tree,
                        p4est_quadrant_t * p4est_quadrants[])
{
  int                 i;
  pdest_p4est_context_t *p4est_cont;

  p4est_cont = (pdest_p4est_context_t *) p4est->user_pointer;
  P4EST_ASSERT (p4est_cont != NULL);
  P4EST_ASSERT (p4est_cont->dim == P4EST_DIM);
  P4EST_ASSERT (p4est_cont->coarsen_fn != NULL);

  if (p4est_quadrants[1] == NULL) {
    P4EST_ASSERT (p4est_cont->quadrants2[1] == NULL);
    pdest_quadrant_from_p4est (p4est_cont->quadrants2[0], p4est_quadrants[0]);
    return p4est_cont->coarsen_fn (p4est_cont->user_pointer,
                                   which_tree, p4est_cont->quadrants2);
  }
  else {
    for (i = 0; i < P4EST_CHILDREN; ++i) {
      pdest_quadrant_from_p4est (p4est_cont->quadrants[i],
                                 p4est_quadrants[i]);
    }
    return p4est_cont->coarsen_fn (p4est_cont->user_pointer,
                                   which_tree, p4est_cont->quadrants);
  }
}

static void
pdest_p4est_replace_cb (p4est_t * p4est, p4est_topidx_t which_tree,
                        int num_outgoing, p4est_quadrant_t * outgoing[],
                        int num_incoming, p4est_quadrant_t * incoming[])
{
  int                 i;
  pdest_p4est_context_t *p4est_cont;

  p4est_cont = (pdest_p4est_context_t *) p4est->user_pointer;
  P4EST_ASSERT (p4est_cont != NULL);
  P4EST_ASSERT (p4est_cont->dim == P4EST_DIM);
  P4EST_ASSERT (p4est_cont->replace_fn != NULL);

  if (p4est_cont->refine_fn != NULL) {
    /* refine */
    P4EST_ASSERT (p4est_cont->coarsen_fn == NULL);
    P4EST_ASSERT (num_outgoing == 1);
    P4EST_ASSERT (num_incoming == P4EST_CHILDREN);

    /* if not done already, fill information on all result quadrants */
    if (p4est_cont->init_fn == NULL) {
      P4EST_ASSERT (p4est_cont->ifamily == 0);
      for (i = 0; i < P4EST_CHILDREN; ++i) {
        pdest_quadrant_from_p4est (p4est_cont->quadrants[i], incoming[i]);
      }
    }
    else {
      P4EST_ASSERT (p4est_cont->ifamily == P4EST_CHILDREN);
    }
    P4EST_ASSERT (pdest_p4est_is_clan (p4est_cont->quadrant,
                                       p4est_cont->quadrants));
    p4est_cont->replace_fn (p4est_cont->user_pointer, which_tree,
                            1, &p4est_cont->quadrant,
                            P4EST_CHILDREN, p4est_cont->quadrants);
  }
  else {
    /* coarsen */
    P4EST_ASSERT (p4est_cont->coarsen_fn != NULL);
    P4EST_ASSERT (p4est_cont->ifamily == 0);
    P4EST_ASSERT (num_outgoing == P4EST_CHILDREN);
    P4EST_ASSERT (num_incoming == 1);

    /* if not done already, fill information on single result quadrant */
    if (p4est_cont->init_fn == NULL) {
      pdest_quadrant_from_p4est (p4est_cont->quadrant, incoming[0]);
    }
    P4EST_ASSERT (pdest_p4est_is_clan (p4est_cont->quadrant,
                                       p4est_cont->quadrants));
    p4est_cont->replace_fn (p4est_cont->user_pointer, which_tree,
                            P4EST_CHILDREN, p4est_cont->quadrants,
                            1, &p4est_cont->quadrant);
  }
}

static int
pdest_p4est_weight_cb (p4est_t * p4est, p4est_topidx_t which_tree,
                       p4est_quadrant_t * p4est_quadrant)
{
  pdest_p4est_context_t *p4est_cont;

  p4est_cont = (pdest_p4est_context_t *) p4est->user_pointer;
  P4EST_ASSERT (p4est_cont != NULL);
  P4EST_ASSERT (p4est_cont->dim == P4EST_DIM);
  P4EST_ASSERT (p4est_cont->init_fn == NULL);
  P4EST_ASSERT (p4est_cont->refine_fn == NULL);
  P4EST_ASSERT (p4est_cont->coarsen_fn == NULL);
  P4EST_ASSERT (p4est_cont->replace_fn == NULL);
  P4EST_ASSERT (p4est_cont->weight_fn != NULL);

  pdest_quadrant_from_p4est (p4est_cont->quadrant, p4est_quadrant);
  return p4est_cont->weight_fn (p4est_cont->user_pointer,
                                which_tree, p4est_cont->quadrant);
}

pdest_object_t     *
pdest_p4est_new_p4est (p4est_t * p4est, int absorb_p4est,
                       pdest_object_t * ob_conn)
{
  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_ASSERT (pdest_p4est_connectivity_access (ob_conn) ==
                p4est->connectivity);

  /* wrap the forest inside an object */
  return pdest_object_new (P4EST_DIM, PDEST_TYPE_FOREST,
                           (void *) p4est, absorb_p4est, ob_conn);
}

pdest_object_t     *
pdest_p4est_new_conn (sc_MPI_Comm mpicomm, pdest_object_t * ob_conn,
                      int min_level, size_t data_size,
                      void *p4est_user_pointer,
                      pdest_init_t init_fn, void *user_pointer)
{
  p4est_connectivity_t *conn;
  p4est_t            *p4est;
  p4est_init_t        p4est_init;
  pdest_p4est_context_t p4cont, *p4est_cont;

  /* access and verify dimension */
  conn = pdest_p4est_connectivity_access (ob_conn);

  if (init_fn == NULL) {
    p4est_init = NULL;
    p4est_cont = NULL;
  }
  else {
    p4est_init = pdest_p4est_init_cb;

    /* create context for callback function */
    memset (&p4cont, 0, sizeof (p4cont));
    p4est_cont = &p4cont;
    p4est_cont->dim = P4EST_DIM;
    p4est_cont->init_fn = init_fn;
    p4est_cont->user_pointer = user_pointer;
    p4est_cont->quadrant = &p4est_cont->squadrant;
#ifndef P4_TO_P8
    p4est_cont->quadrant->z = -1;
#endif
  }

  /* defer to original p4est function */
  p4est = p4est_new_ext (mpicomm, conn, 0, min_level, 1,
                         data_size, p4est_init, p4est_cont);
  p4est->user_pointer = p4est_user_pointer;

  /* wrap result inside an object */
  return pdest_object_new (P4EST_DIM, PDEST_TYPE_FOREST,
                           (void *) p4est, 1, ob_conn);
}

pdest_object_t     *
pdest_p4est_new_copy (pdest_object_t * ob_p4est, int copy_data)
{
  p4est_t            *p4est, *copy;
  pdest_object_t     *ob_conn;

  /* access and verify dimension */
  p4est = pdest_p4est_access (ob_p4est);
  ob_conn = ob_p4est->refed;
  P4EST_ASSERT (pdest_p4est_connectivity_access (ob_conn) ==
                p4est->connectivity);

  /* defer to original p4est function */
  copy = p4est_copy_ext (p4est, copy_data, 0);

  /* wrap result inside an object */
  return pdest_object_new (P4EST_DIM, PDEST_TYPE_FOREST,
                           (void *) copy, 1, ob_conn);
}

void
pdest_p4est_refine (pdest_object_t * ob_p4est,
                    int refine_recursive, int maxlevel,
                    pdest_refine_t refine_fn, pdest_init_t init_fn,
                    pdest_replace_t replace_fn, void *user_pointer)
{
  int                 i;
  void               *p4est_user_pointer;
  p4est_t            *p4est;
  p4est_refine_t      p4est_refine;
  p4est_init_t        p4est_init;
  p4est_replace_t     p4est_replace;
  pdest_p4est_context_t p4cont, *p4est_cont;

  /* access and verify dimension */
  p4est = pdest_p4est_access (ob_p4est);

  P4EST_ASSERT (refine_fn != NULL);
  p4est_refine = pdest_p4est_refine_cb;
  p4est_init = init_fn == NULL ? NULL : pdest_p4est_init_cb;
  p4est_replace = replace_fn == NULL ? NULL : pdest_p4est_replace_cb;

  /* create context for callback functions */
  memset (&p4cont, 0, sizeof (p4cont));
  p4est_cont = &p4cont;
  p4est_cont->dim = P4EST_DIM;
  p4est_cont->init_fn = init_fn;
  p4est_cont->refine_fn = refine_fn;
  p4est_cont->replace_fn = replace_fn;
  p4est_cont->user_pointer = user_pointer;
  p4est_cont->quadrant = &p4est_cont->squadrant;
  for (i = 0; i < P4EST_CHILDREN; ++i) {
    p4est_cont->quadrants[i] = &p4est_cont->squadrants[i];
#ifndef P4_TO_P8
    p4est_cont->quadrants[i]->z = -1;
#endif
  }

  /* defer to original p4est function */
  p4est_user_pointer = p4est->user_pointer;
  p4est->user_pointer = p4est_cont;
  p4est_refine_ext (p4est, refine_recursive, maxlevel,
                    p4est_refine, p4est_init, p4est_replace);
  p4est->user_pointer = p4est_user_pointer;
}

void
pdest_p4est_coarsen (pdest_object_t * ob_p4est,
                     int coarsen_recursive, int callback_orphans,
                     pdest_coarsen_t coarsen_fn, pdest_init_t init_fn,
                     pdest_replace_t replace_fn, void *user_pointer)
{
  int                 i;
  void               *p4est_user_pointer;
  p4est_t            *p4est;
  p4est_coarsen_t     p4est_coarsen;
  p4est_init_t        p4est_init;
  p4est_replace_t     p4est_replace;
  pdest_p4est_context_t p4cont, *p4est_cont;

  /* access and verify dimension */
  p4est = pdest_p4est_access (ob_p4est);

  P4EST_ASSERT (coarsen_fn != NULL);
  p4est_coarsen = pdest_p4est_coarsen_cb;
  p4est_init = init_fn == NULL ? NULL : pdest_p4est_init_cb;
  p4est_replace = replace_fn == NULL ? NULL : pdest_p4est_replace_cb;

  /* create context for callback functions */
  memset (&p4cont, 0, sizeof (p4cont));
  p4est_cont = &p4cont;
  p4est_cont->dim = P4EST_DIM;
  p4est_cont->init_fn = init_fn;
  p4est_cont->coarsen_fn = coarsen_fn;
  p4est_cont->replace_fn = replace_fn;
  p4est_cont->user_pointer = user_pointer;
  p4est_cont->quadrant = &p4est_cont->squadrant;
  for (i = 0; i < P4EST_CHILDREN; ++i) {
    p4est_cont->quadrants[i] = &p4est_cont->squadrants[i];
#ifndef P4_TO_P8
    p4est_cont->quadrants[i]->z = -1;
#endif
    p4est_cont->quadrants2[i] = NULL;
  }
  p4est_cont->quadrants2[0] = &p4est_cont->squadrant2;
#ifndef P4_TO_P8
  p4est_cont->quadrants2[0]->z = -1;
#endif

  /* defer to original p4est function */
  p4est_user_pointer = p4est->user_pointer;
  p4est->user_pointer = p4est_cont;
  p4est_coarsen_ext (p4est, coarsen_recursive, callback_orphans,
                     p4est_coarsen, p4est_init, p4est_replace);
  p4est->user_pointer = p4est_user_pointer;
}

void
pdest_p4est_partition (pdest_object_t * ob_p4est,
                       int partition_for_coarsening, pdest_weight_t weight_fn,
                       void *user_pointer)
{
  void               *p4est_user_pointer;
  p4est_t            *p4est;
  p4est_weight_t      p4est_weight;
  pdest_p4est_context_t p4cont, *p4est_cont;

  /* access and verify dimension */
  p4est = pdest_p4est_access (ob_p4est);

  if (weight_fn == NULL) {
    p4est_weight = NULL;
    p4est_cont = NULL;
  }
  else {
    p4est_weight = pdest_p4est_weight_cb;

    /* create context for callback function */
    memset (&p4cont, 0, sizeof (p4cont));
    p4est_cont = &p4cont;
    p4est_cont->dim = P4EST_DIM;
    p4est_cont->weight_fn = weight_fn;
    p4est_cont->user_pointer = user_pointer;
    p4est_cont->quadrant = &p4est_cont->squadrant;

    /* bend user pointer to internal context */
    p4est_user_pointer = p4est->user_pointer;
    p4est->user_pointer = p4est_cont;
  }

  /* defer to original p4est function */
  p4est_partition_ext (p4est, partition_for_coarsening, p4est_weight);
  if (p4est_cont != NULL) {
    p4est->user_pointer = p4est_user_pointer;
  }
}

p4est_t            *
pdest_p4est_access (pdest_object_t * ob_p4est)
{
  p4est_t            *p4est;

  P4EST_ASSERT (pdest_object_is_valid (ob_p4est));
  P4EST_ASSERT (ob_p4est->typ == PDEST_TYPE_FOREST);
  P4EST_ASSERT (ob_p4est->dim == P4EST_DIM);

  p4est = (p4est_t *) ob_p4est->ject;

  P4EST_ASSERT (pdest_p4est_connectivity_access (ob_p4est->refed) ==
                p4est->connectivity);

  return p4est;
}

void               *
pdest_p4est_user_pointer (pdest_object_t * ob_p4est)
{
  return pdest_p4est_access (ob_p4est)->user_pointer;
}

void
pdest_p4est_ject_destroy (pdest_type_t typ, void *ject)
{
  switch (typ) {
  case PDEST_TYPE_CONNECTIVITY:
    p4est_connectivity_destroy ((p4est_connectivity_t *) ject);
    break;
  case PDEST_TYPE_FOREST:
    p4est_destroy ((p4est_t *) ject);
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}
