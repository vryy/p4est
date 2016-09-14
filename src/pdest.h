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

#ifndef PDEST_H
#define PDEST_H

#include <sc_refcount.h>
#include <p4est_base.h>

SC_EXTERN_C_BEGIN;

typedef enum pdest_type
{
  PDEST_TYPE_FIRST,
  PDEST_TYPE_CONNECTIVITY = PDEST_TYPE_FIRST,
  PDEST_TYPE_FOREST,
#if 0
  PDEST_TYPE_GHOST,
  PDEST_TYPE_LNODES,
#endif
  PDEST_TYPE_LAST,
  PDEST_TYPE_INVALID = PDEST_TYPE_LAST
}
pdest_type_t;

typedef struct pdest_object
{
  sc_refcount_t       rc;
  int                 dim;
  pdest_type_t        typ;
  void               *ject;
  int                 jown;
  struct pdest_object *refed;
}
pdest_object_t;

/*********************** generic object handling **********************/

int                 pdest_object_is_valid (const pdest_object_t * ob);

int                 pdest_object_is_last (const pdest_object_t * ob);

void                pdest_object_ref (pdest_object_t * ob);

int                 pdest_object_unref (pdest_object_t * ob);

pdest_object_t     *pdest_object_new (int dim, pdest_type_t typ,
                                      void *ject, int jown,
                                      pdest_object_t * refed);

void                pdest_object_destroy (pdest_object_t * ob);

int                 pdest_object_dim (const pdest_object_t * ob);

/************************** quadrant handling *************************/

int                 pdest_qmaxlevel (const pdest_object_t * ob);

typedef struct pdest_quadrant
{
  p4est_qcoord_t      x;
  p4est_qcoord_t      y;
  p4est_qcoord_t      z;
  int                 level;
  void               *user_data;
}
pdest_quadrant_t;

int                 pdest_quadrant_is_valid (const pdest_quadrant_t *
                                             quadrant);

int                 pdest_quadrant_is_equal (const pdest_quadrant_t * q1,
                                             const pdest_quadrant_t * q2);

/** Only returns true for unequal quadrants */
int                 pdest_quadrant_is_ancestor (const pdest_quadrant_t * q,
                                                const pdest_quadrant_t * r);

int                 pdest_quadrant_dim (const pdest_quadrant_t * quadrant);

void                pdest_quadrant_length (const pdest_quadrant_t * quadrant,
                                           int *dim, p4est_qcoord_t *
                                           len, p4est_qcoord_t * root_len);

/*************************** forest handling **************************/

typedef struct pdest_info
{
  int                 dim;
  long                revision;
  p4est_topidx_t      first_local_tree, last_local_tree;
  p4est_locidx_t      local_num_quadrants;
  p4est_gloidx_t      global_num_quadrants;
  p4est_gloidx_t      global_quadrants_offset;
  size_t              data_size;
}
pdest_info_t;

typedef struct pdest_tree_info
{
  int                 dim;
  p4est_locidx_t      tree_quadrants_offset;
  p4est_locidx_t      tree_num_quadrants;
  sc_array_t         *quadrants;        /**< treat this as opaque */
}
pdest_tree_info_t;

typedef void        (*pdest_init_t) (void *user_pointer,
                                     p4est_topidx_t which_tree,
                                     pdest_quadrant_t * quadrant);

typedef int         (*pdest_refine_t) (void *user_pointer,
                                       p4est_topidx_t which_tree,
                                       pdest_quadrant_t * quadrant);

typedef int         (*pdest_coarsen_t) (void *user_pointer,
                                        p4est_topidx_t which_tree,
                                        pdest_quadrant_t * quadrants[]);

typedef int         (*pdest_replace_t) (void *user_pointer,
                                        p4est_topidx_t which_tree,
                                        int num_outgoing,
                                        pdest_quadrant_t * outgoing[],
                                        int num_incoming,
                                        pdest_quadrant_t * incoming[]);

typedef int         (*pdest_weight_t) (void *user_pointer,
                                       p4est_topidx_t which_tree,
                                       pdest_quadrant_t * quadrant);

int                 pdest_is_valid (pdest_object_t * ob_p4est);

long                pdest_revision (pdest_object_t * ob_p4est);

sc_MPI_Comm         pdest_mpicomm (pdest_object_t * ob_p4est);

int                 pdest_info_is_equal (const pdest_info_t * i1,
                                         const pdest_info_t * i2);

void                pdest_info (pdest_object_t * ob_p4est,
                                pdest_info_t * info);

void                pdest_tree_info (pdest_object_t * ob_p4est,
                                     p4est_topidx_t which_tree,
                                     pdest_tree_info_t * tree_info);

/** Access a quadrant relative to a tree. */
void                pdest_quadrant_index (pdest_quadrant_t * quadrant,
                                          const pdest_tree_info_t * tree_info,
                                          p4est_locidx_t n);

/** Test whether a quadrant is fully contained in a ranks' owned regien. */
int                 pdest_quadrant_is_contained (pdest_object_t * ob_p4est,
                                                 p4est_topidx_t which_tree,
                                                 const pdest_quadrant_t
                                                 * quadrant, int rank);

pdest_object_t     *pdest_new_conn (sc_MPI_Comm mpicomm,
                                    pdest_object_t * ob_conn, int min_level,
                                    size_t data_size,
                                    void *p4est_user_pointer,
                                    pdest_init_t init_fn, void *user_pointer);

pdest_object_t     *pdest_new_copy (pdest_object_t * ob_p4est, int copy_data);

void               *pdest_user_pointer (pdest_object_t * ob_p4est);

void                pdest_refine (pdest_object_t * ob_p4est,
                                  int refine_recursive, int maxlevel,
                                  pdest_refine_t refine_fn,
                                  pdest_init_t init_fn,
                                  pdest_replace_t replace_fn,
                                  void *user_pointer);

void                pdest_coarsen (pdest_object_t * ob_p4est,
                                   int coarsen_recursive,
                                   int callback_orphans,
                                   pdest_coarsen_t coarsen_fn,
                                   pdest_init_t init_fn,
                                   pdest_replace_t replace_fn,
                                   void *user_pointer);

void                pdest_partition (pdest_object_t * ob_p4est,
                                     int partition_for_coarsening,
                                     pdest_weight_t weight_fn,
                                     void *user_pointer);

SC_EXTERN_C_END;

#endif /* !PDEST_H */
