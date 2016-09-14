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

#ifndef PDEST_P8EST_H
#define PDEST_P8EST_H

#include <pdest.h>
#include <p8est.h>

SC_EXTERN_C_BEGIN;

pdest_object_t     *pdest_p8est_connectivity_new
  (p8est_connectivity_t * conn, int absorb_conn);

p8est_connectivity_t *pdest_p8est_connectivity_access
  (pdest_object_t * ob_conn);

void                pdest_quadrant_from_p8est (pdest_quadrant_t * quadrant,
                                               const p8est_quadrant_t *
                                               p8est_quadrant);

void                pdest_quadrant_to_p8est (p8est_quadrant_t *
                                             p8est_quadrant,
                                             const pdest_quadrant_t *
                                             quadrant);

int                 pdest_p8est_quadrant_is_valid (const pdest_quadrant_t *
                                                   quadrant);

pdest_object_t     *pdest_p8est_new_p8est (p8est_t * p8est, int absorb_p8est,
                                           pdest_object_t * ob_conn);

pdest_object_t     *pdest_p8est_new_conn (sc_MPI_Comm mpicomm,
                                          pdest_object_t * ob_conn,
                                          int min_level, size_t data_size,
                                          void *p4est_user_pointer,
                                          pdest_init_t init_fn,
                                          void *user_pointer);

pdest_object_t     *pdest_p8est_new_copy (pdest_object_t * ob_p8est,
                                          int copy_data);

void                pdest_p8est_refine (pdest_object_t * ob_p8est,
                                        int refine_recursive,
                                        int maxlevel,
                                        pdest_refine_t refine_fn,
                                        pdest_init_t init_fn,
                                        pdest_replace_t replace_fn,
                                        void *user_pointer);

void                pdest_p8est_coarsen (pdest_object_t * ob_p8est,
                                         int coarsen_recursive,
                                         int callback_orphans,
                                         pdest_coarsen_t coarsen_fn,
                                         pdest_init_t init_fn,
                                         pdest_replace_t replace_fn,
                                         void *user_pointer);

void                pdest_p8est_partition (pdest_object_t * ob_p8est,
                                           int partition_for_coarsening,
                                           pdest_weight_t weight_fn,
                                           void *user_pointer);

p8est_t            *pdest_p8est_access (pdest_object_t * ob_p8est);

void               *pdest_p8est_user_pointer (pdest_object_t * ob_p8est);

void                pdest_p8est_ject_destroy (pdest_type_t typ, void *ject);

SC_EXTERN_C_END;

#endif /* !PDEST_P8EST_H */
