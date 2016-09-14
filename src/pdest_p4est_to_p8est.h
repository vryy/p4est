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

#ifndef PDEST_P4EST_TO_P8EST_H
#define PDEST_P4EST_TO_P8EST_H

#include <p4est_to_p8est.h>

/* functions */
#define pdest_p4est_connectivity_new    pdest_p8est_connectivity_new
#define pdest_p4est_connectivity_access pdest_p8est_connectivity_access
#define pdest_p4est_new_p4est           pdest_p8est_new_p8est
#define pdest_p4est_new_conn            pdest_p8est_new_conn
#define pdest_p4est_new_copy            pdest_p8est_new_copy
#define pdest_p4est_refine              pdest_p8est_refine
#define pdest_p4est_coarsen             pdest_p8est_coarsen
#define pdest_p4est_partition           pdest_p8est_partition
#define pdest_p4est_access              pdest_p8est_access
#define pdest_p4est_user_pointer        pdest_p8est_user_pointer
#define pdest_p4est_ject_destroy        pdest_p8est_ject_destroy
#define pdest_p4est_quadrant_is_valid   pdest_p8est_quadrant_is_valid
#define pdest_quadrant_from_p4est       pdest_quadrant_from_p8est
#define pdest_quadrant_to_p4est         pdest_quadrant_to_p8est

#endif /* !PDEST_P4EST_TO_P8EST_H */
