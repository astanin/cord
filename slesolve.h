/* Tumour cord growth model */

/*
 * Copyright (C) 2005-2007 Sergey Astanin, Luigi Preziosi, Andrea Tosin
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef SLESOLVE_H
#define SLESOLVE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "spmat.h"
#include "amesh2d.h"
#include "params.h"

ASparseMatrix*
build_sle_solver_matrix(int const dim, vector<double> x0, double accuracy=-1.0)
throw(MeshException);

void
update_dirichlet_points(AMesh2D& m, const BCSet& bcs, string const var);

#endif

