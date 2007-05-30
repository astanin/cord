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
 *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef RDSOLVE_H
#define RDSOLVE_H

#include "amesh2d.h"
#include "spmat.h"
#include "boundary.h"
#include "params.h"

/** time integration step for reaction--diffusion equation
 *
 * \DD{var}{t} = \div (Dvar(x,y)\nabla(var)) + Rvar(x,y)
 *
 * @param bcs	set of boundary conditions
 * @param dt	time step
 * @param m1	initial state of the mesh
 * @param var	name of the mesh variable
 * @param Dvar	variable with diffusion coefficient
 * @param Rvar	variable with reaction term coefficient
 * @return	state of the mesh after time step @c dt
 **/
AMesh2D*
reaction_diffusion_step(BCSet const& bcs, double dt,
	AMesh2D const& m1, string const var,
	string const Dvar, string const Rvar, MP::rd_solver_t solver=MP::RDS_ADI)
throw(MeshException);

/** assume diffusion coefficient and reaction term to be constant everywhere */
AMesh2D*
reaction_diffusion_step(BCSet const& bcs, double dt,
	AMesh2D const& m1, string const var,
	double const D, double const R, MP::rd_solver_t solver=MP::RDS_ADI)
throw(MeshException);

/** @brief fill @c A and @c rhs to construct equation for boundary point
 * @c k_boundary, which satisfies boundary condition @c bc and uses given inner
 * point @c k_inner and point distance @c dx
 */
void
build_boundary_point_eq(ASparseMatrix& A, vector<double>& rhs,
	int const k_boundary, int const k_inner,
	BoundaryCondition const& bc, double const dx);

#endif

