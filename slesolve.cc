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

#include "slesolve.h"

#include <stdlib.h>

#include <memory>
using std::auto_ptr;
#include <vector>
using std::vector;

ASparseMatrix*
build_sle_solver_matrix(int const dim, vector<double> x0, double accuracy)
throw(MeshException) {
	if (accuracy < 0) {
		accuracy=Method::it().sle_solver_accuracy;
	}
	auto_ptr<ASparseMatrix> pA;
#ifdef HAVE_LIBUMFPACK
	if (Method::it().sle_solver == MP::SLES_UMFPACK) {
		pA.reset(new UMFPACKMatrix(dim,dim));
	} else
#endif
	if (Method::it().sle_solver == MP::SLES_CG) {
		pA.reset(new LSolverMatrix(dim, x0,
			accuracy,
			Method::it().sle_solver_max_iters,
			LSolverMatrix::CG));
	} else if (Method::it().sle_solver == MP::SLES_BICG) {
		pA.reset(new LSolverMatrix(dim, x0,
			accuracy,
			Method::it().sle_solver_max_iters,
			LSolverMatrix::BICG));
	} else if (Method::it().sle_solver == MP::SLES_BICGSTAB) {
		pA.reset(new LSolverMatrix(dim, x0,
			accuracy,
			Method::it().sle_solver_max_iters,
			LSolverMatrix::BICGstab));
	} else if (Method::it().sle_solver == MP::SLES_GMRES) {
		pA.reset(new LSolverMatrix(dim, x0,
			accuracy,
			Method::it().sle_solver_max_iters,
			LSolverMatrix::GMRES,
			Method::it().sle_solver_gmres_restart_after));
	} else {
		throw MeshException("unknown Method::it().sle_solver");
	}
	return pA.release();
}

template <typename T>
const T& max(const T& x, const T& y)
{
    return x > y ? x : y;
}

template<class fid_t>
void
update_dirichlet_points(AMesh2D<fid_t>& m, const BCSet& bcs, fid_t const var) {
	int xdim=m.get_xdim();
	int ydim=m.get_ydim();
	for (int i=1; i<(xdim-1); ++i) {
		BoundaryCondition bc;
		// north
		bc=bcs.get_north();
		if (bc.get_type() == BoundaryCondition::DIRICHLET_BC) {
			m.set(var,i,ydim-1,bc.c()/bc.a());
		}
		// south
		bc=bcs.get_south();
		if (bc.get_type() == BoundaryCondition::DIRICHLET_BC) {
			m.set(var,i,0,bc.c()/bc.a());
		}
	}
	for (int j=1; j<(ydim-1); ++j) {
		BoundaryCondition bc;
		// west
		bc=bcs.get_west();
		if (bc.get_type() == BoundaryCondition::DIRICHLET_BC) {
			m.set(var,0,j,bc.c()/bc.a());
		}
		// east
		bc=bcs.get_east();
		if (bc.get_type() == BoundaryCondition::DIRICHLET_BC) {
			m.set(var,xdim-1,j,bc.c()/bc.a());
		}
	}
	// corners:
	m.set(var,0,0,max(m[var](1,0),m[var](0,1)));
	m.set(var,xdim-1,0,max(m[var](xdim-2,0),m[var](xdim-1,1)));
	m.set(var,0,ydim-1,max(m[var](1,ydim-1),m[var](0,ydim-2)));
	m.set(var,xdim-1,ydim-1,max(m[var](xdim-2,ydim-1)
						,m[var](xdim-1,ydim-2)));
}

// templates
template
void
update_dirichlet_points<int>(AMesh2D<int>& m, const BCSet& bcs, int const var);

