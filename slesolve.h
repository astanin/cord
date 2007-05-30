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

#include <memory>
using std::auto_ptr;
#include <vector>
using std::vector;

ASparseMatrix*
build_sle_solver_matrix(int const dim, vector<double> x0, double const accuracy)
throw(MeshException) {
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

#endif

