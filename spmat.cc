/*
 * Copyright (C) 2007 Sergey Astanin, Luigi Preziosi, Andrea Tosin
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>

#include <math.h>
#include <limits>
#include "spmat.h"

#ifdef HAVE_LIBLAPACK
// LAPACK's tridiagonal solver
extern "C" {
extern void
dgtsv_(int *Np, int *NRHSp, double *DL, double *D, double *DU, double *B,
	int *LDBp, int *INFOp);
}

int
dgtsv(int N, int NRHS, double *DL, double *D, double *DU, double *B, int ldb) {
	int info;
	dgtsv_ (&N, &NRHS, DL, D, DU, B, &ldb, &info);
	return info;
}
#endif

bool
is_symmetric(ASparseMatrix const& A) {
	if (A.get_rows() != A.get_cols()) {
		return false;
	}
	for (int i=1; i<A.get_rows(); ++i) {
		for (int j=0; j<i; ++j) {
			if (fabs(A.get(i,j)-A.get(j,i)) >
				std::numeric_limits<double>::epsilon()) {
				return false;
			}
		}
	}
	return true;
}


