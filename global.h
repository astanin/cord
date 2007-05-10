/* Tumour cord growth model */

/* 
 * Copyright (C) 2005-2006 Sergey Astanin, Luigi Preziosi, Andrea Tosin
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

#ifndef GLOBAL_H
#define GLOBAL_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

extern int verbose; ///< print verbose output if @c verbose > 0

extern int show_version; ///< print version info

extern int use_euler_explicit; ///< use Euler explicit method instead of ADI

typedef enum {
	SOLVER_UMFPACK = 0,
	SOLVER_ITERATIVE_EXPLICIT = 1,
	SOLVER_ITERATIVE_IMPLICIT = 2,
	SOLVER_CG = 3,
	SOLVER_BICG = 4,
	SOLVER_BICGSTAB =5,
	SOLVER_GMRES = 6
} poisson_solver_method;

/// method used to solve Poisson equation
extern poisson_solver_method poisson_solver;

/// time step for iterative Poisson solvers (automatic if less than zero)
extern double poisson_solver_iteration_step;

/// accuracy of the convergence condition of Poisson equation solver
extern double poisson_solver_accuracy;

/// maximum number of iterations in iterative method
extern int poisson_solver_max_iterations; 

/// number of GMRES iterations before method restart
extern int gmres_restart_after;

#endif
