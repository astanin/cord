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

#ifndef TUMEVAL_H
#define TUMEVAL_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "dmesh.h"
#include "params.h"

/** @brief take @c m1 and make one step of explicit Euler method
 * @param dt	time step
 * @param m1	previous state of the mesh
 * @param var	name of the mesh variable to be used as $\phi$
 * @param p	model parameters
 * @return	the next state of the mesh
 *
 * The method evaluates
 * \DD{\phi}{t} = \div ( \mu \phi \nabla (\phi \Sigma(\phi))) */
template<class fid_t>
AMesh2D<fid_t>*
step_euler_explicit(const Params& p,
	double dt, const AMesh2D<fid_t>& m1, fid_t const var)
throw(MeshException);

/** @brief take @c m1 and make one full time step dt
 * The method is chosen based on global setting Method::it().rd_solver
 *
 * @param dt	time step
 * @param m1	previous state of the mesh
 * @param var	name of the mesh variable to be used as $\phi$
 * @param p	model parameters
 * @return	the next state of the mesh
 *
 * The method evaluates
 * \DD{\phi}{t} = \div ( \mu \phi \nabla (\phi \Sigma(\phi))) */
template<class fid_t>
AMesh2D<fid_t>*
phi_step(const Params& p, double dt, const AMesh2D<fid_t>& m1, fid_t const var)
throw(MeshException);

/** @brief integrate function @f through out the domain */
template<class fid_t>
double
integrate(const AMesh2D<fid_t>& m, fid_t const f=PHI);

template<class fid_t>
AMesh2D<fid_t>*
step_level_set(double dt, const AMesh2D<fid_t>& m1)
throw(MeshException);

template<class fid_t>
AMesh2D<fid_t>*
extrapolate_var(const AMesh2D<fid_t>& m, fid_t var, fid_t var_ls, double v);

template<class fid_t>
AMesh2D<fid_t>*
extrapolate_subphases(const AMesh2D<fid_t>& m, fid_t var_ls,
	fid_t var_t1, fid_t var_t2, fid_t var_h);

template<class fid_t>
void
reconstruct_total_density(AMesh2D<fid_t>& m, fid_t var_ls,
	fid_t var, fid_t var_t1, fid_t var_t2, fid_t var_h);

#endif

