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
AMesh2D*
step_euler_explicit(const Params& p,
	double dt, const AMesh2D& m1, string const var)
throw(MeshException);

/** @brief take @c m1 and make one full step of Peaceman-Rachford
 * alternative directions implicit method (both substeps)
 *
 * @param dt	time step
 * @param m1	previous state of the mesh
 * @param var	name of the mesh variable to be used as $\phi$
 * @param p	model parameters
 * @return	the next state of the mesh
 *
 * The method evaluates
 * \DD{\phi}{t} = \div ( \mu \phi \nabla (\phi \Sigma(\phi))) */
AMesh2D*
phi_step_adi(const Params& p, double dt,
	const AMesh2D& m1, string const var)
throw(MeshException);

/** @brief integrate function @f through out the domain */
double
integrate(const AMesh2D& m, string const f="phi");

AMesh2D*
step_level_set(double dt, const AMesh2D& m1)
throw(MeshException);

#endif

