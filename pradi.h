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

 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef PRADI_H
#define PRADI_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/* Peaceman-Rachford Alternative Directions Implicit method
for quasilinear parabolic equations */

#include "spmat.h"
#include "function.h"
#include "dmesh.h"
#include "params.h"

/** @brief the first half-step of Peaceman-Rachford alternative directions
 * implicit method
 *
 * @param dt	time step (half of the full time step)
 * @param m1	previous state of the mesh
 * @param var	name of the variable in the mesh
 * @param D_coef_var mesh variable with the diffusion coefficient values
 * @param reaction_term_var mesh variable with the reaction term values
 *
 * This is the first half-step in resolving 2D reaction-diffusion equation
 * 
 * \DD{var}{t} = \div (D_coev_var(x,y) \nable(var)) + reaction_term_var(x,y)
 */
AMesh2D*
step_peaceman_rachford_adi_x(const BCSet& bcs,
	double dt, const AMesh2D& m1, string const var,
	string const D_coef_var, string const reaction_term_var)
throw(MeshException);

/** @brief the second half-step of Peaceman-Rachford alternative directions
 * implicit method
 *
 * @param dt	time step (half of the full time step)
 * @param m1	previous state of the mesh
 * @param var	name of the variable in the mesh
 * @param D_coef_var mesh variable with the diffusion coefficient values
 * @param reaction_term_var mesh variable with the reaction term values
 *
 * This is the second half-step in resolving 2D reaction-diffusion equation
 * 
 * \DD{var}{t} = \div (D_coev_var(x,y) \nable(var)) + reaction_term_var(x,y)
 */
AMesh2D*
step_peaceman_rachford_adi_y(const BCSet& bcs,
	double dt, const AMesh2D& m1, string const var,
	string const D_coef_var, string const reaction_term_var)
throw(MeshException);

#endif

