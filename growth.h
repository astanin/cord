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

#ifndef GROWTH_H
#define GROWTH_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "dmesh.h"
#include "utils.h"

inline
double f_atp_per_cell(double const phi) {
	return (1-phi);
}

// TODO: update nutrient solver if g != c.
inline
double g_atp_per_oxygen(double const c) {
	return c;
}

inline
double atp_balance(double const phi, double const c, double const theta) {
	return phi*f_atp_per_cell(phi)*g_atp_per_oxygen(c)-theta*phi;
}

inline
double growth_term(double const phi, double const c,
	double const theta, double const psi=1.0, double const gamma=1.0) {
	double atp=atp_balance(phi,c,theta);
	return H(atp)*atp*H(psi);
}

inline
double death_term(double const phi, double const c,
	double const theta, double const psi=1.0, double const epsilon=1.0,
	double const host_activity=0.0) {
	double atp=atp_balance(phi,c,theta);
	return epsilon*H(-atp)*(-atp)*H(psi) // tumour death
		+ epsilon*H(-psi)*H(-atp)*(-atp)*host_activity; // host death
}

template<class fid_t>
double net_growth_term
(AMesh2D<fid_t> const& m, array2d const& phi, array2d const& c,
	array2d const& psi, int const i, int const j);

/** @brief take @c m1 and evaluate growth-death (ODE)
 * @param dt	time step
 * @param m1	previous state of the mesh
 * @param m2	next state of the mesh (may be the same as m1)
 * @param fid	function ID of phi variable
 * @return	the next state of the mesh
 *
 * The method is C-style rather than class. It evaluates
 * \DD{\Phi}{t} = \Phi*(1-\Phi)*(c - c_{crit}) */
template<class fid_t>
AMesh2D<fid_t>*
step_growth_death(const double dt, const AMesh2D<fid_t>& m1, const fid_t& fid)
throw(MeshException);

/** @brief ODE step in bicomponent tissue */
template<class fid_t>
AMesh2D<fid_t>*
step_bc_tumour_growth_death(const double dt, const AMesh2D<fid_t>& m1,
	const fid_t& fid1, const fid_t& fid2,
	const fid_t& phase, const double where)
throw(MeshException);

#endif
