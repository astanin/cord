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

double f_atp_per_cell(double const phi);

double g_atp_per_oxygen(double const c);

template<class fid_t>
double net_growth_term
(AMesh2D<fid_t> const& m, int const i, int const j, const fid_t& fid);

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

#endif
