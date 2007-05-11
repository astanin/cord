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

#ifndef NUTRIENT_H
#define NUTRIENT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "params.h"
#include "dmesh.h"

/** @brief take @c m1 and evaluate new distribution of the nutrient
 * @param m1		previous state of the mesh
 * @param epsilon	precision (equation residual) to be achived
 * @return		the next state of the mesh
 *
 * The method is C-style rather than class. It evaluates
 * \nabla^2 c - \alpha c \Phi = 0 */
AMesh2D*
eval_nutrient(const Params& p, const AMesh2D& m1,
	double const epsilon=1e-3, double const dt=0.0)
throw(MeshException);

#endif

