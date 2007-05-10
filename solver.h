/* Tumour cord growth model */

/*
 * Copyright (C) 2006 Sergey Astanin, Luigi Preziosi, Andrea Tosin
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

#ifndef SOLVER_H
#define SOLVER_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "dmesh.h"
#include "params.h"
#include "global.h"

/** Resolve the model with given set of parameters @c p and given initial
 * condition @c initial. The time to be evaluated is the model parameter.
 * See _Params.
 * TODO: to be implemented as ConcreteSolver class of more general Solver
 * super class, this way alternative models and/or methods may be run for
 * the same initial conditions.
 * @return	final state of the system */
AMesh2D*
solve(const Params& p, const AMesh2D& initial);

#endif

