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

#ifndef DMESHOPS_H
#define DMESHOPS_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <iomanip>
#include <sstream>
#include "params.h"

using std::string;
using std::vector;
using std::ostringstream;

string timestamp2str(int timestamp);

// forward declaration
class DMesh;

/** @brief construct DMesh satisfying given parameters */
DMesh
build_mesh(const Params& p);

int
hdf2dx(const Params& p);

int
hdf2gp(const Params& p);

#endif /* DMESHOPS_H */

