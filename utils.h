/* Tumour cord growth model */

/*
 * Copyright (C) 2005-2006 Sergey Astanin, Luigi Preziosi
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
#ifndef UTILS_H
#define UTILS_H

#include "amesh2d.h"
#include "boundary.h"
#include "params.h"
#include <string>

using std::string;

string
dbg_stamp(double t);

string
timestamp2str(int timestamp);

BoundaryCondition
create_boundary_condition(const string bc_type, const double bc_val);

void
read_conf_file(Params& p, const string filename);

int
init_params(Params& p, int argc, const char *argv[]);

/** @brief dump mesh stat to gnuplot file */
void
dump2gp(AMesh2D& m, string const& filename);

/** @brief dump mesh stat to gnuplot file */
void
dump2gp(AMesh2D& m, int const timestamp);

/** @brief dump mesh stat to OpenDX file */
void
dump2dx(AMesh2D const& m, string const& filename);

/** @brief dump mesh stat to OpenDX file */
void
dump2dx(AMesh2D const& m, int const timestamp);

void
dump_mesh(AMesh2D& m, int timestamp);

double get_y_size(const AMesh2D& m);
double get_x_size(const AMesh2D& m);

#endif

