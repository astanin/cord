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
template<class fid_t>
void
dump2gp(AMesh2D<fid_t>& m, string const& filename);

/** @brief dump mesh stat to gnuplot file */
template<class fid_t>
void
dump2gp(AMesh2D<fid_t>& m, int const timestamp);

/** @brief dump mesh stat to OpenDX file */
template<class fid_t>
void
dump2dx(AMesh2D<fid_t> const& m, string const& filename);

/** @brief dump mesh stat to OpenDX file */
template<class fid_t>
void
dump2dx(AMesh2D<fid_t> const& m, int const timestamp);

template<class fid_t>
void
dump_mesh(AMesh2D<fid_t>& m, int timestamp);

template<class fid_t>
double get_y_size(const AMesh2D<fid_t>& m);
template<class fid_t>
double get_x_size(const AMesh2D<fid_t>& m);

template<class fid_t>
double norm_1(AMesh2D<fid_t> const& m, fid_t const fid);
template<class fid_t>
double norm_2(AMesh2D<fid_t> const& m, fid_t const fid);

/** @brief Heaviside funcion */
int H(double const x);

#endif

