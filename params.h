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

#ifndef PARAMS_H
#define PARAMS_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include "boundary.h"

using std::string;

typedef BoundaryCondition BC;

class BCSet {
private:
	BoundaryCondition condition[4];
public:
	enum {
		NORTH=0,
		EAST=1,
		SOUTH=2,
		WEST=3
	} BOUNDARIES;

	BCSet(BoundaryCondition const north, BoundaryCondition const east,
		BoundaryCondition const south, BoundaryCondition const west) {
		condition[NORTH]=north;
		condition[EAST]=east;
		condition[SOUTH]=south;
		condition[WEST]=west;
	}
	BoundaryCondition get_north() const { return condition[NORTH];}
	void set_north(BoundaryCondition const bc) { condition[NORTH]=bc; }
	BoundaryCondition get_east() const {return condition[EAST];}
	void set_east(BoundaryCondition const bc) { condition[EAST]=bc; }
	BoundaryCondition get_south() const {return condition[SOUTH];}
	void set_south(BoundaryCondition const bc) { condition[SOUTH]=bc; }
	BoundaryCondition get_west() const {return condition[WEST];}
	void set_west(BoundaryCondition const bc) { condition[WEST]=bc; }
	void set_bc(const string side, BoundaryCondition const bc) {
		if (side == "north") {
			set_north(bc);
		} else if (side == "east") {
			set_east(bc);
		} else if (side == "south") {
			set_south(bc);
		} else if (side == "west") {
			set_west(bc);
		}
	}
};

/// set of command line parameters
class Params {
public:
#ifdef HAVE_LIBHDF5
	int hdf2dx; ///< whether convert to OpenDX format or not
	int hdf2gp; ///< whether convert to gnuplot format or not
	string inputfile; ///< filename to read data from
	string outputfile; ///< filename to save the final data
#endif
	double dt; ///< time step to be force (dt < 0 implies automatic step)
	double eval_t; ///< evaluate system behaviour for time @c eval_t
	int xdim; ///< grid resolution along x-axis (grid points)
	int ydim; ///< grid resolution along y-axis (grid points)
	double xsize; ///< length of the physical domain
	double ysize; ///< width of the physical domain
	double dump_every; ///< time between dumping system state to file
	double phi_stress_free; ///< stress free packing density phi0
	double tk1; ///< sigma=tk1*(phi-phi0), for psi>0, phi>phi0
	double ts1; ///< sigma=ts1*(phi-phi0), for psi>0, phi<phi0
	double hk1; ///< sigma=hk1*(phi-phi0), for psi<0, phi>phi0
	double hs1; ///< sigma=hs1*(phi-phi0), for psi<0, phi<phi0
	double initial_cord_length; ///< initial length of cord (psi>0 tissue)
	double initial_cord_width; ///< initial width of cord (psi>0 tissue)
	double initial_cord_x; ///< initial "center" of the tumour seed
	double initial_cord_y; ///< initial "center" of the tumour seed
	double cell_motility; ///< motility of the cellular phase
	double c_critical; ///< critical concentration of nutrient/oxygen
	double consumption_c; ///< density-related nutrient/oxygen consumption
	double gconsumption_c; ///< growth-related nutrient/oxygen consumption
	BCSet phi_bc; ///< boundary conditions for phi
	BCSet c_bc; ///< boundary conditions for c
	/// default parameters
	Params() :
#ifdef HAVE_LIBHDF5
		hdf2dx(0), hdf2gp(0), inputfile(""), outputfile(""),
#endif
		dt(-1.0), eval_t(100.0), xdim(20), ydim(19),
		xsize(2.0), ysize(2.0), dump_every(1.0), phi_stress_free(0.75),
		tk1(1.0), ts1(1.0), hk1(1.0), hs1(1.0),
		initial_cord_length(0.5), initial_cord_width(0.5),
		initial_cord_x(0.0), initial_cord_y(0.0),
		cell_motility(1e-2), c_critical(0.8), consumption_c(6.5),
		gconsumption_c(0.0),
		phi_bc(BC::createDirichletBC(this->phi_stress_free),
			BC::createDirichletBC(this->phi_stress_free),
			BC::createNeumannBC(),BC::createNeumannBC()),
		c_bc(BC::createNeumannBC(),BC::createNeumannBC(),
			BC::createDirichletBC(1.0),BC::createNeumannBC())
		{}
};

#endif

