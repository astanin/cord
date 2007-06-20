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

#include "singleton.h" 

/// set of command line parameters
class Params {
public:
	typedef enum {
		EQ_POISSON=0,
		EQ_DIFFUSION=1
	} eq_type;
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
	double o2_uptake; ///< rate of oxygen consumption
	double upkeep_per_cell; ///< minimal energy requirement per cell
	double death_rate; ///< rate of cellular death
	BCSet phi_bc; ///< boundary conditions for phi
	BCSet c_bc; ///< boundary conditions for c
	eq_type c_equation; ///< type of equation for c
	bool host_active; ///< host tissue dies and consumes
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
		cell_motility(1e-2),
		o2_uptake(200), upkeep_per_cell(0.15), death_rate(0.8),
		phi_bc(BC::createDirichletBC(this->phi_stress_free),
			BC::createDirichletBC(this->phi_stress_free),
			BC::createNeumannBC(),BC::createNeumannBC()),
		c_bc(BC::createNeumannBC(),BC::createNeumannBC(),
			BC::createDirichletBC(1.0),BC::createNeumannBC()),
		c_equation(EQ_POISSON),
		host_active(false)
		{}
};

struct MethodParams {
	/// solver method for systems of linear equations
	typedef enum {
		SLES_UMFPACK,
		SLES_CG,
		SLES_BICG,
		SLES_BICGSTAB,
		SLES_GMRES
	} sle_solver_t;
	/// solver method for reaction-diffusion equation
	typedef enum {
		RDS_EXPLICIT,
		RDS_ADI,
		RDS_IMPLICIT
	} rd_solver_t;
	/// solver method for Poisson equation
	typedef enum {
		PS_RELAX,
		PS_SLE
	} p_solver_t;
	/// SLE solver
	sle_solver_t sle_solver;
	double sle_solver_accuracy;
	int sle_solver_max_iters;
	int sle_solver_gmres_restart_after;
	/// reaction-diffusion solver
	rd_solver_t rd_solver;
	/// Poisson solver
	p_solver_t p_solver;
	double p_solver_accuracy;
	double p_solver_relax_step;
	int p_solver_relax_max_iters;
	/// level set method
	bool level_set_reset;
	/// default methods
	MethodParams() :
		sle_solver(SLES_CG),
		sle_solver_accuracy(1e-5),
		sle_solver_max_iters(5000),
		sle_solver_gmres_restart_after(1000),
		rd_solver(RDS_ADI),
		p_solver(PS_SLE),
		p_solver_accuracy(1e-5),
		p_solver_relax_step(-1),
		p_solver_relax_max_iters(5000),
		level_set_reset(true)
		{} 
};

typedef MethodParams MP;
typedef Singleton< MethodParams > Method;

#endif

