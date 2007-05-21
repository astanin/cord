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
 *
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

#include "rdsolve.h"
#include "pradi.h"

#include <memory>
using std::auto_ptr;

#ifdef ENABLE_ADI_ISO_FIX
#include <stdlib.h>
#endif

AMesh2D*
rd_step_adi(BCSet const& bcs, double dt, const AMesh2D& m1, string const var,
	string const Dvar, string const Rvar)
throw(MeshException) {
	try {
	auto_ptr<AMesh2D> m2(m1.clone());
	// ADI
#ifdef ENABLE_ADI_ISO_FIX
	if (rand()%2) { // choose sequence of directions: 50% x,y - 50% y,x
#endif
		m2.reset(step_peaceman_rachford_adi_x
			(bcs, 0.5*dt, *m2, var, Dvar, Rvar));
		m2.reset(step_peaceman_rachford_adi_y
			(bcs, 0.5*dt, *m2, var, Dvar, Rvar));
#ifdef ENABLE_ADI_ISO_FIX
	} else {
		m2.reset(step_peaceman_rachford_adi_y
			(bcs, 0.5*dt, *m2, var, Dvar, Rvar));
		m2.reset(step_peaceman_rachford_adi_x
			(bcs, 0.5*dt, *m2, var, Dvar, Rvar));
	}
#endif
	return m2.release();
	} catch(MeshException& e) {
		ostringstream ss;
		ss << "rd_step_adi: " << e.what();
		throw MeshException(ss.str());
	}
}


AMesh2D*
reaction_diffusion_step(BCSet const& bcs, double dt,
	AMesh2D const& m1, string const var,
	string const Dvar, string const Rvar,
	MP::rd_solver_t solver)
throw(MeshException) {
	auto_ptr<AMesh2D> m2(0);
	switch(solver) {
	case MP::RDS_EXPLICIT:
		throw MeshException("reaction_diffusion_step: "
			"RDS_EXPLICIT not implemented");
		break;
	case MP::RDS_ADI:
		m2.reset(rd_step_adi(bcs,dt,m1,var,Dvar,Rvar));
		break;
	case MP::RDS_IMPLICIT:
		throw MeshException("reaction_diffusion_step: "
			"RDS_IMPLICIT not implemented");
		break;
	default:
		throw MeshException("reaction_diffusion_step: unknown solver");
		break;
	}
	return m2.release();
}

AMesh2D*
reaction_diffusion_step(BCSet const& bcs, double dt,
	AMesh2D const& m1, string const var,
	double const D, double const R,
	MP::rd_solver_t solver)
throw(MeshException) {
	auto_ptr<AMesh2D> m2(m1.clone());
	m2->remove_function_ifdef("Dvar_tmp");
	m2->remove_function_ifdef("Rvar_tmp");
	m2->add_function("Dvar_tmp",D);
	m2->add_function("Rvar_tmp",R);
	m2.reset(reaction_diffusion_step(bcs,dt,m1,var,
				"Dvar_tmp","Rvar_tmp",solver));
	m2->remove_function_ifdef("Dvar_tmp");
	m2->remove_function_ifdef("Rvar_tmp");
	return m2.release();
}

