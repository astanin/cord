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

#include "growth.h"
#include "dmesh.h"
#include <gsl/gsl_odeiv.h>

#include <iostream>
#include <errno.h>
#include <memory>

using std::auto_ptr;

int H(double const x) {
	if (x >= 0.0) {
		return 1;
	} else {
		return 0;
	}
}

#include <iostream>
using std::cerr;

double growth_term(double const phi, double const c,
	double const phi0, double const c0,
	double const psi=1.0, double const gamma=1.0) {
	if ((c-phi*c0) >= 0) {
		return H(psi)*phi*(1-phi);
	} else {
		return 0.0;
	}
}

double death_term(double const phi, double const c,
	double const phi0, double const c0,
	double const psi=1.0, double const gamma=1.0) {
	if ((c-phi*c0) < 0) {
		return H(psi)*phi*(1-phi);
	} else {
		return 0.0;
	}
}

double
growth_term(AMesh2D const& m, int const i, int const j, const string& fid) {
	double phi=m.get(fid,i,j);
	double c=m.get("c",i,j);
	double psi=m.get("psi",i,j);
	double c0=m.get_attr("c_critical");
	double phi0=m.get_attr("phi0");
	return growth_term(phi,c,phi0,c0,psi);
}

double
net_growth_term(AMesh2D const& m, int const i, int const j, const string& fid) {
	double phi=m.get(fid,i,j);
	double c=m.get("c",i,j);
	double psi=m.get("psi",i,j);
	double c0=m.get_attr("c_critical");
	double phi0=m.get_attr("phi0");
	return growth_term(phi,c,phi0,c0,psi)-death_term(phi,c,phi0,c0,psi);
}

/** right hand side function for GSL ode-solver
 * y[0] := phi
 * dydt[0] := gamma*phi*(1-phi)*(c-c_critical), 
 * 	where c_critical is a critical value of c
 * params[0] := gamma (growth rate) (==1.0 for non-dimensional model)
 * params[1] := c   (nutrient concentration)
 * params[2] := c_critical  (critical value of c) */
int rhs_f(double t, const double y[], double dydt[], void *params) {
	double c=*(((double*)(params))+0);
	double phi0=*(((double*)(params))+1);
	double c0=*(((double*)(params))+2);
	double psi=*(((double*)(params))+3);
	dydt[0]=growth_term(y[0],c,phi0,c0,psi)-death_term(y[0],c,phi0,c0,psi);
	return 0;
}

AMesh2D*
step_growth_death(const double dt, const AMesh2D& m1, const string& fid)
throw(MeshException) {
	if (!m1.defined(fid)) {
		ostringstream ss;
		ss << "step_growth_death: " << fid << " not defined";
		throw MeshException(ss.str());
	}
	if (!m1.defined("c")) {
		throw MeshException("step_growth_death: c not defined");
	}
	auto_ptr<AMesh2D> m2(m1.clone());
	string gfid=fid+"_growth";
	m2->add_function_ifndef(gfid);
	for (int i=0; i<m2->get_xdim(); ++i) {
		for (int j=0; j<m2->get_ydim(); ++j) {
			m2->set(gfid,i,j,net_growth_term(*m2,i,j,fid));
		}
	}
	double params[4];
	gsl_odeiv_system sys;
	sys.function=rhs_f;
	sys.jacobian=0;
	sys.dimension=1;
	sys.params=params;
	gsl_odeiv_step *step=0;
	step=gsl_odeiv_step_alloc(gsl_odeiv_step_rk4, 1);
	if (step == (void*)ENOMEM) {
		throw MeshException
			("step_growth_death: cannot alloc gsl_odeiv_step");
	}
	for (int i=0; i<m1.get_xdim(); ++i) {
		for (int j=0; j<m1.get_ydim(); ++j) {
			double err=0.0;
			double phi=m1.get(fid,i,j);
			params[0]=m1.get("c",i,j);
			params[1]=m1.get_attr("phi0");
			params[2]=m1.get_attr("c_critical");
			params[3]=m1.get("psi",i,j);
			gsl_odeiv_step_reset(step);
			gsl_odeiv_step_apply(step, m1.get_time(), 
					dt, &phi, &err, 0, 0, &sys);
			m2->set(fid,i,j,phi);
		}
	}
	gsl_odeiv_step_free(step);
	return m2.release();
}

