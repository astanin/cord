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

#include <stdlib.h>

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

double f_atp_per_cell(double const phi) {
	return (1-phi);
}

// TODO: update nutrient solver if g != c.
double g_atp_per_oxygen(double const c) {
	return c;
}

double atp_balance(double const phi, double const c, double const theta) {
	return phi*f_atp_per_cell(phi)*g_atp_per_oxygen(c)-theta*phi;
}

#include <iostream>
using std::cerr;

double growth_term(double const phi, double const c,
	double const theta, double const psi=1.0, double const gamma=1.0) {
	double atp=atp_balance(phi,c,theta);
	return H(atp)*atp*H(psi);
}

double death_term(double const phi, double const c,
	double const theta, double const psi=1.0, double const epsilon=1.0) {
	double atp=atp_balance(phi,c,theta);
	return epsilon*H(-atp)*(-atp)*H(psi);
}

#if 0
double
growth_term(AMesh2D const& m, int const i, int const j, const string& fid) {
	double phi=m.get(fid,i,j);
	double c=m.get("c",i,j);
	double psi=m.get("psi",i,j);
	double theta=m.get_attr("upkeep_per_cell");
	return growth_term(phi,c,theta,psi);
}
#endif

double
net_growth_term(AMesh2D const& m, int const i, int const j, const string& fid) {
	double phi=m.get(fid,i,j);
	double c=m.get("c",i,j);
	double psi=m.get("psi",i,j);
	double theta=m.get_attr("upkeep_per_cell");
	double epsilon=m.get_attr("death_rate");
	return growth_term(phi,c,theta,psi)-death_term(phi,c,theta,psi,epsilon);
}

/** right hand side function for GSL ode-solver
 * y[0] := phi
 * dydt[0] := \Gamma(phi,c), see the paper for details
 * params[0] := c (oxygen concentration)
 * params[1] := theta (upkeep per cell)
 * params[2] := psi (kind of tissue)
 * params[3] := epsilon (death rate) */
int rhs_f(double t, const double y[], double dydt[], void *params) {
	double c=*(((double*)(params))+0);
	double theta=*(((double*)(params))+1);
	double psi=*(((double*)(params))+2);
	double eps=*(((double*)(params))+3);
	dydt[0]=growth_term(y[0],c,theta,psi)-death_term(y[0],c,theta,psi,eps);
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
			params[1]=m1.get_attr("upkeep_per_cell");
			params[2]=m1.get("psi",i,j);
			params[3]=m1.get_attr("death_rate");
			gsl_odeiv_step_reset(step);
			gsl_odeiv_step_apply(step, m1.get_time(), 
					dt, &phi, &err, 0, 0, &sys);
			m2->set(fid,i,j,phi);
		}
	}
	gsl_odeiv_step_free(step);
	return m2.release();
}

