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
#include "utils.h"
#include <gsl/gsl_odeiv.h>

#include <iostream>
#include <errno.h>
#include <memory>

using std::auto_ptr;

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
	double const theta, double const psi=1.0, double const epsilon=1.0,
	double const host_activity=0.0) {
	double atp=atp_balance(phi,c,theta);
	return epsilon*H(-atp)*(-atp)*H(psi) // tumour death
		+ epsilon*H(-psi)*H(-atp)*(-atp)*host_activity; // host death
}

template<class fid_t>
double net_growth_term
(AMesh2D<fid_t> const& m, int const i, int const j, const fid_t& fid) {
	double phi=m[fid](i,j);
	double c=m[CO2](i,j);
	double psi=m[PSI](i,j);
	double theta=m.get_attr("upkeep_per_cell");
	double epsilon=m.get_attr("death_rate");
	double host_activity=m.get_attr("host_activity");
	return growth_term(phi,c,theta,psi)
		-death_term(phi,c,theta,psi,epsilon,host_activity);
}

/** right hand side function for GSL ode-solver
 * y[0] := phi
 * dydt[0] := \Gamma(phi,c), see the paper for details
 * params[0] := c (oxygen concentration)
 * params[1] := theta (upkeep per cell)
 * params[2] := psi (kind of tissue)
 * params[3] := epsilon (death rate)
 * params[4] := host_activity (1 if host dies, 0 if not) */
int rhs_f(double t, const double y[], double dydt[], void *params) {
	double c=*(((double*)(params))+0);
	double theta=*(((double*)(params))+1);
	double psi=*(((double*)(params))+2);
	double eps=*(((double*)(params))+3);
	double host_activity=*(((double*)(params))+4);
	dydt[0]=growth_term(y[0],c,theta,psi)
		-death_term(y[0],c,theta,psi,eps,host_activity);
	return 0;
}

template<class fid_t>
AMesh2D<fid_t>*
step_growth_death(const double dt, const AMesh2D<fid_t>& m1, const fid_t& fid)
throw(MeshException) {
	if (!m1.defined(fid)) {
		ostringstream ss;
		ss << "step_growth_death: " << id2str(fid) << " not defined";
		throw MeshException(ss.str());
	}
	if (!m1.defined(CO2)) {
		throw MeshException("step_growth_death: c not defined");
	}
	auto_ptr<AMesh2D<fid_t> > m2(m1.clone());
	fid_t gfid=PHI_GROWTH;
	m2->add_function_ifndef(gfid);
	for (int i=0; i<m2->get_xdim(); ++i) {
		for (int j=0; j<m2->get_ydim(); ++j) {
			m2->set(gfid,i,j,net_growth_term(*m2,i,j,fid));
		}
	}
	double params[5];
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
	array2d phi_arr=m1[fid];
	array2d c_arr=m1[CO2];
	array2d psi_arr=m1[PSI];
	for (int i=0; i<m1.get_xdim(); ++i) {
		for (int j=0; j<m1.get_ydim(); ++j) {
			double err=0.0;
			double phi=phi_arr(i,j);
			params[0]=c_arr(i,j);
			params[1]=m1.get_attr("upkeep_per_cell");
			params[2]=psi_arr(i,j);
			params[3]=m1.get_attr("death_rate");
			params[4]=m1.get_attr("host_activity");
			gsl_odeiv_step_reset(step);
			gsl_odeiv_step_apply(step, m1.get_time(), 
					dt, &phi, &err, 0, 0, &sys);
			m2->set(fid,i,j,phi);
		}
	}
	gsl_odeiv_step_free(step);
	return m2.release();
}

// template instantiations
template
AMesh2D<int>*
step_growth_death<int>(const double dt, const AMesh2D<int>& m1, const int& fid);

