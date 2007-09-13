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

#include <iostream>
using std::cerr;

double net_growth_term
(array2d const& phi, array2d const& c, array2d const& psi,
	double const theta, double const epsilon, double const host_activity,
	int const i, int const j) {
	return growth_term(phi(i,j),c(i,j),theta,psi(i,j))
	-death_term(phi(i,j),c(i,j),theta,psi(i,j),epsilon,host_activity);
}

inline
double
pos(double const value) {
	return (value>=0)?value:0.0;
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

/** right hand side function for GSL ode-solver (bicomponent tissue)
 * y[0] := phi1
 * y[1] := phi2
 * params[0] := c_o (oxygen)
 * params[1] := c_g (glucose)
 * params[2] := psi (level set function)
 * params[3] := psi_where (subdomain selector)
 * params[4] := gamma (growth rate)
 * params[5] := epsilon (death rate)
 * params[6] := theta (upkeep per cell)
 * params[7] := k (relative rate of oxidation in anaerobic catabolic pathway)
 * params[8] := nu (rate of conversion from aerobic to anaerobic cells)
 */
int bc_tumour_rhs_f(double t, const double y[], double dydt[], void *params) {
	double phi1=y[0];
	double phi2=y[1];
	double *param=(double*)params;
	double c_o=param[0];
	double c_g=param[1];
	double psi=param[2];
	double where=param[3];
	double gamma=param[4];
	double eps=param[5];
	double theta=param[6];
	double k=param[7];
	double nu=param[8];
	if (psi*where >= 0) {
		double atp_per_aerobic=(1-phi1)*c_o-theta;
		double atp_per_anaerobic=k*(1-phi2)*c_g-theta;
		double switch_rate=nu*phi1*(atp_per_aerobic<0?1.0:0.0);
		dydt[0]=//gamma*phi1*pos(atp_per_aerobic)
			//-eps*phi1*pos(-atp_per_aerobic)
			-switch_rate;
		dydt[1]=//gamma*phi2*pos(atp_per_anaerobic)
			//-eps*phi2*pos(-atp_per_anaerobic)
			+switch_rate;
	} else {
		dydt[0]=0.0;
		dydt[1]=0.0;
	}
	return 0;
}

template<class fid_t>
AMesh2D<fid_t>*
step_bc_tumour_growth_death(const double dt, const AMesh2D<fid_t>& m1,
	const fid_t& fid1, const fid_t& fid2,
	const fid_t& phase, const double where)
throw(MeshException) {
	if (!m1.defined(fid1)) {
		ostringstream ss;
		ss << "step_growth_death: " << id2str(fid1) << " not defined";
		throw MeshException(ss.str());
	}
	if (!m1.defined(fid2)) {
		ostringstream ss;
		ss << "step_growth_death: " << id2str(fid2) << " not defined";
		throw MeshException(ss.str());
	}
	if (!m1.defined(CO2)) {
		throw MeshException("step_growth_death: c_o not defined");
	}
	if (!m1.defined(GLC)) {
		throw MeshException("step_growth_death: c_g not defined");
	}
	auto_ptr<AMesh2D<fid_t> > m2(m1.clone());
	double params[9];
	gsl_odeiv_system sys;
	sys.function=bc_tumour_rhs_f;
	sys.jacobian=0;
	sys.dimension=2;
	sys.params=params;
	gsl_odeiv_step *step=0;
	step=gsl_odeiv_step_alloc(gsl_odeiv_step_rk4, 2);
	if (step == (void*)ENOMEM) {
		throw MeshException("step_bc_tumour_growth_death: "
				"cannot alloc gsl_odeiv_step");
	}
	array2d m1_phi1=m1[fid1];
	array2d m1_phi2=m1[fid2];
	array2d m1_c_o=m1[CO2];
	array2d m1_c_g=m1[GLC];
	array2d m1_psi=m1[PSI];
	double theta=m1.get_attr("upkeep_per_cell");
	double eps=m1.get_attr("death_rate");
	double k=m1.get_attr("anaerobic_rate");
	double nu=m1.get_attr("conversion_rate");
	double gamma=1.0;
	int xdim=m1.get_xdim();
	int ydim=m1.get_ydim();
	// results
	array2d m2_phi1=(*m2)[fid1];
	array2d m2_phi2=(*m2)[fid2];
	// these are the same for all points:
	params[3]=where;
	params[4]=gamma;
	params[5]=eps;
	params[6]=theta;
	params[7]=k;
	params[8]=nu;
	double phi[2];
	// now solve ODE in every point
	for (int i=0; i<xdim; ++i) {
		for (int j=0; j<ydim; ++j) {
			double err=0.0;
			phi[0]=m1_phi1(i,j);
			phi[1]=m1_phi2(i,j);
			params[0]=m1_c_o(i,j);
			params[1]=m1_c_g(i,j);
			params[2]=m1_psi(i,j);
			gsl_odeiv_step_reset(step);
			gsl_odeiv_step_apply(step, m1.get_time(), 
					dt, phi, &err, 0, 0, &sys);
			m2_phi1(i,j)=phi[0];
			m2_phi2(i,j)=phi[1];
		}
	}
	gsl_odeiv_step_free(step);
	return m2.release();
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
	// shortcuts to arrays
	array2d m2_phi=(*m2)[fid];
	array2d m2_c=(*m2)[CO2];
	array2d m2_psi=(*m2)[PSI];
	array2d m2_gphi=(*m2)[gfid];
	double theta=m2->get_attr("upkeep_per_cell");
	double epsilon=m2->get_attr("death_rate");
	double h_a=m2->get_attr("host_activity");
	int xdim=m2->get_xdim();
	int ydim=m2->get_ydim();
	for (int i=0; i<xdim; ++i) {
		for (int j=0; j<ydim; ++j) {
			m2_gphi(i,j)=net_growth_term(m2_phi,m2_c,m2_psi,
						theta,epsilon,h_a,i,j);
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
	array2d m1_phi=m1[fid];
	array2d m1_c=m1[CO2];
	array2d m1_psi=m1[PSI];
	double upkeep_per_cell=m1.get_attr("upkeep_per_cell");
	double death_rate=m1.get_attr("death_rate");
	double host_activity=m1.get_attr("host_activity");
	for (int i=0; i<xdim; ++i) {
		for (int j=0; j<ydim; ++j) {
			double err=0.0;
			double phi=m1_phi(i,j);
			params[0]=m1_c(i,j);
			params[1]=upkeep_per_cell;
			params[2]=m1_psi(i,j);
			params[3]=death_rate;
			params[4]=host_activity;
			gsl_odeiv_step_reset(step);
			gsl_odeiv_step_apply(step, m1.get_time(), 
					dt, &phi, &err, 0, 0, &sys);
			m2_phi(i,j)=phi;
		}
	}
	gsl_odeiv_step_free(step);
	return m2.release();
}

// template instantiations
template
AMesh2D<int>*
step_growth_death<int>(const double dt, const AMesh2D<int>& m1, const int& fid);

template
AMesh2D<int>*
step_bc_tumour_growth_death(const double dt, const AMesh2D<int>& m1,
	const int& fid1, const int& fid2,
	const int& phase, const double where);

