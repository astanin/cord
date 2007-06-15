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

 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <stdlib.h>

#include "tumeval.h"
#include "global.h"
#include <iostream>
#include <gsl/gsl_integration.h>
#include "gfm.h"
#include "rdsolve.h"

using std::cerr;
using std::string;
using std::ostringstream;
using std::setprecision;
using std::ios;
using std::scientific;

#include <memory>
using std::auto_ptr;

string dbg_stamp(double t);

blitz::Array<double,1>
grad(const AMesh2D& m, string const fid, const int i, const int j,
	const int i1, const int j1, const int i2, const int j2)
throw(MeshException) {
	double df1=m.get(fid,i,j)-m.get(fid,i1,j1);
	double df2=m.get(fid,i,j)-m.get(fid,i2,j2);
	double dx1=m.x(i,j)-m.x(i1,j1);
	double dy1=m.y(i,j)-m.y(i1,j1);
	double dx2=m.x(i,j)-m.x(i2,j2);
	double dy2=m.y(i,j)-m.y(i2,j2);
	double denom=dx1*dy2-dx2*dy1;
	if (denom == 0) {
		ostringstream sstr;
		sstr << "grad: degenerated mesh, cannot eval gradient: ";
		sstr << "(" << m.x(i,j) << "," << m.y(i,j) << "), ";
		sstr << "(" << m.x(i1,j1) << "," << m.y(i1,j1) << "), ";
		sstr << "(" << m.x(i2,j2) << "," << m.y(i2,j2) << "), ";
		sstr << "i=" << i << ", j=" << j << ", ";
		sstr << "i1=" << i1 << ", j1=" << j1 << ", ";
		sstr << "i2=" << i2 << ", j2=" << j2 ;
		string errmsg=sstr.str();
		throw MeshException(errmsg);
	}
	blitz::Array<double,1> g(2);
	g(0)= -(df1*dy2-df2*dy1)/denom;
	g(1)= -(df2*dx1-df1*dx2)/denom;
	return g;
}

blitz::Array<double,1>
smart_grad(const AMesh2D& m, string const fid, const int i, const int j)
throw(MeshException) {
	try {
	if (!m.is_inner(i,j)) {
		ostringstream ss;
		ss << "smart_grad: outer point (" << i << "," << j << ")";
		throw MeshException(ss.str());
	}
	// differentiation points
	int i1=-1, i2=-1, j1=-1, j2=-1;
	// WARNING: assuming uniform rectangular grid
	if (i==0) {
		i1=0;
	} else {
		i1=i-1;
	}
	if (i==(m.get_xdim()-1)) {
		i2=i;
	} else {
		i2=i+1;
	}
	if (j==0) {
		j1=0;
	} else {
		j1=j-1;
	}
	if (j==(m.get_ydim()-1)) {
		j2=j;
	} else {
		j2=j+1;
	}
	// evaluating gradient
	blitz::Array<double,1> g(2);
	g(0)=(m.get(fid,i2,j)-m.get(fid,i1,j))/(m.x(i2,j)-m.x(i1,j));
	g(1)=(m.get(fid,i,j2)-m.get(fid,i,j1))/(m.y(i,j2)-m.y(i,j1));
	return g;
	} catch (MeshException& e) {
		ostringstream ss;
		ss << "smart_grad: " << e.what();
		throw MeshException(ss.str());
	}
}

// {vx,vy} -- velocity field of the material (tissue)
void
eval_v(AMesh2D& m, string const vx="vx", string const vy="vy")
throw(MeshException) {
	m.add_function_ifndef("term_phi_sigma");
	m.add_function_ifndef(vx);
	m.add_function_ifndef(vy);
	TumourSigma tumour(m);
	HostSigma host(m);
	for (int i=0; i<m.get_xdim(); ++i) {
		for (int j=0; j<m.get_ydim(); ++j) {
			double phi=m.get("phi",i,j);
			double psi=m.get("psi",i,j);
			auto_ptr<ADoubleFunction> sigma;
			if (psi>0) {
				sigma.reset(tumour.build_sigma());
			} else {
				sigma.reset(host.build_sigma());
			}
			m.set("term_phi_sigma",i,j,phi*sigma->eval(phi));
		}
	}
	for (int i=0; i<m.get_xdim(); ++i) {
		for (int j=0; j<m.get_ydim(); ++j) {
			blitz::Array<double,1> g=
				smart_grad(m,"term_phi_sigma",i,j);
			double mu=m.get_attr("cell_motility");
			double vxval=(i>0)?(-g(0)*mu):0;
			double vyval=(j>0)?(-g(1)*mu):0;
			m.set(vx,i,j,vxval);
			m.set(vy,i,j,vyval);
		}
	}
	m.remove_function_ifdef("term_phi_sigma");
}

struct _K_derivative_params {
	const AMesh2D *pm;
	ADoubleFunction* sigma;
	ADoubleFunction* sigma_prime;
	_K_derivative_params() :
		pm((AMesh2D*)0), sigma(0), sigma_prime(0) {}
	~_K_derivative_params() {
		if (sigma) {
			delete sigma;
		}
		if (sigma_prime) {
			delete sigma_prime;
		}
	}
};

double
integrate(const AMesh2D& m, string const f) {
	double integral=0.0;
	for (int i=0; i<m.get_xdim(); ++i) {
		for (int j=0; j<m.get_ydim(); ++j) {
			integral+=m.get(f,i,j)*m.area(i,j);
		}
	}
	return integral;
}

double
estimate_optimal_dt(const AMesh2D& m) {
	double max_dt=0.0;
	double max_D=0.0;
	double min_h2=1.0e99;
	double mu=m.get_attr("cell_motility");
	// an estimation for effective diffusion coefficient
	// as $\mu ( \Phi^2\Sigma\prime(\Phi)+\Phi\Sigma(\Phi) )$
	for (int i=0; i<m.get_xdim(); ++i) {
		for (int j=0; j<m.get_ydim(); ++j) {
			double phi=m.get("phi",i,j);
			double psi=m.get("psi",i,j);
			auto_ptr<ADoubleFunction> ps;
			auto_ptr<ADoubleFunction> psp;
			if (psi > 0) {
				ps.reset(TumourSigma(m).build_sigma());
				psp.reset(TumourSigma(m).build_sigma_prime());
			} else {
				ps.reset(HostSigma(m).build_sigma());
				psp.reset(HostSigma(m).build_sigma_prime());
			}
			double s=ps->eval(phi);
			double sp=psp->eval(phi);
			double eff_D=fabs((phi*phi*sp+phi*s)*mu*1.0);
			if (eff_D > max_D) {
				max_D=eff_D;
			}
		}
	}
	// an estimation for min(h^2) for uniform grid
	double dx2=m.get_dx()*m.get_dx();
	double dy2=m.get_dy()*m.get_dy();
	min_h2=(dx2>dy2)?dy2:dx2;
	max_dt=0.005*min_h2/max_D; // 0.1*min(h**2)/max(effective D)
	if (verbose > 1) { // extra verbose
		cerr << dbg_stamp(m.get_time()) << "estimate_optimal_dt: "
			<< "min_h2=" << min_h2 << " max_D=" << max_D << "\n";
	}
	return max_dt;
}

// return maximum value of sqrt(vx^2+vy^2) on the mesh
double v_max(const AMesh2D& m) throw(MeshException) {
	if (!m.defined("vx") || !m.defined("vy")) {
		throw MeshException("v_max: vx or vy not defined");
	}
	double vmax=0.0;
	for (int i=0; i < (m.get_xdim()); ++i) {
		for (int j=0; j < (m.get_ydim()); ++j) {
			double vx, vy, v;
			vx=m.get("vx",i,j);
			vy=m.get("vy",i,j);
			v=sqrt(vx*vx+vy*vy);
			if (v > vmax) {
				vmax=v;
			}
		}
	}
	return vmax;
}

/** Pseudo diffusion coefficient for phi */
double
phi_pseudo_diff_coef(const AMesh2D& m, const int i, const int j,
		string const var) {
	auto_ptr<ADoubleFunction> sigma;
	auto_ptr<ADoubleFunction> sigma_prime;
	double psi=m.get("psi",i,j);
	if (psi > 0) {
		sigma.reset(TumourSigma(m).build_sigma());
		sigma_prime.reset(TumourSigma(m).build_sigma_prime());
	} else {
		sigma.reset(HostSigma(m).build_sigma());
		sigma_prime.reset(HostSigma(m).build_sigma_prime());
	}
	double mu=m.get_attr("cell_motility");
	double phi=m.get(var,i,j);
	double pseudoD=mu*(phi*sigma->eval(phi)+phi*phi*sigma_prime->eval(phi));
	return pseudoD;
}

void
phi_init_pseudo_D(AMesh2D& m, string const var, string const Dvar) {
	// init diffusion coefficient values
	m.remove_function_ifdef(Dvar);
	m.add_function_ifndef(Dvar,0.0);
	for (int i=0; i < (m.get_xdim()); ++i) {
		for (int j=0; j < (m.get_ydim()); ++j) {
			m.set(Dvar,i,j,phi_pseudo_diff_coef(m,i,j,var));
		}
	}
	// debug output
	if (verbose > 1) {
		cerr << dbg_stamp(m.get_time())
			<<setprecision(3)<<setiosflags(ios::scientific)
			<< "max(" << Dvar << ")= " << max(m[Dvar]) << " "
			<< "min(" << Dvar << ")= " << min(m[Dvar]) << "\n";
	}
}

AMesh2D*
phi_step(const Params& p, double dt, const AMesh2D& m1, string const var)
throw(MeshException) {
	try {
	auto_ptr<AMesh2D> m2(m1.clone());
	phi_init_pseudo_D(*m2,var,"pseudo_D");
	switch (Method::it().rd_solver) {
	case MethodParams::RDS_EXPLICIT:
		m2.reset(reaction_diffusion_step(p.phi_bc,dt,*m2,var,
			"pseudo_D","",MP::RDS_EXPLICIT));
		break;
	case MethodParams::RDS_ADI:
		m2.reset(reaction_diffusion_step(p.phi_bc,dt,*m2,var,
			"pseudo_D","",MP::RDS_ADI));
		break;
	case MethodParams::RDS_IMPLICIT:
		m2.reset(reaction_diffusion_step(p.phi_bc,dt,*m2,var,
			"pseudo_D","",MP::RDS_IMPLICIT));
		break;
	default:
		throw MeshException("phi_step: method unknown");
		break;
	}
	m2->remove_function_ifdef("pseudo_D");
	return m2.release();
	} catch(MeshException& e) {
		ostringstream ss;
		ss << "phi_step_adi: " << e.what();
		throw MeshException(ss.str());
	}
}

void
substep_level_set(double dt, AMesh2D& m)
throw(MeshException) {
	try {
	m.remove_function_ifdef("newpsi");
	m.add_function_ifndef("newpsi",0.0);
	// resolve for inner points
	// WARNING: will work for uniform rectangular grid only
	for (int i=1; i < (m.get_xdim()-1); ++i) {
		for (int j=1; j < (m.get_ydim()-1); ++j) {
			int i1, j1, i2, j2;
			if (m.get("vx",i,j) < 0) { // TODO: check > or < here
				i1=i+1;
				j1=j;
			} else {
				i1=i-1;
				j1=j;
			}
			if (m.get("vy",i,j) < 0) { // TODO: check > or < here
				i2=i;
				j2=j+1;
			} else {
				i2=i;
				j2=j-1;
			}
			blitz::Array<double,1> gpsi=
					grad(m,"psi",i,j,i1,j1,i2,j2);
			double psi, vx, vy;
			psi=m.get("psi",i,j);
			vx=m.get("vx",i,j);
			vy=m.get("vy",i,j);
			double vdotgradpsi;
			vdotgradpsi=vx*gpsi(0)+vy*gpsi(1);
			m.set("newpsi",i,j,psi+dt*vdotgradpsi);
		}
	}
	// apply boundary conditions
	// WARNING: valid for rectangular grid only
	for (int j=1; j<(m.get_ydim()-1); ++j) {
		double psi_in;
		psi_in=m.get("newpsi",1,j);
		m.set("newpsi",0,j,psi_in);
		psi_in=m.get("newpsi",m.get_xdim()-2,j);
		m.set("newpsi",m.get_xdim()-1,j,psi_in);
	}
	for (int i=0; i<m.get_xdim(); ++i) {
		double psi_in;
		psi_in=m.get("newpsi",i,1);
		m.set("newpsi",i,0,psi_in);
		psi_in=m.get("newpsi",i,m.get_ydim()-2);
		m.set("newpsi",i,m.get_ydim()-1,psi_in);
	}
	// overwrite old psi
	for (int i=0; i < (m.get_xdim()); ++i) {
		for (int j=0; j < (m.get_ydim()); ++j) {
			double newpsi;
			newpsi=m.get("newpsi",i,j);
			m.set("psi",i,j,newpsi);
		}
	}
	m.remove_function_ifdef("newpsi");
	return;
	} catch (MeshException& e) {
		ostringstream ss;
		ss << "substep_level_set: " << e.what();
		throw MeshException(ss.str());
	}
}

double level_set_function_reset_step(AMesh2D& m,
	string const& var, string const& dvar) {
	using namespace blitz;
	m.remove_function_ifdef(dvar);
	m.add_function_ifndef(dvar,0.0);
	int xdim=m.get_xdim();
	int ydim=m.get_ydim();
	double dx=m.get_dx();
	double dy=m.get_dy();
	for (int i=1; i<(xdim-1); ++i) {
		for (int j=1; j<(ydim-1); ++j) {
			double dpx=m.get(var,i+1,j)-m.get(var,i,j);
			double dmx=m.get(var,i,j)-m.get(var,i-1,j);
			double dpy=m.get(var,i,j+1)-m.get(var,i,j);
			double dmy=m.get(var,i,j)-m.get(var,i,j-1);
			// grad(psi2)
			double gx=0.5*(dpx+fabs(dpx)+dmx-fabs(dmx))/dx;
			double gy=0.5*(dpy+fabs(dpy)+dmy-fabs(dmy))/dy;
			// d(psi2)/dt
			double dpsi2=m.get("Spsi",i,j)*(1-sqrt(gx*gx+gy*gy));
			m.set(dvar,i,j,dpsi2);
		}
	}
	double maxdvar=max(fabs(m[dvar]));
	double dt=0.2*(dx+dy)/maxdvar;
	for (int i=1; i<(xdim-1); ++i) {
		for (int j=1; j<(ydim-1); ++j) {
			m.set(var,i,j,m.get(var,i,j)+dt*m.get(dvar,i,j));
		}
	}
	// zero flux on all boundaries
	for (int i=1; i<(xdim-1); ++i) {
		m.set(var,i,0,m.get(var,i,1));
		m.set(var,i,ydim-1,m.get(var,i,ydim-2));
	}
	for (int j=1; j<(ydim-1); ++j) {
		m.set(var,0,j,m.get(var,1,j));
		m.set(var,xdim-1,j,m.get(var,xdim-2,j));
	}
	// avoid extremums in corner points
	m.set(var,0,0,0.5*(m.get(var,1,0)+m.get(var,0,1)));
	m.set(var,xdim-1,0,0.5*(m.get(var,xdim-2,0)+m.get(var,xdim-1,1)));
	m.set(var,0,ydim-1,0.5*(m.get(var,1,ydim-1)+m.get(var,0,ydim-2)));
	m.set(var,xdim-1,ydim-1,0.5*(m.get(var,xdim-2,ydim-1)+m.get(var,xdim-1,ydim-2)));
	m.remove_function_ifdef(dvar);
	return maxdvar;
}

/** find steady state solution of
 * d psi2 / dt = S(psi) (1- |grad(psi2)|), S(psi) = psi/\sqrt{psi^2+h_x^2+h_y^2}
 */
int level_set_function_reset(AMesh2D& m) {
	m.remove_function_ifdef("psi2");
	m.add_function_ifndef("psi2");
	m.remove_function_ifdef("Spsi");
	m.add_function_ifndef("Spsi");
	double psi;
	// eps is heuristic parameter
	double eps=16*(m.get_dx()*m.get_dx()+m.get_dy()*m.get_dy());
	int xdim=m.get_xdim();
	int ydim=m.get_ydim();
	for (int i=0; i<xdim; ++i) {
		for (int j=0; j<ydim; ++j) {
			psi=m.get("psi",i,j);
			m.set("Spsi",i,j,psi/sqrt(psi*psi+eps));
		}
	}
	int count=0;
	do {
		eps=level_set_function_reset_step(m,"psi2","dpsi2_dt");
		++count;
	} while (eps < Method::it().sle_solver_accuracy);
	for (int i=0; i<xdim; ++i) {
		for (int j=0; j<ydim; ++j) {
			m.set("psi",i,j,m.get("psi2",i,j));
		}
	}
	m.remove_function_ifdef("psi2");
	m.remove_function_ifdef("Spsi");
	return count;
}

/** time step for dpsi/dt - div(psi*v) = 0
 *  with uniform boundary condition for psi (zero flux);
 *  WARNING: this implementation works for uniform rectangular grid only */
AMesh2D*
step_level_set(double dt, const AMesh2D& m1)
throw(MeshException) {
	try {
	auto_ptr<AMesh2D> m2(m1.clone());
	eval_v(*m2);
	if (!m2->defined("vx") || !m2->defined("vy")) {
		throw MeshException("step_level_set: vx or vy not defined");
	}
	double vmax=v_max(*m2);
	// WARNING: for uniform rectangular grid only
	double dx=m2->get_dx();
	double dy=m2->get_dy();
	double hmin=(dx<dy)?dx:dy;
	double dtmax=hmin/vmax;
	int nsubsteps=(int)floor((dt/dtmax))+1;
	double subdt=dt/nsubsteps;
	if (verbose > 1) { // extra verbose
		cerr << dbg_stamp(m1.get_time())
			<< "step_level_set: "
			<< "dt=" << dt
			<< " subdt=" << subdt << "\n";
	}
	for (int i=0; i<nsubsteps; ++i) {
		substep_level_set(subdt,*m2);
	}
	if (Method::it().level_set_reset) {
		int count;
		count=level_set_function_reset(*m2);
		if (verbose > 1) {
			cerr << dbg_stamp(m1.get_time())
				<< "step_level_set: "
				"reset in " << count
				<< " iters\n";
		}
	}
	return m2.release();
	} catch (MeshException& e) {
		ostringstream ss;
		ss << "step_level_set: " << e.what();
		throw MeshException(ss.str());
	}
}

