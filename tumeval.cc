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

template<class fid_t>
blitz::Array<double,1>
grad(const AMesh2D<fid_t>& m, fid_t const fid, const int i, const int j,
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

template<class fid_t>
blitz::Array<double,1>
smart_grad(const AMesh2D<fid_t>& m, fid_t const fid, const int i, const int j)
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
template<class fid_t>
void
eval_v(AMesh2D<fid_t>& m, fid_t const vx=VX, fid_t const vy=VY)
throw(MeshException) {
	m.add_function_ifndef(TERM_PHI_SIGMA);
	m.add_function_ifndef(vx);
	m.add_function_ifndef(vy);
	TumourSigma<fid_t> tumour(m);
	HostSigma<fid_t> host(m);
	for (int i=0; i<m.get_xdim(); ++i) {
		for (int j=0; j<m.get_ydim(); ++j) {
			double phi=m[PHI](i,j);
			double psi=m[PSI](i,j);
			auto_ptr<ADoubleFunction> sigma;
			if (psi>0) {
				sigma.reset(tumour.build_sigma());
			} else {
				sigma.reset(host.build_sigma());
			}
			m.set(TERM_PHI_SIGMA,i,j,phi*sigma->eval(phi));
		}
	}
	for (int i=0; i<m.get_xdim(); ++i) {
		for (int j=0; j<m.get_ydim(); ++j) {
			blitz::Array<double,1> g=
				smart_grad<int>(m,TERM_PHI_SIGMA,i,j);
			double mu=m.get_attr("cell_motility");
			double vxval=(i>0)?(-g(0)*mu):0;
			double vyval=(j>0)?(-g(1)*mu):0;
			m.set(vx,i,j,vxval);
			m.set(vy,i,j,vyval);
		}
	}
	m.remove_function_ifdef(TERM_PHI_SIGMA);
}

template<class fid_t>
struct _K_derivative_params {
	const AMesh2D<fid_t> *pm;
	ADoubleFunction* sigma;
	ADoubleFunction* sigma_prime;
	_K_derivative_params() :
		pm((AMesh2D<fid_t>*)0), sigma(0), sigma_prime(0) {}
	~_K_derivative_params() {
		if (sigma) {
			delete sigma;
		}
		if (sigma_prime) {
			delete sigma_prime;
		}
	}
};

template<class fid_t>
double
integrate(const AMesh2D<fid_t>& m, fid_t const f) {
	double integral=0.0;
	for (int i=0; i<m.get_xdim(); ++i) {
		for (int j=0; j<m.get_ydim(); ++j) {
			integral+=m[f](i,j)*m.area(i,j);
		}
	}
	return integral;
}

template<class fid_t>
double
estimate_optimal_dt(const AMesh2D<fid_t>& m) {
	double max_dt=0.0;
	double max_D=0.0;
	double min_h2=1.0e99;
	double mu=m.get_attr("cell_motility");
	// an estimation for effective diffusion coefficient
	// as $\mu ( \Phi^2\Sigma\prime(\Phi)+\Phi\Sigma(\Phi) )$
	for (int i=0; i<m.get_xdim(); ++i) {
		for (int j=0; j<m.get_ydim(); ++j) {
			double phi=m[PHI](i,j);
			double psi=m[PSI](i,j);
			auto_ptr<ADoubleFunction> ps;
			auto_ptr<ADoubleFunction> psp;
			if (psi > 0) {
				ps.reset(TumourSigma<fid_t>(m).build_sigma());
				psp.reset(TumourSigma<fid_t>(m).build_sigma_prime());
			} else {
				ps.reset(HostSigma<fid_t>(m).build_sigma());
				psp.reset(HostSigma<fid_t>(m).build_sigma_prime());
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
template<class fid_t>
double v_max(const AMesh2D<fid_t>& m)
throw(MeshException) {
	if (!m.defined(VX) || !m.defined(VY)) {
		throw MeshException("v_max: vx or vy not defined");
	}
	double vmax=0.0;
	for (int i=0; i < (m.get_xdim()); ++i) {
		for (int j=0; j < (m.get_ydim()); ++j) {
			double vx, vy, v;
			vx=m[VX](i,j);
			vy=m[VY](i,j);
			v=sqrt(vx*vx+vy*vy);
			if (v > vmax) {
				vmax=v;
			}
		}
	}
	return vmax;
}

/** Pseudo diffusion coefficient for phi */
template<class fid_t>
double
phi_pseudo_diff_coef
(const AMesh2D<fid_t>& m, const int i, const int j, fid_t const var) {
	auto_ptr<ADoubleFunction> sigma;
	auto_ptr<ADoubleFunction> sigma_prime;
	double psi=m[PSI](i,j);
	if (psi > 0) {
		sigma.reset(TumourSigma<fid_t>(m).build_sigma());
		sigma_prime.reset(TumourSigma<fid_t>(m).build_sigma_prime());
	} else {
		sigma.reset(HostSigma<fid_t>(m).build_sigma());
		sigma_prime.reset(HostSigma<fid_t>(m).build_sigma_prime());
	}
	double mu=m.get_attr("cell_motility");
	double phi=m.get(var,i,j);
	double pseudoD=mu*(phi*sigma->eval(phi)+phi*phi*sigma_prime->eval(phi));
	return pseudoD;
}

template<class fid_t>
void
phi_init_pseudo_D(AMesh2D<fid_t>& m, fid_t const var, fid_t const Dvar) {
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
			<< "max("<< id2str(Dvar) <<")= "<< max(m[Dvar]) << " "
			<< "min("<< id2str(Dvar) <<")= "<< min(m[Dvar]) << "\n";
	}
}

template<class fid_t>
AMesh2D<fid_t>*
phi_step(const Params& p, double dt, const AMesh2D<fid_t>& m1, fid_t const var)
throw(MeshException) {
	try {
	auto_ptr<AMesh2D<fid_t> > m2(m1.clone());
	phi_init_pseudo_D<fid_t>(*m2,var,PSEUDO_D);
	switch (Method::it().rd_solver) {
	case MethodParams::RDS_EXPLICIT:
		m2.reset(reaction_diffusion_step<fid_t>(p.phi_bc,dt,*m2,var,
			PSEUDO_D,NONE,MP::RDS_EXPLICIT));
		break;
	case MethodParams::RDS_ADI:
		m2.reset(reaction_diffusion_step<fid_t>(p.phi_bc,dt,*m2,var,
			PSEUDO_D,NONE,MP::RDS_ADI));
		break;
	case MethodParams::RDS_IMPLICIT:
		m2.reset(reaction_diffusion_step<fid_t>(p.phi_bc,dt,*m2,var,
			PSEUDO_D,NONE,MP::RDS_IMPLICIT));
		break;
	default:
		throw MeshException("phi_step: method unknown");
		break;
	}
	m2->remove_function_ifdef(PSEUDO_D);
	return m2.release();
	} catch(MeshException& e) {
		ostringstream ss;
		ss << "phi_step_adi: " << e.what();
		throw MeshException(ss.str());
	}
}

template<class fid_t>
void
substep_level_set(double dt, AMesh2D<fid_t>& m)
throw(MeshException) {
	try {
	m.remove_function_ifdef(NEWPSI);
	m.add_function_ifndef(NEWPSI,0.0);
	// resolve for inner points
	// WARNING: will work for uniform rectangular grid only
	double dx=m.get_dx();
	double dy=m.get_dy();
	for (int i=1; i < (m.get_xdim()-1); ++i) {
		for (int j=1; j < (m.get_ydim()-1); ++j) {
			double psi, vx, vy;
			psi=m[PSI](i,j);
			vx=m[VX](i,j);
			vy=m[VY](i,j);
			double vdotgradpsi;
			vdotgradpsi=0.5*(vx+fabs(vx))*
				(m[PSI](i,j)-m[PSI](i-1,j))/dx
				+0.5*(vx-fabs(vx))*
				(m[PSI](i+1,j)-m[PSI](i,j))/dx
				+0.5*(vy+fabs(vy))*
				(m[PSI](i,j)-m[PSI](i,j-1))/dy
				+0.5*(vy-fabs(vy))*
				(m[PSI](i,j+1)-m[PSI](i,j))/dy;
			m.set(NEWPSI,i,j,psi-dt*vdotgradpsi);
		}
	}
	// apply boundary conditions
	// WARNING: valid for rectangular grid only
	for (int j=1; j<(m.get_ydim()-1); ++j) {
		double psi_in;
		psi_in=m.get(NEWPSI,1,j);
		m.set(NEWPSI,0,j,psi_in);
		psi_in=m.get(NEWPSI,m.get_xdim()-2,j);
		m.set(NEWPSI,m.get_xdim()-1,j,psi_in);
	}
	for (int i=0; i<m.get_xdim(); ++i) {
		double psi_in;
		psi_in=m.get(NEWPSI,i,1);
		m.set(NEWPSI,i,0,psi_in);
		psi_in=m.get(NEWPSI,i,m.get_ydim()-2);
		m.set(NEWPSI,i,m.get_ydim()-1,psi_in);
	}
	// overwrite old psi
	for (int i=0; i < (m.get_xdim()); ++i) {
		for (int j=0; j < (m.get_ydim()); ++j) {
			double newpsi;
			newpsi=m[NEWPSI](i,j);
			m.set(PSI,i,j,newpsi);
		}
	}
	m.remove_function_ifdef(NEWPSI);
	return;
	} catch (MeshException& e) {
		ostringstream ss;
		ss << "substep_level_set: " << e.what();
		throw MeshException(ss.str());
	}
}

template<class fid_t>
double level_set_function_reset_step(AMesh2D<fid_t>& m,
	fid_t const& var, fid_t const& dvar) {
	using namespace blitz;
	m.remove_function_ifdef(dvar);
	m.add_function_ifndef(dvar,0.0);
	int xdim=m.get_xdim();
	int ydim=m.get_ydim();
	double dx=m.get_dx();
	double dy=m.get_dy();
	array2d var_a=m[var];
	array2d spsi_a=m[SPSI];
	array2d dpsi_a=m[dvar];
	for (int i=1; i<(xdim-1); ++i) {
		for (int j=1; j<(ydim-1); ++j) {
			double dpx=var_a(i+1,j)-var_a(i,j);
			double dmx=var_a(i,j)-var_a(i-1,j);
			double dpy=var_a(i,j+1)-var_a(i,j);
			double dmy=var_a(i,j)-var_a(i,j-1);
			// grad(psi2)
			double gx=0.5*(dpx+fabs(dpx)+dmx-fabs(dmx))/dx;
			double gy=0.5*(dpy+fabs(dpy)+dmy-fabs(dmy))/dy;
			// d(psi2)/dt
			dpsi_a(i,j)=spsi_a(i,j)*(1-sqrt(gx*gx+gy*gy));
		}
	}
	double maxdvar=max(fabs(m[dvar]));
	double dt=0.2*(dx+dy)/maxdvar;
	for (int i=1; i<(xdim-1); ++i) {
		for (int j=1; j<(ydim-1); ++j) {
			var_a(i,j)=var_a(i,j)+dt*dpsi_a(i,j);
		}
	}
	// zero flux on all boundaries
	for (int i=1; i<(xdim-1); ++i) {
		var_a(i,0)=var_a(i,1);
		var_a(i,ydim-1)=var_a(i,ydim-2);
	}
	for (int j=1; j<(ydim-1); ++j) {
		var_a(0,j)=var_a(1,j);
		var_a(xdim-1,j)=var_a(xdim-2,j);
	}
	// avoid extremums in corner points
	var_a(0,0)=0.5*(var_a(1,0)+var_a(0,1));
	var_a(xdim-1,0)=0.5*(var_a(xdim-2,0)+var_a(xdim-1,1));
	var_a(0,ydim-1)=0.5*(var_a(1,ydim-1)+var_a(0,ydim-2));
	var_a(xdim-1,ydim-1)=0.5*(var_a(xdim-2,ydim-1)+var_a(xdim-1,ydim-2));
	m.remove_function_ifdef(dvar);
	return maxdvar;
}

/** find steady state solution of
 * d psi2 / dt = S(psi) (1- |grad(psi2)|), S(psi) = psi/\sqrt{psi^2+h_x^2+h_y^2}
 */
template<class fid_t>
int level_set_function_reset(AMesh2D<fid_t>& m) {
	m.remove_function_ifdef(PSI2);
	m.add_function_ifndef(PSI2);
	m.remove_function_ifdef(SPSI);
	m.add_function_ifndef(SPSI);
	double psi;
	// eps is heuristic parameter
	double eps=16*(m.get_dx()*m.get_dx()+m.get_dy()*m.get_dy());
	int xdim=m.get_xdim();
	int ydim=m.get_ydim();
	array2d psi_a=m[PSI];
	array2d spsi_a=m[SPSI];
	for (int i=0; i<xdim; ++i) {
		for (int j=0; j<ydim; ++j) {
			psi=psi_a(i,j);
			spsi_a(i,j)=psi/sqrt(psi*psi+eps);
//			psi=m.get(PSI,i,j);
//			m.set(SPSI,i,j,psi/sqrt(psi*psi+eps));
		}
	}
	int count=0;
	do {
		eps=level_set_function_reset_step<fid_t>(m,PSI2,DPSI2_DT);
		++count;
	} while (eps < Method::it().sle_solver_accuracy);
	array2d psi2_a=m[PSI2];
	for (int i=0; i<xdim; ++i) {
		for (int j=0; j<ydim; ++j) {
//			m.set(PSI,i,j,m.get(PSI2,i,j));
			psi_a(i,j)=psi2_a(i,j);
		}
	}
	m.remove_function_ifdef(PSI2);
	m.remove_function_ifdef(SPSI);
	return count;
}

/** time step for dpsi/dt + v*grad(psi) = 0
 *  with uniform boundary condition for psi (zero flux);
 *  WARNING: this implementation works for uniform rectangular grid only */
template<class fid_t>
AMesh2D<fid_t>*
step_level_set(double dt, const AMesh2D<fid_t>& m1)
throw(MeshException) {
	try {
	auto_ptr<AMesh2D<fid_t> > m2(m1.clone());
	eval_v<fid_t>(*m2);
	if (!m2->defined(VX) || !m2->defined(VY)) {
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
		count=level_set_function_reset<fid_t>(*m2);
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

// tepmlates
template
AMesh2D<int>*
phi_step<int>
(const Params& p, double dt, const AMesh2D<int>& m1, int const var);

template
AMesh2D<int>*
step_level_set<int>(double dt, const AMesh2D<int>& m1);

template
double
estimate_optimal_dt<int>(const AMesh2D<int>& m);
