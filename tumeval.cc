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
smart_grad(const AMesh2D<fid_t>& m, array2d const& u, const int i, const int j)
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
	g(0)=(u(i2,j)-u(i1,j))/(m.x(i2,j)-m.x(i1,j));
	g(1)=(u(i,j2)-u(i,j1))/(m.y(i,j2)-m.y(i,j1));
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
	array2d m_phi=m[PHI];
	array2d m_psi=m[PSI];
	array2d m_phi_sigma=m[TERM_PHI_SIGMA];
	array2d m_vx=m[VX];
	array2d m_vy=m[VY];
	TumourSigma<fid_t> tumour(m);
	auto_ptr<ADoubleFunction> t_sigma(tumour.build_sigma());
	HostSigma<fid_t> host(m);
	auto_ptr<ADoubleFunction> h_sigma(host.build_sigma());
	int xdim=m.get_xdim();
	int ydim=m.get_ydim();
	for (int i=0; i<xdim; ++i) {
		for (int j=0; j<ydim; ++j) {
			double phi=m_phi(i,j);
			double psi=m_psi(i,j);
			ADoubleFunction* sigma;
			if (psi>0) {
				sigma=t_sigma.get();
			} else {
				sigma=h_sigma.get();
			}
			m_phi_sigma(i,j)=phi*sigma->eval(phi);
		}
	}
	for (int i=0; i<xdim; ++i) {
		for (int j=0; j<ydim; ++j) {
			blitz::Array<double,1> g=
				smart_grad<int>(m,m_phi_sigma,i,j);
			double mu=m.get_attr("cell_motility");
			double vxval=(i>0)?(-g(0)*mu):0;
			double vyval=(j>0)?(-g(1)*mu):0;
			m_vx(i,j)=vxval;
			m_vy(i,j)=vyval;
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
	array2d m_vx=m[VX];
	array2d m_vy=m[VY];
	int xdim=m.get_xdim();
	int ydim=m.get_ydim();
	for (int i=0; i < xdim; ++i) {
		for (int j=0; j < ydim; ++j) {
			double vx, vy, v;
			vx=m_vx(i,j);
			vy=m_vy(i,j);
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
(AMesh2D<fid_t> const& m,
	array2d const& phi_a, array2d const& psi_a, double mu,
	const int i, const int j) {
	// WARNING: assuming that Sigmas do not change with time
	static ADoubleFunction *t_s=TumourSigma<fid_t>(m).build_sigma();
	static ADoubleFunction *t_s_p=TumourSigma<fid_t>(m).build_sigma_prime();
	static ADoubleFunction *h_s=HostSigma<fid_t>(m).build_sigma();
	static ADoubleFunction *h_s_p=HostSigma<fid_t>(m).build_sigma_prime();
	ADoubleFunction *sigma;
	ADoubleFunction *sigma_prime;
	double phi=phi_a(i,j);
	double psi=psi_a(i,j);
	if (psi > 0) {
		sigma=t_s;
		sigma_prime=t_s_p;
	} else {
		sigma=h_s;
		sigma_prime=h_s_p;
	}
	double pseudoD=mu*(phi*sigma->eval(phi)+phi*phi*sigma_prime->eval(phi));
	return pseudoD;
}

template<class fid_t>
void
phi_init_pseudo_D(AMesh2D<fid_t>& m, fid_t const var, fid_t const Dvar) {
	// init diffusion coefficient values
	m.remove_function_ifdef(Dvar);
	m.add_function_ifndef(Dvar,0.0);
	array2d m_D=m[Dvar];
	array2d m_phi=m[var];
	array2d m_psi=m[PSI];
	double mu=m.get_attr("cell_motility");
	for (int i=0; i < (m.get_xdim()); ++i) {
		for (int j=0; j < (m.get_ydim()); ++j) {
			m_D(i,j)=phi_pseudo_diff_coef(m,m_phi,m_psi,mu,i,j);
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
		ss << "phi_step: " << e.what();
		throw MeshException(ss.str());
	}
}

/** Pseudo diffusion coefficient for phi */
template<class fid_t>
double
bc_phi_pseudo_diff_coef(AMesh2D<fid_t>& m,
	array2d& phi1, array2d& phi2, double mu, int i, int j) {
	// WARNING: assuming that Sigmas do not change with time
	static ADoubleFunction *t_s=TumourSigma<fid_t>(m).build_sigma();
	static ADoubleFunction *t_s_p=TumourSigma<fid_t>(m).build_sigma_prime();
	double p1=phi1(i,j);
	double p2=phi2(i,j);
	double phi=p1+p2;
	double pseudoD=mu*p1*(t_s->eval(phi)+(phi)*t_s_p->eval(phi));
	return pseudoD;
}

template<class fid_t>
void
bc_phi_init_pseudo_D(AMesh2D<fid_t>& m, fid_t const phi1, fid_t const phi2,
	fid_t const Dvar) {
	// init diffusion coefficient values
	m.remove_function_ifdef(Dvar);
	m.add_function_ifndef(Dvar,0.0);
	array2d m_D=m[Dvar];
	array2d m_phi1=m[phi1];
	array2d m_phi2=m[phi2];
	array2d m_psi=m[PSI];
	double mu=m.get_attr("cell_motility");
	for (int i=0; i < (m.get_xdim()); ++i) {
		for (int j=0; j < (m.get_ydim()); ++j) {
			m_D(i,j)=bc_phi_pseudo_diff_coef<fid_t>
						(m,m_phi1,m_phi2,mu,i,j);
		}
	}
}

template<class fid_t>
AMesh2D<fid_t>*
bc_phi_step
(const Params& p, double dt, const AMesh2D<fid_t>& m1, fid_t const phi1,
	fid_t const phi2, fid_t const phase, double const where)
throw(MeshException) {
	try {
	auto_ptr<AMesh2D<fid_t> > m2(m1.clone());
	m2->remove_function_ifdef(PSEUDO_D);
	m2->add_function_ifndef(PSEUDO_D,0.0);
	bc_phi_init_pseudo_D<fid_t>(*m2,phi1,phi2,PSEUDO_D);
	m2.reset(reaction_diffusion_step<fid_t>(p.phi_bc,dt,*m2,phi1,
					PSEUDO_D,NONE,Method::it().rd_solver));
	bc_phi_init_pseudo_D<fid_t>(*m2,phi2,phi1,PSEUDO_D);
	m2.reset(reaction_diffusion_step<fid_t>(p.phi_bc,dt,*m2,phi2,
					PSEUDO_D,NONE,Method::it().rd_solver));
	return m2.release();
	} catch(MeshException& e) {
		ostringstream ss;
		ss << "bc_phi_step: " << e.what();
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
	array2d psi_a=m[PSI];
	array2d vx_a=m[VX];
	array2d vy_a=m[VY];
	array2d newpsi_a=m[NEWPSI];
	for (int i=1; i < (m.get_xdim()-1); ++i) {
		for (int j=1; j < (m.get_ydim()-1); ++j) {
			double psi, vx, vy;
			psi=psi_a(i,j);
			vx=vx_a(i,j);
			vy=vy_a(i,j);
			double vdotgradpsi;
			vdotgradpsi=0.5*(vx+fabs(vx))*
				(psi_a(i,j)-psi_a(i-1,j))/dx
				+0.5*(vx-fabs(vx))*
				(psi_a(i+1,j)-psi_a(i,j))/dx
				+0.5*(vy+fabs(vy))*
				(psi_a(i,j)-psi_a(i,j-1))/dy
				+0.5*(vy-fabs(vy))*
				(psi_a(i,j+1)-psi_a(i,j))/dy;
			newpsi_a(i,j)=psi-dt*vdotgradpsi;
		}
	}
	// apply boundary conditions
	// WARNING: valid for rectangular grid only
	for (int j=1; j<(m.get_ydim()-1); ++j) {
		double psi_in;
		psi_in=newpsi_a(1,j);
		newpsi_a(0,j)=psi_in;
		psi_in=newpsi_a(m.get_xdim()-2,j);
		newpsi_a(m.get_xdim()-1,j)=psi_in;
	}
	for (int i=0; i<m.get_xdim(); ++i) {
		double psi_in;
		psi_in=newpsi_a(i,1);
		newpsi_a(i,0)=psi_in;
		psi_in=newpsi_a(i,m.get_ydim()-2);
		newpsi_a(i,m.get_ydim()-1)=psi_in;
	}
	// overwrite old psi
	for (int i=0; i < (m.get_xdim()); ++i) {
		for (int j=0; j < (m.get_ydim()); ++j) {
			psi_a(i,j)=newpsi_a(i,j);
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

/**
 * Extrapolate var using
 *
 * d var
 * ----- + v * grad(var_ls)/|grad(var_ls)| * grad(var) = 0
 *  d t
 *
 * for t=1
 */
template<class fid_t>
AMesh2D<fid_t>*
extrapolate_var(const AMesh2D<fid_t>& m, fid_t var, fid_t var_ls, double v) {
	auto_ptr<AMesh2D<fid_t> > m2(m.clone());
//	cerr << dbg_stamp(m.get_time()) << "extrapolate_var: " << var << "\n";
	m2->remove_function_ifdef(TMP1);
	m2->add_function_ifndef(TMP1); // v*n_x
	m2->remove_function_ifdef(TMP2);
	m2->add_function_ifndef(TMP2); // v*n_y
	// init vector of normal to the interface
	int xdim=m.get_xdim();
	int ydim=m.get_ydim();
	array2d nx=(*m2)[TMP1];
	array2d ny=(*m2)[TMP2];
	array2d psi=(*m2)[var_ls];
	// init normal vecotr nx, ny
//	cerr << dbg_stamp(m.get_time()) << "extrapolate_var: init nx, ny\n";
	for (int i=1; i<(xdim-1); ++i) {
		for (int j=1; j<(ydim-1); ++j) {
			double gx=psi(i+1,j)-psi(i-1,j);
			double gy=psi(i,j+1)-psi(i,j-1);
			nx(i,j)=v*gx/sqrt(gx*gx+gy*gy);
			nx(i,j)=v*gy/sqrt(gx*gx+gy*gy);
		}
	}
	// evaluate
	double dx=m.get_dx();
	double dy=m.get_dy();
	double dt=0.5*fabs(std::min(m.get_dx(),m.get_dy())/v);
	double t=0.0;
	array2d u=m[var];
	array2d u2=(*m2)[var];
//	cerr << dbg_stamp(m.get_time()) << "extrapolate_var: dt=" << dt << "\n";
	while (t < 1.0) {
		// inner points
		for (int i=1; i<(xdim-1); ++i) {
			for (int j=1; j<(ydim-1); ++j) {
				if (v*psi(i,j) > 0) { // ignore the other domain
					double ax=nx(i,j);
					double ay=ny(i,j);
					double g= 0.5*(ax+fabs(ax))*
							(u(i,j)-u(i-1,j))/dx
						+ 0.5*(ax-fabs(ax))*
							(u(i+1,j)-u(i,j))/dx
						+ 0.5*(ay+fabs(ay))*
							(u(i,j)-u(i,j-1))/dy
						+ 0.5*(ay-fabs(ay))*
							(u(i,j+1)-u(i,j))/dy;
					u2(i,j)=u(i,j)+g*dt;
				}
			}
		}
		// boundaries (zero flux)
		for (int i=1; i<(xdim-1); ++i) {
			if (psi(i,0)*v > 0) {
				u2(i,0)=u2(i,1);
			}
			if (psi(i,ydim-1)*v > 0) {
				u2(i,ydim-1)=u2(i,ydim-2);
			}
		}
		for (int j=1; j<(ydim-1); ++j) {
			if (psi(0,j)*v > 0) {
				u2(0,j)=u2(1,j);
			}
			if (psi(xdim-1,j)*v > 0) {
				u2(xdim-1,j)=u2(xdim-2,j);
			}
		}
		// avoid extremums in corner points
		u2(0,0)=0.5*(u2(1,0)+u2(0,1));
		u2(xdim-1,0)=0.5*(u2(xdim-2,0)+u2(xdim-1,1));
		u2(0,ydim-1)=0.5*(u2(1,ydim-1)+u2(0,ydim-2));
		u2(xdim-1,ydim-1)=0.5*(u2(xdim-2,ydim-1)+u2(xdim-1,ydim-2));
		// increment time counter
		t+=dt;
	}
	m2->remove_function_ifdef(TMP1);
	m2->remove_function_ifdef(TMP2);
	return m2.release();
}

template<class fid_t>
AMesh2D<fid_t>*
extrapolate_subphases(const AMesh2D<fid_t>& m, fid_t var_ls,
	fid_t var_t1, fid_t var_t2, fid_t var_h) {
//	cerr << dbg_stamp(m.get_time()) << "extrapolate_subphases\n";
	auto_ptr<AMesh2D<fid_t> > m2(m.clone());
	array2d psi=(*m2)[var_ls];
	array2d phi_h=(*m2)[var_h];
	array2d phi2=(*m2)[var_t2];
	array2d phi1=(*m2)[var_t1];
	m2.reset(extrapolate_var(*m2,var_h, var_ls,+1.0));
	m2.reset(extrapolate_var(*m2,var_t1,var_ls,-1.0));
	m2.reset(extrapolate_var(*m2,var_t2,var_ls,-1.0));
	return m2.release();
}

template<class fid_t>
void
reconstruct_total_density(AMesh2D<fid_t>& m, fid_t var_ls,
	fid_t var, fid_t var_t1, fid_t var_t2, fid_t var_h) {
//	cerr << dbg_stamp(m.get_time()) << "reconstruct_total_density\n";
	int xdim=m.get_xdim();
	int ydim=m.get_ydim();
	array2d phi=m[var];
	array2d psi=m[var_ls];
	array2d phi1=m[var_t1];
	array2d phi2=m[var_t2];
	array2d phih=m[var_h];
	for (int i=0; i<xdim; ++i) {
		for (int j=0; j<ydim; ++j) {
			if (psi(i,j) > 0) { // tumour
				phi(i,j)=phi1(i,j)+phi2(i,j);
			} else {
				phi(i,j)=phih(i,j);
			}
		}
	}
}

// tepmlates
template
AMesh2D<int>*
phi_step<int>
(const Params& p, double dt, const AMesh2D<int>& m1, int const var);

template
AMesh2D<int>*
bc_phi_step<int>
(const Params& p, double dt, const AMesh2D<int>& m1, int const phi1,
	int const phi2, int const phase, double const where);

template
AMesh2D<int>*
step_level_set<int>(double dt, const AMesh2D<int>& m1);

template
double
estimate_optimal_dt<int>(const AMesh2D<int>& m);

template
AMesh2D<int>*
extrapolate_var(const AMesh2D<int>& m, int var, int var_ls, double v);

template
AMesh2D<int>*
extrapolate_subphases(const AMesh2D<int>& m, int var_ls,
	int var_t1, int var_t2, int var_h);

template
void
reconstruct_total_density(AMesh2D<int>& m, int var_ls,
	int var, int var_t1, int var_t2, int var_h);

template
double
bc_phi_pseudo_diff_coef(AMesh2D<int>& m,
	array2d& phi1, array2d& phi2,
	double mu, int i, int j);

