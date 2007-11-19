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

#include <stdlib.h>

#include "global.h"
#include "solver.h"
#include "tumeval.h"
#include "growth.h"
#include "nutrient.h"
#include "gfm.h"
#include "pradi.h"

using std::ostringstream;
using std::setw;
using std::setiosflags;
using std::ios_base;
using std::cerr;
using std::setprecision;
using std::ios;
using std::scientific;

#include <limits>
using std::numeric_limits;

#include <memory>
using std::auto_ptr;

#include "dmeshops.h"
#include "utils.h"

template<class fid_t>
double
estimate_optimal_dt(const AMesh2D<fid_t>& m);

/// @brief multicomponent tissue evaluation
/// evaluate behaviour of the tissiu compomnent @c density
/// in the subdomain where @c phase * @c where  > 0
template<class fid_t>
AMesh2D<fid_t>*
eval_bc_tissue(const Params& p, double const dt, const AMesh2D<fid_t>& m1,
	const fid_t& phi1, const fid_t& phi2, const fid_t& phi_h,
	const fid_t& phase, double const where) {
	// growth-death (ODE)
	auto_ptr<AMesh2D<fid_t> > mt(step_bc_tumour_growth_death<fid_t>
					(dt,m1, phi1,phi2,phase,where));
//	return mt.release();
	// PDE step
	auto_ptr<AMesh2D<fid_t> > m2(bc_phi_step<fid_t>
					(p,dt,*mt,phi1,phi2,phase,where));
	m2.reset(phi_step<fid_t>(p, dt, *m2, phi_h));
	return m2.release();
}

template<class fid_t>
AMesh2D<fid_t>*
eval_tissue(const Params& p, double const dt, const AMesh2D<fid_t>& m1,
	fid_t density_var) {
	// growth-death (ODE)
	auto_ptr<AMesh2D<fid_t> > tmp(step_growth_death<fid_t>(dt,m1,density_var));
	auto_ptr<AMesh2D<fid_t> > m2(phi_step<fid_t>(p, dt, *tmp, density_var));
	return m2.release();
}

void
validate_var(const string varname, array2d const& var,
	double const min, double const max, double const err=1e-6) {
	if ((blitz::max(var) > (max+err)) || (blitz::min(var) < (min-err))) {
		ostringstream ss;
		ss << "solution out of valid range: "
			<< "max(" << varname << ")= " << blitz::max(var) << " "
			<< "min(" << varname << ")= " << blitz::min(var);
		throw MeshException(ss.str());
	}
}

template<class fid_t>
AMesh2D<fid_t>*
solve(const Params& p, const AMesh2D<fid_t>& initial) {
	auto_ptr<AMesh2D<fid_t> > m1(initial.clone());
	double double_epsilon=numeric_limits<double>::epsilon();
	double last_dump_t=m1->get_time();
	double final_t=m1->get_time()+p.eval_t;
	double prevxsize=0.0;
	double prevysize=0.0;
	double eff_dt_initial=-1;
	bool use_ghostfluidmethod; // true if we need to use ghost fluid method
	bool use_bicomponenttissue; // true if tumour tissue is bicomponent
	if (m1->get_attr("conversion_rate") > double_epsilon) {
		use_bicomponenttissue=true;
	} else {
		use_bicomponenttissue=false;
	}
	if ((fabs(m1->get_attr("tk1")-m1->get_attr("hk1"))+
		fabs(m1->get_attr("ts1")-m1->get_attr("hs1"))) > 1e-99) {
		use_ghostfluidmethod=true;
		cerr << "Using Ghost Fluid method\n";
	} else {
		use_ghostfluidmethod=false;
		// we need these functions anyway to maintain save compatibility
		m1->add_function_ifndef(PHI_T);
		m1->add_function_ifndef(PHI_H);
		array2d phi_h=(*m1)[PHI_H];
		phi_h=(*m1)[PHI];
	}
	if (use_ghostfluidmethod && use_bicomponenttissue) {
		// unsupported combination of parameters
		throw MeshException("solve: different Sigmas are unsupported "
			"for glucose switch model");
	}
	TumourHostSplitter<fid_t> gfm(*m1);
	while (m1->get_time() < final_t) {
		double eff_dt=0.0;
		// if ((dt == 0) || (dt < 0)) use automatic time step
		if (p.dt < 1e-99) {
			eff_dt=estimate_optimal_dt<fid_t>(*m1);
			if (verbose) {
				cerr << dbg_stamp(m1->get_time())
					<< "autostep=" << eff_dt << "\n";
			}
		} else {
			eff_dt=p.dt;
		}
		if (eff_dt_initial < 0) {
			eff_dt_initial=eff_dt;
		}
		auto_ptr<AMesh2D<fid_t> > m2;
		if (use_ghostfluidmethod) {
			// update "real" (tumour) and "ghost" (host) fluids
			gfm.split(*m1,PHI,PSI,PHI_T,PHI_H);
			// evaluate tissue behaviour
			m2.reset(eval_tissue<fid_t>(p,eff_dt,*m1,PHI_T));
			m2.reset(eval_tissue<fid_t>(p,eff_dt,*m2,PHI_H));
			// re-construct solution
			gfm.merge(*m2,PHI,PSI,PHI_T,PHI_H);
		} else if (use_bicomponenttissue) {
			// 1. extrapolate PHI_H, PHI1, PHI2
			// 2. run modified eval_tissue for all of them
			// 3. construct new PHI from PHI_H, PHI1, PHI2
			m2.reset(extrapolate_subphases<fid_t>
					(*m1,PSI,PHI1,PHI2,PHI_H));
//			cerr << "phi1+phi2=" << m2->get(PHI1,16,0)+m2->get(PHI2,16,0) << " phi_h=" << m2->get(PHI_H,16,0) << " phi=" << m2->get(PHI,16,0) << " psi=" << m2->get(PSI,16,0) << "\n";
			m2.reset(eval_bc_tissue<fid_t>
					(p,eff_dt,*m2,PHI1,PHI2,PHI_H,PSI,+1));
//			cerr << "phi1+phi2=" << m2->get(PHI1,16,0)+m2->get(PHI2,16,0) << " phi_h=" << m2->get(PHI_H,16,0) << " phi=" << m2->get(PHI,16,0) << " psi=" << m2->get(PSI,16,0) << "\n";
			reconstruct_total_density<fid_t>
					(*m2,PSI,PHI,PHI1,PHI2,PHI_H);
//			cerr << "phi1+phi2=" << m2->get(PHI1,16,0)+m2->get(PHI2,16,0) << " phi_h=" << m2->get(PHI_H,16,0) << " phi=" << m2->get(PHI,16,0) << " psi=" << m2->get(PSI,16,0) << "\n";
		} else {
			// evaluate tissue behaviour
			m2.reset(eval_tissue<fid_t>(p,eff_dt,*m1,PHI));
		}
		// level set (interface tracking)
		m2.reset(step_level_set<fid_t>(eff_dt,*m2));
		// nutrient
		if (use_bicomponenttissue) {
			// glucose
			m2.reset(eval_bc_glc<fid_t>(p,*m2));
			// oxygen
			m2.reset(eval_bc_o2<fid_t>(p,*m2));
		} else {
			// only oxygen
			m2.reset(eval_nutrient<fid_t>(p,*m2,Method::it().
						p_solver_accuracy,eff_dt));
		}
		// done all sub steps
		m2->inc_time(eff_dt);
		m1.reset(m2->clone());
		if (verbose) {
			cerr << dbg_stamp(m1->get_time())
				<<setprecision(5)<<setiosflags(ios::scientific)
				<< "max(phi)= " << max((*m1)[PHI]) << " "
				<< "min(phi)= " << min((*m1)[PHI])
				<< "\n";
		}
		if (verbose > 1) { // extra verbose
			if (use_bicomponenttissue) {
			cerr << dbg_stamp(m1->get_time())
				<<setprecision(5)<<setiosflags(ios::scientific)
				<< "max(phi1)= " << max((*m1)[PHI1]) << " "
				<< "min(phi1)= " << min((*m1)[PHI1])
				<< "\n";
			cerr << dbg_stamp(m1->get_time())
				<<setprecision(5)<<setiosflags(ios::scientific)
				<< "max(phi2)= " << max((*m1)[PHI2]) << " "
				<< "min(phi2)= " << min((*m1)[PHI2])
				<< "\n";
			cerr << dbg_stamp(m1->get_time())
				<<setprecision(5)<<setiosflags(ios::scientific)
				<< "max(phi_h)= " << max((*m1)[PHI_H]) << " "
				<< "min(phi_h)= " << min((*m1)[PHI_H])
				<< "\n";
			}
			cerr << dbg_stamp(m1->get_time())
				<<setprecision(5)<<setiosflags(ios::scientific)
				<< "max(c)= " << max((*m1)[CO2]) << " "
				<< "min(c)= " << min((*m1)[CO2])
				<< "\n";
			if (use_bicomponenttissue) {
				cerr << dbg_stamp(m1->get_time())
				<<setprecision(5)<<setiosflags(ios::scientific)
				<< "max(c_g)= " << max((*m1)[GLC]) << " "
				<< "min(c_g)= " << min((*m1)[GLC])
				<< "\n";
			}
			cerr << dbg_stamp(m1->get_time())
				<<setprecision(5)<<setiosflags(ios::scientific)
				<< "max(psi)= " << max((*m1)[PSI]) << " "
				<< "min(psi)= " << min((*m1)[PSI])
				<< "\n";
			double xsize=get_x_size(*m1);
			double ysize=get_y_size(*m1);
			double vx=(xsize-prevxsize)/eff_dt;
			double vy=(ysize-prevysize)/eff_dt;
			cerr << dbg_stamp(m1->get_time())
				<<setprecision(2)<<setiosflags(ios::scientific)
				<< "xsize=" << xsize << " vx=" << vx << " "
				<< "ysize=" << ysize << " vy=" << vy << "\n";
			prevxsize=xsize;
			prevysize=ysize;
		}
		// TODO: reset variables nicely outside of their subdomains
		// dump state
		if (m1->get_time()>=(last_dump_t+p.dump_every-double_epsilon)){
			if (verbose) {
			cerr << dbg_stamp(m1->get_time()) << "dumping state\n";
			dump_mesh<fid_t>(*m1,(int)(floor(m1->get_time()*1e5)));
			last_dump_t=m1->get_time();
			}
		}
		// validate solution range
		validate_var("phi",(*m1)[PHI],0.0,1.0);
		validate_var("c",(*m1)[CO2],0.0,1.0,1e-2);
		if (use_bicomponenttissue) {
			validate_var("c_g",(*m1)[GLC],0.0,1.0,0.05);
			validate_var("phi1",(*m1)[PHI1],0.0,1.0,1e-2);
			validate_var("phi2",(*m1)[PHI2],0.0,1.0,1e-2);
		}
	}
	return m1.release();
}

template
AMesh2D<int>*
eval_tissue<int>
(const Params& p, double const dt, const AMesh2D<int>& m1, int density_var);

template
AMesh2D<int>*
solve<int>(const Params& p, const AMesh2D<int>& initial);
