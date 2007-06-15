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

double
estimate_optimal_dt(const AMesh2D& m);

AMesh2D*
eval_tissue(const Params& p, double const dt, const AMesh2D& m1,
	string density_var) {
	// growth-death (ODE)
	auto_ptr<AMesh2D> tmp(step_growth_death(dt,m1,density_var));
	if (verbose > 1) { // extraverbose
		cerr << dbg_stamp(m1.get_time()) << "average_"
			<< density_var << "_growth="
			<< mean((*tmp)[density_var+"_growth"])
			<< "\n";
	}
	auto_ptr<AMesh2D> m2(phi_step(p, dt, *tmp, density_var));
	return m2.release();
}

AMesh2D*
solve(const Params& p, const AMesh2D& initial) {
	auto_ptr<AMesh2D> m1(initial.clone());
	double last_dump_t=m1->get_time();
	double final_t=m1->get_time()+p.eval_t;
	double prevxsize=0.0;
	double prevysize=0.0;
	double eff_dt_initial=-1;
	bool use_ghostfluidmethod; // true if we need to use ghost fluid method
	if ((fabs(m1->get_attr("tk1")-m1->get_attr("hk1"))+
		fabs(m1->get_attr("ts1")-m1->get_attr("hs1"))) > 1e-99) {
		use_ghostfluidmethod=true;
		cerr << "Using Ghost Fluid method\n";
	} else {
		use_ghostfluidmethod=false;
		// we need these functions anyway to maintain save compatibility
		m1->add_function_ifndef("phi_t");
		m1->add_function_ifndef("phi_h");
	}
	TumourHostSplitter gfm(*m1);
	while (m1->get_time() < final_t) {
		double eff_dt=0.0;
		// if ((dt == 0) || (dt < 0)) use automatic time step
		if (p.dt < 1e-99) {
			eff_dt=estimate_optimal_dt(*m1);
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
		auto_ptr<AMesh2D> m2;
		if (use_ghostfluidmethod) {
			// update "real" (tumour) and "ghost" (host) fluids
			gfm.split(*m1,"phi","psi","phi_t","phi_h");
			// evaluate tissue behaviour
			m2.reset(eval_tissue(p,eff_dt,*m1,"phi_t"));
			m2.reset(eval_tissue(p,eff_dt,*m2,"phi_h"));
			// re-construct solution
			gfm.merge(*m2,"phi","psi","phi_t","phi_h");
		} else {
			// evaluate tissue behaviour
			m2.reset(eval_tissue(p,eff_dt,*m1,"phi"));
		}
		// level set (interface tracking)
		m2.reset(step_level_set(eff_dt,*m2));
		// nutrient (eliptic)
		m2.reset(eval_nutrient(p,*m2,Method::it().p_solver_accuracy,eff_dt));
		// done all sub steps
		m2->inc_time(eff_dt);
		m1.reset(m2->clone());
		if (verbose) {
			cerr << dbg_stamp(m1->get_time())
				<<setprecision(5)<<setiosflags(ios::scientific)
				<< "max(phi)= " << max((*m1)["phi"]) << " "
				<< "min(phi)= " << min((*m1)["phi"])
				<< "\n";
		}
		if (verbose > 1) { // extra verbose
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
			cerr << dbg_stamp(m1->get_time())
				<<setprecision(5)<<setiosflags(ios::scientific)
				<< "max(c)= " << max((*m1)["c"]) << " "
				<< "min(c)= " << min((*m1)["c"])
				<< "\n";
			cerr << dbg_stamp(m1->get_time())
				<<setprecision(5)<<setiosflags(ios::scientific)
				<< "max(psi)= " << max((*m1)["psi"]) << " "
				<< "min(psi)= " << min((*m1)["psi"])
				<< "\n";
		}
		// dump state
		if (m1->get_time() >= (last_dump_t+p.dump_every
				-numeric_limits<double>::epsilon())) {
			if (verbose) {
			cerr << dbg_stamp(m1->get_time()) << "dumping state\n";
			dump_mesh(*m1,(int)(floor(m1->get_time()*1e5)));
			last_dump_t=m1->get_time();
			}
		}
		// validate solution range
		if ((max((*m1)["phi"]) > (1.0+1e-9)) ||
			(min((*m1)["phi"]) < (-1e-9)) ||
			(min((*m1)["c"]) < (-1e-9))) {
			ostringstream ss;
			ss << "solve: solution out of valid range: "
				<< "max(phi)= " << max((*m1)["phi"]) << " "
				<< "min(phi)= " << min((*m1)["phi"]) << " "
				<< "max(c)= " << max((*m1)["c"]) << " "
				<< "min(c)= " << min((*m1)["c"]);
			throw MeshException(ss.str());
		}
	}
	return m1.release();
}

