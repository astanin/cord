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

#include "params.h"
#include "dmesh.h"
#include "dmeshops.h"
#include "solver.h"

#include <iostream>

#include "utils.h"

using std::string;
using std::cerr;
using std::endl;
using std::ostringstream;
using std::setw;
using std::setfill;
using std::setiosflags;
using std::ios_base;


// global flags
int show_version=0;
int verbose=0;

void
print_model_params(const DMesh& m, const Params& p,
	int argc, const char *argv[]) {
	if (!verbose) {
		return;
	}
	// print command line (useful for logging)
	for (int i=0; i<argc; ++i) {
		cerr << argv[i] << " ";
	}
	cerr << "\n";
  	cerr<<"Evaluating from t0="<<setw(10)
		<<setiosflags(ios_base::left)<<m.get_time()<<"\n";
	cerr<<"\tfor t="<<p.eval_t<<" and dumping every "<<p.dump_every <<"\n";
	cerr<<"\twith dt: "<<p.dt<<"\n";
	cerr<<"\txdim x ydim: "<<m.get_xdim()<<"x"<<m.get_ydim()<<"\n";
	cerr<<"\tstress-free density: "<<m.get_attr("phi0")<<"\n";
	cerr<<"\ttumour sigma="<<m.get_attr("tk1")<<"*(phi-phi0),phi>phi0\n";
	cerr<<"\ttumour sigma="<<m.get_attr("ts1")<<"*(phi-phi0),phi<phi0\n";
	cerr<<"\thost sigma="<<m.get_attr("hk1")<<"*(phi-phi0),phi>phi0\n";
	cerr<<"\thost sigma="<<m.get_attr("hs1")<<"*(phi-phi0),phi<phi0\n";
	cerr<<"\tcell motility: "<<m.get_attr("cell_motility")<<"\n";
	cerr<<"\tcritical nutrient concentration: "
		<< m.get_attr("c_critical")<<"\n";
	cerr<<"\tbasic nutrient consumption: "
		<<m.get_attr("consumption_c")<<"\n";
	cerr<<"\tgrowth-related nutrient consumption: "
		<<m.get_attr("gconsumption_c")<<"\n";
	cerr<<"\tnutrient equation: ";
	switch(p.c_equation) {
		case Params::EQ_POISSON: cerr << "Poisson\n"; break;
		case Params::EQ_DIFFUSION: cerr << "diffusion\n"; break;
		default: cerr << "unknown\n"; break;
	}
#ifdef HAVE_LIBHDF5
	if (p.inputfile.empty()) {
#endif
		cerr<<"\tinitial cord length : "<< p.initial_cord_length<<"\n";
		cerr<<"\tinitial cord width : "<< p.initial_cord_width<<"\n";
#ifdef HAVE_LIBHDF5
	}
#endif
	vector<string> type_name(4);
	type_name[0]="unknown";
	type_name[1]="Dirichlet";
	type_name[2]="Neumann";
	type_name[3]="Robin";
	cerr<<"\tBoundary conditions:\n";
	cerr<<"\t\tnorth: phi: "<<type_name.at(p.phi_bc.get_north().get_type())
		<<" level=" << p.phi_bc.get_north().c() << "\n";
	cerr<<"\t\tnorth: c: "<<type_name.at(p.c_bc.get_north().get_type())
		<<" level=" << p.c_bc.get_north().c() << "\n";
	cerr<<"\t\teast: phi: "<<type_name.at(p.phi_bc.get_east().get_type())
		<<" level=" << p.phi_bc.get_east().c() << "\n";
	cerr<<"\t\teast: c: "<<type_name.at(p.c_bc.get_east().get_type())
		<<" level=" << p.c_bc.get_east().c() << "\n";
	cerr<<"\t\tsouth: phi: "<<type_name.at(p.phi_bc.get_south().get_type())
		<<" level=" << p.phi_bc.get_south().c() << "\n";
	cerr<<"\t\tsouth: c: "<<type_name.at(p.c_bc.get_south().get_type())
		<<" level=" << p.c_bc.get_south().c() << "\n";
	cerr<<"\t\twest: phi: "<<type_name.at(p.phi_bc.get_west().get_type())
		<<" level=" << p.phi_bc.get_west().c() << "\n";
	cerr<<"\t\twest: c: "<<type_name.at(p.c_bc.get_west().get_type())
		<<" level=" << p.c_bc.get_west().c() << "\n";
	cerr<<"\tNumerical methods:\n";
	cerr<<"\t\treaction-diffusion solver: ";
	switch (Method::it().rd_solver) {
	case MethodParams::RDS_IMPLICIT:
		cerr << "Euler implicit\n";
		break;
	case MethodParams::RDS_ADI:
		cerr << "alternative directions implicit\n";
		break;
	case MethodParams::RDS_EXPLICIT:
		cerr << "Euler explicit\n";
		break;
	default:
		cerr << "unknown\n";
		break;
	}
	cerr << "\t\tPoisson solver: ";
	switch (Method::it().p_solver) {
	case MethodParams::PS_RELAX:
		cerr << "relaxation to the steady state solution\n";
		break;
	case MethodParams::PS_SLE:
		cerr << "solution of the SLE\n";
		break;
	default:
		cerr << "unknown\n";
		break;
	}
	cerr << "\t\tSLE solver: ";
	switch(Method::it().sle_solver) {
	case MethodParams::SLES_UMFPACK:
		cerr<<"direct solution with UMFPACK\n";
		break;
	case MethodParams::SLES_CG:
		cerr<<"CG\n";
		break;
	case MethodParams::SLES_BICG:
		cerr<<"BICG\n";
		break;
	case MethodParams::SLES_BICGSTAB:
		cerr<<"BICGstab\n";
		break;
	case MethodParams::SLES_GMRES:
		cerr<<"GMRES (reset after "
			<< Method::it().sle_solver_gmres_restart_after <<")\n";
		break;
	default:
		cerr<<"unknown\n";
		break;
	};
}

int
main(int argc, const char *argv[]) {
	Params p;
	int ret=init_params(p,argc,argv);
	if (ret) {
		exit(EXIT_FAILURE);
	}
	if (show_version) {
		using std::cout;
		cout << PACKAGE << " " << VERSION "\n";
		cout << "Copyright (C) 2005-2007 Sergey Astanin, Luigi Preziosi, Andrea Tosin\n";
		cout << "This is free software; see the source for copying conditions.  There is NO\nwarranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n";
		exit(EXIT_SUCCESS);
	}
#ifdef HAVE_LIBHDF5
	if (p.hdf2dx) { // HDF5 to OpenDX conversion, no calculation
		ret=hdf2dx(p);
		if (ret) {
			exit(EXIT_FAILURE);
		} else {
			exit(EXIT_SUCCESS);
		}
	}
	if (p.hdf2gp) { // HDF5 to gnuplot conversion, no calculation
		ret=hdf2gp(p);
		if (ret) {
			exit(EXIT_FAILURE);
		} else {
			exit(EXIT_SUCCESS);
		}
	}
#endif
	try {
		DMesh m=build_mesh(p);
		print_model_params(m,p,argc,argv);
		DMesh* pm2=static_cast<DMesh*>(solve(p,m));
		m=*pm2;
#ifdef HAVE_LIBHDF5
		if (!p.outputfile.empty()) {
			m.save(p.outputfile);
		}
#endif
	} catch (MeshException &e) {
		cerr << "Exception: " << e.what() << endl;
		exit(EXIT_FAILURE);
	}
	exit(EXIT_SUCCESS);
}

