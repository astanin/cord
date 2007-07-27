/* Tumour cord growth model */

/*
 * Copyright (C) 2005-2006 Sergey Astanin, Luigi Preziosi
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

#include "utils.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using std::ios_base;
using std::ofstream;

// string case conversion
#include <string>
#include <algorithm>
#include <cctype>

#include <popt.h>
#include <iniparser.h>
#include "global.h"

using std::cerr;
using std::endl;
using std::ostringstream;
using std::setw;
using std::setfill;
using std::setiosflags;
using std::ios_base;

template<class fid_t>
string id2str(fid_t& id);

/// generates the string with timestamp which debug lines should begin with
string
dbg_stamp(double t) {
	using namespace std;
	ostringstream ss;
	ss << "t=" << setw(10) << setiosflags(ios_base::left) << t << ": ";
	return ss.str();
}

string
timestamp2str(int timestamp) {
	using namespace std;
	ostringstream ss;
	ss << setw(9) << setfill('0') << (timestamp);
	return ss.str();
}

BoundaryCondition
create_boundary_condition(const string bc_type, const double bc_val) {
	BoundaryCondition bc;
	if (bc_type == "Dirichlet") {
		bc=BoundaryCondition::createDirichletBC(bc_val);
	} else if (bc_type == "Neumann") {
		bc=BoundaryCondition::createNeumannBC(bc_val);
	} else {
		if (verbose) {
			cerr << "Unknown boundary condition type "
				<< bc_type << " ignored\n";
		}
	}
	return bc;
}

void
read_boundary_condition(dictionary *ini, Params& p, const string side) {
	using std::numeric_limits;
	char *type;
	double val;
	double nanv;
	if (numeric_limits<double>::has_quiet_NaN) {
		nanv=numeric_limits<double>::quiet_NaN();
	} else {
		nanv=1.0e199;
	}
	type=iniparser_getstring(ini,(side+":phi_condition_type").c_str(),0);
	val=iniparser_getdouble(ini,(side+":phi_condition_value").c_str(),nanv);
	if (type && (val != nanv)) {
		p.phi_bc.set_bc(side,create_boundary_condition(type,val));
	} else {
		cerr << "read_boundary_condition: phi conditions at "
			<< side << " not defined\n";
	}
	type=iniparser_getstring(ini,(side+":c_condition_type").c_str(),0);
	val=iniparser_getdouble(ini,(side+":c_condition_value").c_str(),nanv);
	if (type && (val != nanv)) {
		p.c_bc.set_bc(side,create_boundary_condition(type,val));
	} else {
		cerr << "read_boundary_condition: c conditions at "
			<< side << " not defined\n";
	}
}

void
read_solvers_config(dictionary *ini) {
	string s;
	// SLE solver
	s=iniparser_getstring(ini,"method:linear_solver","n/a");
#ifdef HAVE_LIBUMFPACK
	if (s == "umfpack") {
		Method::it().sle_solver=MP::SLES_UMFPACK;
	} else
#endif
	if (s == "gmres") {
		Method::it().sle_solver=MP::SLES_GMRES;
	} else if (s == "cg") {
		Method::it().sle_solver=MP::SLES_CG;
	} else if (s == "bicg") {
		Method::it().sle_solver=MP::SLES_BICG;
	} else if (s == "bicgstab") {
		Method::it().sle_solver=MP::SLES_BICGSTAB;
	}
	Method::it().sle_solver_accuracy=iniparser_getdouble(ini,
		"method:linear_solver_accuracy",
		Method::it().sle_solver_accuracy);
	Method::it().sle_solver_max_iters=iniparser_getint(ini,
		"method:linear_solver_max_iterations",
		Method::it().sle_solver_max_iters);
	Method::it().sle_solver_gmres_restart_after=iniparser_getint(ini,
		"method:gmres_restart_after",
		Method::it().sle_solver_gmres_restart_after);
	s=iniparser_getstring(ini,"method:reaction_diffusion_solver","n/a");
	if (s == "adi") {
		Method::it().rd_solver=MP::RDS_ADI;
	} else if (s == "implicit") {
		Method::it().rd_solver=MP::RDS_IMPLICIT;
	} else if (s == "explicit") {
		Method::it().rd_solver=MP::RDS_EXPLICIT;
	}
	s=iniparser_getstring(ini,"method:poisson_solver","n/a");
	if (s == "relax") {
		Method::it().p_solver=MP::PS_RELAX;
	} else if (s == "sle") {
		Method::it().p_solver=MP::PS_SLE;
	}
	Method::it().p_solver_accuracy=iniparser_getdouble(ini,
		"method:poisson_solver_accuracy",
		Method::it().p_solver_accuracy);
	Method::it().p_solver_relax_step=iniparser_getdouble(ini,
		"method:poisson_solver_relax_step",
		Method::it().p_solver_relax_step);
	Method::it().p_solver_relax_max_iters=iniparser_getint(ini,
		"method:poisson_solver_relax_max_iterations",
		Method::it().p_solver_relax_max_iters);
	Method::it().level_set_reset=iniparser_getboolean(ini,
		"method:level_set_reset",
		Method::it().level_set_reset);
}

// string case conversion
struct ToLower {
	char operator()(char const c) { return std::tolower(c); }
};

void
read_params(dictionary *ini, Params& p) {
	if (!ini) {
		return ;
	}
	p.xdim=iniparser_getint(ini,"params:xdim",p.xdim);
	p.ydim=iniparser_getint(ini,"params:ydim",p.ydim);
	p.xsize=iniparser_getdouble(ini,"params:domain_length",p.xsize);
	p.ysize=iniparser_getdouble(ini,"params:domain_width",p.ysize);
	p.initial_cord_length=iniparser_getdouble(ini,
		"params:initial_cord_length",p.initial_cord_length);
	p.initial_cord_width=iniparser_getdouble(ini,
		"params:initial_cord_width",p.initial_cord_width);
	p.phi_stress_free=iniparser_getdouble(ini,
		"params:phi_stress_free", p.phi_stress_free);
	p.cell_motility=iniparser_getdouble(ini,
		"params:phi_motility",p.cell_motility);
	p.o2_uptake=iniparser_getdouble(ini,
		"params:o2_uptake",p.o2_uptake);
	p.upkeep_per_cell=iniparser_getdouble(ini,
		"params:upkeep_per_cell",p.upkeep_per_cell);
	p.death_rate=iniparser_getdouble(ini,
		"params:death_rate",p.death_rate);
	p.host_active=iniparser_getboolean(ini,
		"params:active_host",p.host_active);
	p.tk1=iniparser_getdouble(ini,"params:tumour_stress_on_compress",p.tk1);
	p.ts1=iniparser_getdouble(ini,"params:tumour_stress_on_stretch",p.ts1);
	p.hk1=iniparser_getdouble(ini,"params:host_stress_on_compress",p.hk1);
	p.hs1=iniparser_getdouble(ini,"params:host_stress_on_stretch",p.hs1);
	// time integration related parameters
	p.dt=iniparser_getdouble(ini,"method:time_step",p.dt);
	p.eval_t=iniparser_getdouble(ini,"method:time_eval",p.eval_t);
	p.dump_every=iniparser_getdouble(ini,"method:dump_period",p.dump_every);
	// nutrient equation type
	string c_eq(iniparser_getstring(ini,"params:nutrient_equation",""));
	std::transform(c_eq.begin(),c_eq.end(),c_eq.begin(),ToLower());
	if (c_eq == "poisson") {
		p.c_equation=Params::EQ_POISSON;
	} else if (c_eq == "diffusion") {
		p.c_equation=Params::EQ_DIFFUSION;
	} else {
		// do not change default value
	}
}

void
read_conf_file(Params& p, const string filename) {
	dictionary *ini;
	ini = iniparser_load(filename.c_str());
	if (!ini) {
		string fallbackfilename=string(CONFDIR)+"/"+filename;
		ini = iniparser_load(fallbackfilename.c_str());
	}
	if (!ini) {
		cerr << "Configuration file \"" << filename << "\" not found\n";
		return;
	}
	// read mesh configuration
	read_params(ini,p);
	// read boundary conditions
	read_boundary_condition(ini,p,"north");
	read_boundary_condition(ini,p,"east");
	read_boundary_condition(ini,p,"south");
	read_boundary_condition(ini,p,"west");
	// read numerical method details
	read_solvers_config(ini);
	iniparser_freedict(ini);
}

int init_params(Params& p, int argc, const char *argv[]) {
#ifdef HAVE_LIBHDF5
	char *ifile=(char*)0;
	char *ofile=(char*)0;
#endif
	char *c_eq=(char*)0;
	char *conffile=(char*)0;
	int arg_val_hidden_arg=(POPT_ARG_VAL|POPT_ARGFLAG_ONEDASH)
					&(~POPT_ARGFLAG_SHOW_DEFAULT);
	int onedash=POPT_ARGFLAG_ONEDASH;
	struct poptOption optionsTable[] = {
		{ "config", 'c', POPT_ARG_STRING|onedash, &conffile, 0,
			"use alternative configuration file", "filename" },
#ifdef HAVE_LIBHDF5
		{ "hdf2dx", 0, arg_val_hidden_arg, &p.hdf2dx, 1,
			"convert HDF5 data to OpenDX format", 0 },
		{ "hdf2gp", 0, arg_val_hidden_arg, &p.hdf2gp, 1,
			"convert HDF5 data to gnuplot format", 0 },
		{ "load", 'i', POPT_ARG_STRING|onedash, &ifile, 0,
			"start with given initial conditions", "filename.hdf" },
		{ "save", 'o', POPT_ARG_STRING|onedash, &ofile, 0,
			"save the final state as HDF5 file","filename.hdf"},
#endif // ifdef HAVE_LIBHDF5
		{ "dt", 'd', POPT_ARG_DOUBLE|onedash, &p.dt, 0,
			"evaluate with specified time step", "timestep" },
		{ "period", 't', POPT_ARG_DOUBLE|onedash,
			&p.eval_t, 0,
			"evaluate during specified time", "time" },
		{ "dump-period", 0, POPT_ARG_DOUBLE|onedash,
			&p.dump_every, 0,
			"save the state after every time interval", "time" },
		{ "xdim", 0, POPT_ARG_INT|onedash, &p.xdim, 0,
			"grid dimension along x-axis", "points" },
		{ "ydim", 0, POPT_ARG_INT|onedash, &p.ydim, 0,
			"grid dimension along y-axis", "points" },
		{ "length", 0, POPT_ARG_DOUBLE|onedash,
			&p.xsize, 0, "length of the growth domain", "xsize" },
		{ "width", 0, POPT_ARG_DOUBLE|onedash,
			&p.ysize, 0, "width of the growth domain", "ysize" },
		{ "cord-length", 0, POPT_ARG_DOUBLE|onedash,
			&p.initial_cord_length, 0,
			"initial length of the cord", "xsize" },
		{ "cord-width", 0, POPT_ARG_DOUBLE|onedash,
			&p.initial_cord_width, 0,
			"initial width of the cord", "ysize" },
		{ "cord-x", 0, POPT_ARG_DOUBLE|onedash,
			&p.initial_cord_x, 0,
			"initial position (x) of the cord", "x" },
		{ "cord-y", 0, POPT_ARG_DOUBLE|onedash,
			&p.initial_cord_y, 0,
			"initial position (y) of the cord", "y" },
		{ "phi-stress-free", 0, POPT_ARG_DOUBLE|onedash,
			&p.phi_stress_free, 0,
			"density of the stress-free tissue", "phi0" },
		{ "phi-motility", 0, POPT_ARG_DOUBLE|onedash,
			&p.cell_motility, 0, "motility (relaxation rate)"
			" of the cellular phase", "mu" },
		{ "tumour-compress", 0, POPT_ARG_DOUBLE|onedash, &p.tk1, 0,
			"tumour stress: sigma=tk1*(phi-phi0), for phi>phi0",
	       		"tk1" },
		{ "tumour-stretch", 0, POPT_ARG_DOUBLE|onedash, &p.ts1, 0,
			"tumour stress: sigma=ts1*(phi-phi0), for phi<phi0",
			"ts1"},
		{ "host-compress", 0, POPT_ARG_DOUBLE|onedash, &p.hk1, 0,
			"host stress: sigma=hk1*(phi-phi0), for phi>phi0",
			"hk1"},
		{ "host-stretch", 0, POPT_ARG_DOUBLE|onedash, &p.hs1, 0,
			"host stress: sigma=hs1*(phi-phi0), for phi<phi0",
			"hs1"},
		{ "nutrient-equation", 0, POPT_ARG_STRING|onedash, &c_eq, 0,
			"equation type: Poisson | diffusion", "type"},
		{ "upkeep-per-cell", 'u', POPT_ARG_DOUBLE|onedash,
			&p.upkeep_per_cell, 0,
			"minimal energy consumption by cells", "theta" },
		{ "death-rate", 'e', POPT_ARG_DOUBLE|onedash,
			&p.death_rate, 0,
			"rate of cellular death", "epsilon" },
		{ "O2-uptake", 'a', POPT_ARG_DOUBLE|onedash,
			&p.o2_uptake, 0,
			"rate of oxygen consumption", "alpha" },
		{ "active-host", 0, arg_val_hidden_arg, &p.host_active, true,
			"host tissue consumes oxygen and may die", 0 },
		{ "verbose", 'v', arg_val_hidden_arg, &verbose, 1,
			"be verbose", 0},
		{ "extraverbose", 0, arg_val_hidden_arg, &verbose, 2,
			"be extra verbose (print really much)", 0},
		{ "version", 0, POPT_ARG_NONE|onedash,
			&show_version, 0, "show version", 0 },
		POPT_AUTOHELP
		{ NULL, 0, 0, NULL, 0, 0 } };
	poptContext con=poptGetContext(NULL, argc, argv, optionsTable, 0);
	int rc=0;
	rc = poptGetNextOpt(con);
	if (conffile) {
		read_conf_file(p,conffile);
	} else {
		read_conf_file(p,"cord.ini");
	}
	// command line options have to override the config file
	poptResetContext(con);
	rc = poptGetNextOpt(con);
	if (rc < -1) {
		cerr << poptBadOption(con, POPT_BADOPTION_NOALIAS) << ": "
			<< poptStrerror(rc) << "\n";
		poptFreeContext(con);
		return 1;
	} else {
#ifdef HAVE_LIBHDF5
		if (ifile) {
			p.inputfile=string(ifile);
		}
		if (ofile) {
			p.outputfile=string(ofile);
		}
#endif
		poptFreeContext(con);
		if (c_eq) {
			string s(c_eq);
			std::transform(s.begin(),s.end(),s.begin(),ToLower());
			if (s == "poisson") {
				p.c_equation=Params::EQ_POISSON;
			} else if (s == "diffusion") {
				p.c_equation=Params::EQ_DIFFUSION;
			} else {
				cerr << "-nutrient-equation=" << c_eq
					<< ": wrong equation type\n";
				return 1;
			}
		}
		return 0;
	}
}

template<class fid_t>
void
dump2gp(AMesh2D<fid_t>& m, string const& filename) {
	ofstream fs;
	fs.open(filename.c_str());
	if (fs.is_open()) {
		fs.setf(ios_base::scientific);
		fs.precision(6);
		fs << "#time: "<<m.get_time()<<"\n";
		fs << "#dimesions: "<<m.get_xdim()<<" "<<m.get_ydim()<<"\n";
		fs << "#x:";
		for (int i=0; i<m.get_xdim(); ++i) {
			fs << " " << m.x(i,0);
		}
		fs << "\n";
		fs << "#y:";
		for (int j=0; j<m.get_ydim(); ++j) {
			fs << " " << m.y(0,j);
		}
		fs << "\n";
		vector<fid_t> fids;
		fids=m.get_fids();
		typename vector<fid_t>::iterator fun;
		fs << "#x #y ";
		for (fun=fids.begin(); fun!=fids.end(); ++fun) {
			if (fun != fids.begin()) {
				fs << " ";
			}
			fs << "#" << id2str(*fun);
		}
		fs << "\n";
		for (int i=0; i<m.get_xdim(); ++i) {
			for (int j=0; j<m.get_ydim(); ++j) {
				fs << m.x(i,j) << " " << m.y(i,j);
				for (fun=fids.begin(); fun!=fids.end(); ++fun) {
					fs << " " << m.get(*fun,i,j);
				}
				fs << "\n";
			}
			fs << "\n";
		}
		fs.close();
	} else {
		throw MeshException("cannot open dump file");
	}
}

template<class fid_t>
void
dump2gp(AMesh2D<fid_t>& m, int const timestamp) {
	ostringstream ss;
	ss << "dmp" << timestamp2str(timestamp) << ".gp";
	dump2gp(m, ss.str());
}

template<class fid_t>
int
count_inner_points(AMesh2D<fid_t> const& m) {
	int count=0;
	for (int i=0; i<m.get_xdim(); ++i) {
		for (int j=0; j<m.get_ydim(); ++j) {
			if (m.is_inner(i,j)) {
				++count;
			}
		}
	}
	return count;
}


/** count quad-connections with bottom-left angle in (i,j)
 * as well as semi-quad-connections (triangular) on the border */
template<class fid_t>
int
count_connections(AMesh2D<fid_t> const& m) {
	int count=0;
	for (int i=0; i<(m.get_xdim()-1); ++i) {
		for (int j=0; j<(m.get_ydim()-1); ++j) {
			if (m.is_inner(i,j) && m.is_inner(i+1,j)
				&& m.is_inner(i,j+1)) {
				// there is either full-quad or triangular
				// connection between (i,j), (i+1,j), and
				// (i,j+1)
				++count;
			}
		}
	}
	return count;
}

template<class fid_t>
void dump2dx_scalar_field
(AMesh2D<fid_t> const& m, fid_t const& fid, ofstream& fs) {
	int xdim=m.get_xdim();
	int ydim=m.get_ydim();
	int npos=xdim*ydim; // WARNING: assuming rectangular grid
	if (m.defined(fid)) {
		// nutrient data
		fs << "# " << fid << " data\n";
		fs << "object \"" << fid << "_data\" class array type float "
			<< "rank 0 items " << npos << " data follows\n";
		for (int i=0; i<xdim; ++i) {
			for (int j=0; j<ydim; ++j) {
				if (m.is_inner(i,j)) {
					fs << "\t" << m.get(fid,i,j) << "\n";
				}
			}
		}
		fs << "attribute \"dep\" string \"positions\"\n";
		// constructing field
		fs << "# data field\n";
		fs << "object \"" << fid << "\" class field\n";
		fs << "component \"positions\" value 1\n";
		fs << "component \"connections\" value 2\n";
		fs << "component \"data\" value \"" << fid << "_data\"\n";
		fs << "attribute \"time\" number " << m.get_time() << "\n";
		vector<string> attrs=m.get_attrs();
		vector<string>::const_iterator ai=attrs.begin();
		for (; ai != attrs.end(); ++ai) {
			double value=m.get_attr(*ai);
			fs << "attribute \"" << *ai
				<< "\" number " << value << "\n";
		}
	}
}

template<class fid_t>
void dump2dx_vector_field
(AMesh2D<fid_t> const& m, fid_t const& vxfid, fid_t const& vyfid,
	string const& fieldname, ofstream& fs) {
	int xdim=m.get_xdim();
	int ydim=m.get_ydim();
	int npos=xdim*ydim; // WARNING: assuming rectangular grid
	string v(fieldname);
	if (m.defined(vxfid) && m.defined(vyfid)) {
		fs << "# ( " << v << "x, " << v << "y ) data\n";
		fs << "object \"" << v << "_data\" class array type float "
			<< "rank 1 shape 2 items " << npos
			<< " data follows\n";
		for (int i=0; i<xdim; ++i) {
			for (int j=0; j<ydim; ++j) {
				if (m.is_inner(i,j)) {
					fs << "\t" << m.get(vxfid,i,j)
						<< "\t"
						<< m.get(vyfid,i,j)
						<< "\n";
				}
			}
		}
		fs << "attribute \"dep\" string \"positions\"\n";
		// constructing field
		fs << "# data field\n";
		fs << "object \"" << v << "\" class field\n";
		fs << "component \"positions\" value 1\n";
		fs << "component \"connections\" value 2\n";
		fs << "component \"data\" value \"v_data\"\n";
		fs << "attribute \"time\" number " << m.get_time() << "\n";
		vector<string> attrs=m.get_attrs();
		vector<string>::const_iterator ai=attrs.begin();
		for (; ai != attrs.end(); ++ai) {
			double value=m.get_attr(*ai);
			fs << "attribute \"" << *ai
				<< "\" number " << value << "\n";
		}
	}
}

template<class fid_t>
void
dump2dx(AMesh2D<fid_t> const& m, string const& filename) {
	ofstream fs;
	fs.open(filename.c_str());
	if (fs.is_open()) {
		int npos=count_inner_points(m);
		int ncon=count_connections(m);
		int xdim=m.get_xdim();
		int ydim=m.get_ydim();
		vector<int> position(xdim*ydim, -1);
		fs.setf(ios_base::fixed);
		fs.precision(6);
		// write positions, save position numbers in @c position
		fs << "# irregular positions - inner cord points\n";
		fs << "object 1 class array type float "
			<< "rank 1 shape 2 items " << npos << " data follows\n";
		int posno=0; // current position number
		for (int i=0; i<xdim; ++i) {
			for (int j=0; j<ydim; ++j) {
				if (m.is_inner(i,j)) {
					fs << "\t" << m.x(i,j)
						<< "\t" << m.y(i,j) << "\n";
					// starting from 0
					position.at(i+j*xdim)=posno;
					++posno;
				}
			}
		}
		// write connections
		fs << "# irregular connections, quads and half-quads\n";
		fs << "object 2 class array type int "
			<< "rank 1 shape 4 items " << ncon << " data follows\n";
		for (int i=0; i<xdim; ++i) {
			for (int j=0; j<ydim; ++j) {
				if (m.is_inner(i,j) && m.is_inner(i+1,j)
					&& m.is_inner(i,j+1)) {
					// list i,j;i+1,j;i+1,j+1 or i+1,j;i,j+1
					fs << "\t" << position[i+j*xdim] << "\t"
						<< position[(i+1)+j*xdim]
						<< "\t";
					fs << position[i+(j+1)*xdim] << "\t";
					if (m.is_inner(i+1,j+1)) {
						fs << position[(i+1)+(j+1)
							*xdim];
					} else { // repeat second vertex twice
						fs << position[(i+1)+j*xdim];
					}
					fs << "\n";
				}
			}
		}
		fs << "attribute \"element type\" string \"quads\"\n";
		fs << "attribute \"ref\" string \"positions\"\n";
		
		dump2dx_scalar_field<fid_t>(m,PHI,fs);
		dump2dx_scalar_field<fid_t>(m,CO2,fs);
		dump2dx_vector_field<fid_t>(m,VX,VY,"v",fs);
		dump2dx_scalar_field<fid_t>(m,PSI,fs);
		dump2dx_scalar_field<fid_t>(m,PHI_T,fs);
		dump2dx_scalar_field<fid_t>(m,PHI_H,fs);

		fs << "end\n";
	} else {
		ostringstream ss;
		ss << "dump2dx: cannot open file `" << filename << "'";
		throw MeshException(ss.str());
	}
}

template<class fid_t>
void
dump2dx(AMesh2D<fid_t> const& m, int const timestamp) {
	ostringstream ss;
	ss << "dmp" << timestamp2str(timestamp) << ".dx";
	dump2dx<fid_t>(m, ss.str());
}

template<class fid_t>
void
dump_mesh(AMesh2D<fid_t>& m, int timestamp) {
#ifdef HAVE_LIBHDF5
	// Save complete HDF file
	ostringstream ss;
	ss << "dmp" << timestamp2str(timestamp) << ".hdf";
	m.save(ss.str());
#else
	// Dump to OpenDX
	dump2dx<fid_t>(m,timestamp);
	// Dump to Gnuplot
	dump2gp<fid_t>(m,timestamp);
#endif
}


/** Evaluate current size of the cord along x-axis (based on psi).
 * Assuming that psi>0 inside the cord. */
template<class fid_t>
double get_x_size(const AMesh2D<fid_t>& m) {
	double psi1=0, psi2=0;
	psi1=m.get(PSI,0,0);
	for (int i=1; i<m.get_xdim(); ++i) {
		psi2=m.get(PSI,i,0);
		if ((psi1>=0) && (psi2<0)) {
			double x1, x2;
			double stepratio;
			double x0;
			x1=m.x(i-1,0);
			x2=m.x(i,0);
			stepratio=psi1/(fabs(psi1)+fabs(psi2));
			x0=x1+stepratio*(x2-x1);
			return x0;
		}
		psi1=psi2;
		psi2=1e99; // reset
	}
	if (fabs(psi1) < 1e-99) { // almost zero
		return m.x(m.get_xdim()-1,0);
	} else {
		return 1e99; // infinite size, out of range value
	}
}

/** Evaluate current size of the cord along y-axis (based on psi).
 * Assuming that psi>0 inside the cord. */
template<class fid_t>
double get_y_size(const AMesh2D<fid_t>& m) {
	double psi1=0, psi2=0;
	psi1=m.get(PSI,0,0);
	for (int j=1; j<(m.get_ydim()); ++j) {
		psi2=m.get(PSI,0,j);
		if ((psi1>=0) && (psi2<0)) {
			double y1, y2;
			double stepratio;
			double y0;
			y1=m.y(0,j-1);
			y2=m.y(0,j);
			stepratio=psi1/(fabs(psi1)+fabs(psi2));
			y0=y1+stepratio*(y2-y1);
			return y0;
		}
		psi1=psi2;
		psi2=1e99; // reset
	}
	if (fabs(psi1) < 1e-99) { // almost zero
		return m.y(0,m.get_ydim()-1);
	} else {
		return 1e99; // infinite size, out of range value
	}
}

/// norm-I: max(abs(f))
template<class fid_t>
double norm_1(AMesh2D<fid_t> const& m, fid_t const fid) {
	using namespace blitz;
	return max(fabs(m[fid]));
}

/// norm-II: sum(abs(f))
template<class fid_t>
double norm_2(AMesh2D<fid_t> const& m, fid_t const fid) {
	using namespace blitz;
	double norm=sum(fabs(m[fid]));
	return norm;
}

// templates

template double get_x_size<int>(const AMesh2D<int>& m);
template double get_y_size<int>(const AMesh2D<int>& m);
template void dump_mesh<int>(AMesh2D<int>& m, int timestamp);
template double norm_1<int>(AMesh2D<int> const& m, int const fid);
template void dump2dx<int>(AMesh2D<int> const& m, string const& filename);
template void dump2gp<int>(AMesh2D<int>& m, string const& filename);

