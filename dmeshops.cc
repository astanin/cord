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

#include <math.h>

#include <vector>
#include <map>

using std::map;

#include "dmesh.h"
#include "dmeshops.h"
#include "global.h"
#include "utils.h"

extern string
dbg_stamp(double t);

DMesh
build_mesh(const Params& p) {
	DMesh m(p.xdim,p.ydim,0.0,p.xsize,0.0,p.ysize);
#ifdef HAVE_LIBHDF5
	if (p.inputfile.empty()) {
#endif
		// initial tissue density
		if (!m.defined("phi")) {
			m.add_function("phi",p.phi_stress_free);
		}
		// levelset potential
		if (!m.defined("psi")) {
			m.add_function("psi");
			for (int i=0; i<m.get_xdim(); ++i) {
				for (int j=0; j<m.get_ydim(); ++j) {
					double x, y;
					x=(m.x(i,j)-p.initial_cord_x)
						/p.initial_cord_length;
					y=(m.y(i,j)-p.initial_cord_y)
						/p.initial_cord_width;
					double levelset=1.0-sqrt(x*x+y*y);
					m.set("psi",i,j,levelset);
				}
			}
		}
		// initial nutrient distribution
		m.add_function("c",1.0);
		// model parameters
		m.set_attr("phi0", p.phi_stress_free);
		m.set_attr("o2_uptake",p.o2_uptake);
		m.set_attr("upkeep_per_cell",p.upkeep_per_cell);
		m.set_attr("death_rate",p.death_rate);
		m.set_attr("cell_motility", p.cell_motility);
		m.set_attr("tk1",p.tk1);
		m.set_attr("ts1",p.ts1);
		m.set_attr("hk1",p.hk1);
		m.set_attr("hs1",p.hs1);
		m.set_attr("host_activity",p.host_active?1.0:0.0);
#ifdef HAVE_LIBHDF5
	} else {
		m.load(p.inputfile);
	}
#endif
	return m;
}

#ifdef HAVE_LIBHDF5
using std::cerr;

int
hdf2dx(const Params& p) {
	DMesh m;
	if ((p.inputfile.empty())) {
		cerr << "should specify file for conversion "
			<< "(use -i filename.hdf)\n";
		return 1;
	}
	string convertedname=p.inputfile;
	int extpos=convertedname.find(".hdf");
	if ((unsigned)extpos != string::npos) {
		convertedname.replace(extpos,4,".dx");
	} else { // not found
		convertedname.append(".dx");
	}
	m.load(p.inputfile);
	dump2dx((AMesh2D&)m, convertedname);
	return 0;
}

int hdf2gp(const Params& p) {
	DMesh m;
	if ((p.inputfile.empty())) {
		cerr << "should specify file for conversion "
			<< "(use -i filename.hdf)\n";
		return 1;
	}
	string convertedname=p.inputfile;
	int extpos=convertedname.find(".hdf");
	if ((unsigned)extpos != string::npos) {
		convertedname.replace(extpos,4,".gp");
	} else { // not found
		convertedname.append(".gp");
	}
	m.load(p.inputfile);
	dump2gp(m, convertedname);
	return 0;
}
#endif // ifdef HAVE_LIBHDF5


