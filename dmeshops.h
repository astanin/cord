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

#ifndef DMESHOPS_H
#define DMESHOPS_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <iomanip>
#include <sstream>
#include "params.h"

using std::string;
using std::vector;
using std::ostringstream;

string timestamp2str(int timestamp);

// forward declaration
template<class fid_t> class DMesh;

/** @brief construct DMesh satisfying given parameters */
template<class fid_t>
DMesh<fid_t>
build_mesh(const Params& p);

template<class fid_t>
int
hdf2dx(const Params& p);

template<class fid_t>
int
hdf2gp(const Params& p);

#include <limits>
using std::numeric_limits;

template<class fid_t>
DMesh<fid_t>
build_mesh(const Params& p) {
	DMesh<fid_t> m(p.xdim,p.ydim,0.0,p.xsize,0.0,p.ysize);
#ifdef HAVE_LIBHDF5
	if (p.inputfile.empty()) {
#endif
		// initial tissue density
		if (!m.defined(PHI)) {
			m.add_function(PHI,p.phi_stress_free);
		}
		if (p.conversion_rate > numeric_limits<double>::epsilon()) {
			if (!m.defined(PHI1)) { // first subpopulation
				m.add_function(PHI1,0.0);
				array2d phi=m[PHI];
				array2d phi1=m[PHI1];
				phi1=phi;
			}
			if (!m.defined(PHI2)) { // second subpopulation
				m.add_function(PHI2,0.0);
			}
			if (!m.defined(GLC)) { // glucose
				m.add_function(GLC,0.0);
			}
		}
		// levelset potential
		if (!m.defined(PSI)) {
			m.add_function(PSI);
			for (int i=0; i<m.get_xdim(); ++i) {
				for (int j=0; j<m.get_ydim(); ++j) {
					double x, y;
					x=(m.x(i,j)-p.initial_cord_x)
						/p.initial_cord_length;
					y=(m.y(i,j)-p.initial_cord_y)
						/p.initial_cord_width;
					double levelset=1.0-sqrt(x*x+y*y);
					m.set(PSI,i,j,levelset);
				}
			}
		}
		// initial nutrient distribution
		m.add_function(CO2,1.0);
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
		m.set_attr("conversion_rate",p.conversion_rate);
		m.set_attr("anaerobic_rate",p.anaerobic_rate);
#ifdef HAVE_LIBHDF5
	} else {
		m.load(p.inputfile);
	}
#endif
	return m;
}

#endif /* DMESHOPS_H */

