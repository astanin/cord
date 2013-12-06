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

#ifdef HAVE_LIBHDF5
using std::cerr;

template<class fid_t>
int
hdf2dx(const Params& p) {
	DMesh<fid_t> m;
	if ((p.inputfile.empty())) {
		cerr << "should specify file for conversion "
			<< "(use -i filename.h5)\n";
		return 1;
	}
	string convertedname=p.inputfile;
	int extpos=convertedname.find(".h5");
	if ((unsigned)extpos != string::npos) {
		convertedname.replace(extpos,4,".dx");
	} else { // not found
		convertedname.append(".dx");
	}
	m.load(p.inputfile);
	dump2dx<int>((AMesh2D<fid_t>&)m, convertedname);
	return 0;
}

template<class fid_t>
int hdf2gp(const Params& p) {
	DMesh<fid_t> m;
	if ((p.inputfile.empty())) {
		cerr << "should specify file for conversion "
			<< "(use -i filename.h5)\n";
		return 1;
	}
	string convertedname=p.inputfile;
	int extpos=convertedname.find(".h5");
	if ((unsigned)extpos != string::npos) {
		convertedname.replace(extpos,4,".gp");
	} else { // not found
		convertedname.append(".gp");
	}
	m.load(p.inputfile);
	dump2gp<int>(m, convertedname);
	return 0;
}

#endif // ifdef HAVE_LIBHDF5


extern "C" {
#ifdef HAVE_LIBNETPBM
#include <netpbm/pgm.h>
#endif
}


#ifdef HAVE_LIBNETPBM
template<class fid_t>
void
load_vasculature(DMesh<fid_t> &m, fid_t vasc_fid, string const filename) {
	FILE* f=fopen(filename.c_str(),"r");
	if (!f) {
		ostringstream ss;
		ss << "Unable to load vasculature from " << filename;
		throw MeshFileException(ss.str());
	}
	int cols, rows;
	gray maxval;
	gray** pgm=pgm_readpgm(f, &cols, &rows, &maxval);
	fclose(f);
	if ((cols != m.get_xdim()) || (rows != (m.get_ydim()))) {
		pgm_freearray(pgm,rows); pgm=(gray**)0;
		ostringstream ss;
		ss << filename << ": " << cols << "x" << rows
			<< "; mesh: " << m.get_xdim()  << "x" << m.get_ydim();
		throw MeshException(ss.str());
	}
	array2d vasc=m[vasc_fid];
	for (int j=0; j<rows; ++j) {
		for (int i=0; i<cols; ++i) {
			double val=1.0-pgm[j][i]/255.0;
			vasc(i,rows-1-j)=val;
//			if (val > 0.5) {
//				std::cerr << "·" ;
//			} else {
//				std::cerr << " " ;
//			}
		}
//		std::cerr << "\n";
	}
	pgm_freearray(pgm,rows); pgm=(gray**)0;
}

template int hdf2gp<int>(const Params& p);
template int hdf2dx<int>(const Params& p);

template void
load_vasculature<int>(DMesh<int> &m, int vasc_fid, string const filename);
#endif


