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

template int hdf2gp<int>(const Params& p);
template int hdf2dx<int>(const Params& p);

#endif // ifdef HAVE_LIBHDF5


