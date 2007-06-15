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

#include <stdlib.h>

#include <math.h>
#include "dmesh.h"
#include "gfm.h"

using std::make_pair;

DMesh::DMesh(int _xdim, int _ydim, double _xmin, double _xmax, double _ymin,
	double _ymax) :
	time(0), xdim(_xdim), ydim(_ydim), xmin(_xmin),
	xmax(_xmax), ymin(_ymin), ymax(_ymax), m_x(_xdim,_ydim),
	m_y(_xdim, _ydim), mf(), attr() {
		init_mesh();
}

DMesh::DMesh(const DMesh& om) :
	time(om.time), xdim(om.xdim), ydim(om.ydim), xmin(om.xmin),
	xmax(om.xmax), ymin(om.ymin), ymax(om.ymax),
	m_x(om.xdim,om.ydim), m_y(om.xdim,om.ydim), mf(), attr() {
	/* copy mesh points */
	m_x.resize(om.m_x.shape());
	m_x=om.m_x.copy();
	m_y.resize(om.m_y.shape());
	m_y=om.m_y.copy();
	/* copy arrays */
	mf.clear();
	map<string,array2d>::const_iterator i = om.mf.begin();
	for (; i != om.mf.end(); ++i) {
		string mfkey=i->first;
		array2d mfc;
		mfc.resize(i->second.shape());
		mfc=i->second.copy();
		mf.insert(make_pair(mfkey,mfc));
	}
	attr=om.attr;
}

DMesh&
DMesh::operator=(const DMesh& rhs) {
	if (this == &rhs) { // do not do anything
		return *this;
	} else {
		time=rhs.time;
		xdim=rhs.xdim;
		ydim=rhs.ydim;
		xmin=rhs.xmin;
		xmax=rhs.xmax;
		ymin=rhs.ymin;
		ymax=rhs.ymax;
		/* copy mesh points */
		m_x.resize(rhs.m_x.shape());
		m_x=rhs.m_x.copy();
		m_y.resize(rhs.m_y.shape());
		m_y=rhs.m_y.copy();
		/* copy arrays */
		mf.clear();
		map<string,array2d>::const_iterator i=rhs.mf.begin();
		for (; i != rhs.mf.end(); ++i) {
			string mfkey=i->first;
			array2d mfc;
			mfc.resize(i->second.shape());
			mfc=i->second.copy();
			mf.insert(make_pair(mfkey,mfc));
		}
		attr=rhs.attr;
		return *this;
	}
}

double
DMesh::x(const int i, const int j)
const throw(MeshException) {
	if (!is_out_of_domain(i,j)) {
		return m_x(i,j);
	} else {
		ostringstream ss;
		ss << "DMesh::x: indices out of range: " << i << ", " << j ;
		throw MeshException(ss.str());
	}
}

double
DMesh::y(const int i, const int j)
const throw(MeshException) {
	if (!is_out_of_domain(i,j)) {
		return m_y(i,j);
	} else {
		ostringstream ss;
		ss << "DMesh::y: indices out of range: " << i << ", " << j ;
		throw MeshException(ss.str());
	}
}

void
DMesh::set_x(const int i, const int j, double const value)
throw(MeshException) {
	if (!is_out_of_domain(i,j)) {
		m_x(i,j)=value;
	} else {
		ostringstream ss;
		ss << "DMesh::set_x: indices out of range: " << i << ", " << j ;
		throw MeshException(ss.str());
	}
}

void
DMesh::set_y(const int i, const int j, double const value)
throw(MeshException) {
	if (!is_out_of_domain(i,j)) {
		m_y(i,j)=value;
	} else {
		ostringstream ss;
		ss << "DMesh::set_y: indices out of range: " << i << ", " << j ;
		throw MeshException(ss.str());
	}
}

void
DMesh::init_mesh(void) {
	if ((xdim < 2) || (ydim < 2)) {
		ostringstream sstr;
		sstr << "init_mesh(): cannot map "
				      << xdim << "*" << ydim << " mesh";
		throw MeshException(sstr.str());
	}
	double dx=(xmax-xmin)/(xdim-1);
	double dy=(ymax-ymin)/(ydim-1);
	for (int i=0; i<xdim; i++) {
		for (int j=0; j<ydim; j++) {
			set_x(i,j,xmin+i*dx);
			set_y(i,j,ymin+j*dy);
		}
	}
}

void
DMesh::resize(int const _xdim, int const _ydim) {
	if ((xdim == _xdim) && (ydim == _ydim)) { // do nothing
		return;
	} else {
		xdim=_xdim;
		ydim=_ydim;
		m_x.resize(_xdim,_ydim);
		m_y.resize(_xdim,_ydim);
		mf.clear();
	}
}

void
DMesh::add_function(string const fid, double const value)
throw() {
	array2d a;
	mf[fid]=a;
	mf[fid].resize(xdim,ydim);
	mf[fid]=value;
}

void
DMesh::add_function_ifndef(string const fid, double const value)
throw() {
	if (!defined(fid)) {
		add_function(fid,value);
	}
}

void
DMesh::remove_function(string const fid)
throw() {
	// TODO: do we have to clean up here?
	mf.erase(fid);
}

void
DMesh::remove_function_ifdef(string const fid)
throw() {
	if (defined(fid)) {
		remove_function(fid);
	}
}


bool
DMesh::defined(string fid)
const throw() {
	if (mf.find(fid) == mf.end()) {
		return false;
	} else {
		return true;
	}
}

vector<string>
DMesh::get_fids(void)
const throw() {
	vector<string> fids;
	map<string,array2d>::const_iterator i = mf.begin();
	for (; i != mf.end(); ++i) {
		fids.push_back(i->first);
	}
	return fids;
}

double
DMesh::get(string const fid, const int i, const int j)
const throw(MeshException) {
	if (!defined(fid)) { // there is no such function
		string msg="DMesh::get: no such mesh function: ";
		msg+=fid;
		throw MeshException(msg);
	}
	if (is_out_of_domain(i,j)) {
		ostringstream ss;
		ss << "DMesh::get: indices out of range: " << i << ", " << j;
		ss << " (fid=" << fid << ")";
		throw MeshException(ss.str());
	}
	return ((map<string,array2d>)mf)[fid](i,j);
}

void
DMesh::set(string const fid, const int i, const int j, const double value)
throw(MeshException) {
	if (!defined(fid)) { // there is no such function
		string msg="DMesh::set: no such mesh function: ";
		msg+=fid;
		throw MeshException(msg);
	}
	if (is_out_of_domain(i,j)) {
		ostringstream ss;
		ss << "DMesh::set: indices out of range: " << i << ", " << j;
		ss << " (fid=" << fid << ")";
		throw MeshException(ss.str());
	}
	mf[fid](i,j)=value;
}

const blitz::Array<double,2>
DMesh::operator[](string const fid)
const throw(MeshException) {
	if (!defined(fid)) { // there is no such function
		string msg="DMesh::operator[]: no such mesh function: ";
		msg+=fid;
		throw MeshException(msg);
	}
	return ((map<string,array2d>)mf)[fid];
}

double
DMesh::triangle_area(double const x1, double const y1, double const x2,
	double const y2, double const x3, double const y3)
const {
	pair<double,double> ab(x2-x1,y2-y1);
	pair<double,double> ac(x3-x1,y3-y1);
	double s=0.5*fabs(vecxvec(ab,ac)); // half the vector product
	return s;
}

bool
DMesh::is_corner(const int i, const int j)
const throw() {
	// if corner point
	if ( ((i==0) || (i==(get_xdim()-1)))
			&& ((j==0) || (j==(get_ydim()-1))) ) {
		return true;
	} else {
		return false;
	}
}

double
DMesh::area(const int i, const int j)
const throw() {
	if (is_corner(i,j)) {
		return 0;
	}
	if (i == 0) {
		double sr=triangle_area(x(i,j-1),y(i,j-1),x(i+1,j),y(i+1,j),
			x(i,j+1),y(i,j+1));
		return 0.5*sr;
	} else if (j == 0) {
		double su=triangle_area(x(i-1,j),y(i-1,j),x(i,j+1),y(i,j+1),
			x(i+1,j),y(i+1,j));
		return 0.5*su;
	} else if (i == (get_xdim()-1)) {
		double sl=triangle_area(x(i,j-1),y(i,j-1),x(i-1,j),y(i-1,j),
			x(i,j+1),y(i,j+1));
		return 0.5*sl;
	} else if (j == (get_ydim()-1)) {
		double sl=triangle_area(x(i-1,j),y(i-1,j),x(i,j-1),y(i,j-1),
			x(i+1,j),y(i+1,j));
		return 0.5*sl;
	} else { // inner point
		double su=triangle_area(x(i-1,j),y(i-1,j),x(i,j+1),y(i,j+1),
			x(i+1,j),y(i+1,j));
		double sl=triangle_area(x(i-1,j),y(i-1,j),x(i,j-1),y(i,j-1),
			x(i+1,j),y(i+1,j));
		return 0.5*(su+sl);
	}
}

vector<string>
DMesh::get_attrs(void) const throw() {
	vector<string> names;
	map<string,double>::const_iterator i = attr.begin();
	for (; i != attr.end(); ++i) {
		names.push_back(i->first);
	}
	return names;
}

void
DMesh::set_attr(string const attribute, double const value)
throw(MeshException) {
	attr[attribute]=value;
}

double
DMesh::get_attr(string const attribute)
const throw(MeshException) {
	map<string,double>::const_iterator i;
	i=attr.find(attribute);
	if (i == attr.end()) {
		ostringstream ss;
		ss << "DMesh::get_attr: attribute `" << attribute
			<< "' not defined";
		throw MeshException(ss.str());
	} else {
		return i->second;
	}
}

bool
DMesh::attr_defined(string const attribute)
const {
	if (attr.find(attribute) == attr.end()) {
		return false;
	} else {
		return true;
	}
}

#ifdef HAVE_LIBHDF5
void
DMesh::save(string const filename)
const throw(MeshException) {
	try {
		H5::H5File file(filename, H5F_ACC_TRUNC);
		H5::Group savegroup=file.createGroup("mesh");
		H5::FloatType doubletype(H5::PredType::NATIVE_DOUBLE);
		H5::IntType inttype(H5::PredType::NATIVE_INT);
		// save 2D arrays (m_x, m_y and all mesh functions)
		hsize_t dims[2];
		dims[0]=get_xdim();
		dims[1]=get_ydim();
		H5::DataSpace ds2d(2, dims, NULL);
		H5::DataSet dataset;
		dataset = savegroup.createDataSet("x", doubletype, ds2d);
		dataset.write(m_x.data(), H5::PredType::NATIVE_DOUBLE);
		dataset = savegroup.createDataSet("y", doubletype, ds2d);
		dataset.write(m_y.data(), H5::PredType::NATIVE_DOUBLE);
		map<string,array2d>::const_iterator i = mf.begin();
		for (; i != mf.end(); ++i) {
			if ((i->first == "phi") || (i->first == "psi")
				|| (i->first == "vx") || (i->first == "vy")
				|| (i->first == "c")
				|| (i->first=="phi_t") || (i->first=="phi_h")
				|| (i->first=="phi_growth")
				|| (i->first=="phi_t_growth")
				|| (i->first=="phi_h_growth")) {
				dataset = savegroup.createDataSet(i->first,
						doubletype, ds2d);
				dataset.write(i->second.data(),
						H5::PredType::NATIVE_DOUBLE);
			}
		}
		// save scalar data
		H5::DataSpace ds0;
		savegroup.createAttribute("time", doubletype, ds0)
			.write(H5::PredType::NATIVE_DOUBLE, &time);
		savegroup.createAttribute("xdim", inttype, ds0)
			.write(H5::PredType::NATIVE_INT, &xdim);
		savegroup.createAttribute("ydim", inttype, ds0)
			.write(H5::PredType::NATIVE_INT, &ydim);
		savegroup.createAttribute("xmin", doubletype, ds0)
			.write(H5::PredType::NATIVE_DOUBLE, &xmin);
		savegroup.createAttribute("xmax", doubletype, ds0)
			.write(H5::PredType::NATIVE_DOUBLE, &xmax);
		savegroup.createAttribute("ymin", doubletype, ds0)
			.write(H5::PredType::NATIVE_DOUBLE, &ymin);
		savegroup.createAttribute("ymax", doubletype, ds0)
			.write(H5::PredType::NATIVE_DOUBLE, &ymax);
		// save all the dynamic attributes
		map<string,double>::const_iterator ai = attr.begin();
		for (; ai != attr.end(); ++ai) {
			string name=ai->first;
			double value=ai->second;
			savegroup.createAttribute(name, doubletype, ds0)
				.write(H5::PredType::NATIVE_DOUBLE, &value);
		}
	}
	catch (H5::Exception h5err) {
		throw MeshException(h5err.getDetailMsg());
	}
}

void
DMesh::load_dataset_2d(H5::Group const g, string const& setname,
	blitz::Array<double,2>& array) {
	try {
	hsize_t dims[2];
	H5::DataSet set=g.openDataSet(setname);
	H5::DataSpace space=set.getSpace();
	if (space.getSimpleExtentNdims() != 2) {
		string errmsg=
			"DMesh::load_dataset_2d: dataset `FOO' is not 2D";
		errmsg.replace(errmsg.find("FOO",0,3),3,setname);
		throw MeshException(errmsg);
	}
	space.getSimpleExtentDims(dims, NULL);
	if ((dims[0] != (hsize_t)xdim) || (dims[1] != (hsize_t)ydim)) {
		string errmsg="DMesh::load_dataset_2d: ";
		ostringstream sout(errmsg);
		sout << "datset `" << setname << "' is "
			<< dims[0] << "x" << dims[1] <<
			"while mesh is " << xdim << "x" << ydim ;
		errmsg=sout.str();
		throw MeshException(errmsg);
	}
	H5::DataType settype=set.getDataType();
	if (!(settype.detectClass(H5T_FLOAT) )) {
		string errmsg="DMesh::load_dataset_2d: ";
		ostringstream sout(errmsg);
		sout << "datset `" << setname
			<< "' is not a floating point type: "
			<< settype.getClass();
		errmsg=sout.str();
		throw MeshException(errmsg);
	}
	double* buffer=new double[xdim*ydim];
	set.read(buffer, settype);
	array.resize(xdim,ydim);
	memcpy(array.data(),buffer,sizeof(double)*xdim*ydim);
	delete[] buffer;
	buffer=(double*)0;
	} catch (H5::Exception h5err) {
		ostringstream ss;
		ss << "DMesh::load_dataset_2d: "
			<< "setname=" << setname << ": HDF5 exception: "
			<< h5err.getDetailMsg();
		throw MeshException(ss.str());
	}
}

void
DMesh::load_dataset_1d(H5::Group const g, string const& setname,
	blitz::Array<int,1>& array) {
	try {
	hsize_t size;
	H5::DataSet set=g.openDataSet(setname);
	H5::DataSpace space=set.getSpace();
	if (space.getSimpleExtentNdims() != 1) {
		string errmsg=
			"DMesh::load_dataset_1d: dataset `FOO' is not 1D";
		errmsg.replace(errmsg.find("FOO",0,3),3,setname);
		throw MeshException(errmsg);
	}
	space.getSimpleExtentDims(&size, NULL);
	H5::DataType settype=set.getDataType();
	if (!(settype.detectClass(H5T_INTEGER) )) {
		string errmsg="DMesh::load_dataset_1d: ";
		ostringstream sout(errmsg);
		sout << "datset `" << setname
			<< "' is not an integer type: " << settype.getClass();
		errmsg=sout.str();
		throw MeshException(errmsg);
	}
	int* buffer=new int[size];
	set.read(buffer, settype);
	array.resize(size);
	for (int i=0; (unsigned)i<size; ++i) {
		array(i)=*(buffer+i);
	}
	delete[] buffer;
	buffer=(int*)0;
	} catch (H5::Exception h5err) {
		ostringstream ss;
		ss << "DMesh::load_dataset_2d: "
			<< "setname=" << setname << ": HDF5 exception: "
			<< h5err.getDetailMsg();
		throw MeshException(ss.str());
	}
}

void
DMesh::load(string const filename)
throw(MeshFileException, MeshException) {
	try {
		H5::Exception::dontPrint();
		H5::FloatType doubletype(H5::PredType::NATIVE_DOUBLE);
		H5::IntType inttype(H5::PredType::NATIVE_INT);
		H5::H5File file(filename, H5F_ACC_RDONLY);
		H5::Group savegroup=file.openGroup("mesh");
		savegroup.openAttribute("time")
			.read(H5::PredType::NATIVE_DOUBLE, &time);
		savegroup.openAttribute("xdim")
			.read(H5::PredType::NATIVE_INT, &xdim);
		savegroup.openAttribute("ydim")
			.read(H5::PredType::NATIVE_INT, &ydim);
		savegroup.openAttribute("xmin")
			.read(H5::PredType::NATIVE_DOUBLE, &xmin);
		savegroup.openAttribute("xmax")
			.read(H5::PredType::NATIVE_DOUBLE, &xmax);
		savegroup.openAttribute("ymin")
			.read(H5::PredType::NATIVE_DOUBLE, &ymin);
		savegroup.openAttribute("ymax")
			.read(H5::PredType::NATIVE_DOUBLE, &ymax);
		// load essential dynamic attributes
		double attr_value;
		savegroup.openAttribute("phi0")
			.read(H5::PredType::NATIVE_DOUBLE, &attr_value);
		set_attr("phi0",attr_value);
		savegroup.openAttribute("cell_motility")
			.read(H5::PredType::NATIVE_DOUBLE, &attr_value);
		set_attr("cell_motility",attr_value);
		savegroup.openAttribute("upkeep_per_cell")
			.read(H5::PredType::NATIVE_DOUBLE, &attr_value);
		set_attr("upkeep_per_cell",attr_value);
		savegroup.openAttribute("death_rate")
			.read(H5::PredType::NATIVE_DOUBLE, &attr_value);
		set_attr("death_rate",attr_value);
		savegroup.openAttribute("o2_uptake")
			.read(H5::PredType::NATIVE_DOUBLE, &attr_value);
		set_attr("o2_uptake",attr_value);
		savegroup.openAttribute("tk1")
			.read(H5::PredType::NATIVE_DOUBLE, &attr_value);
		set_attr("tk1",attr_value);
		savegroup.openAttribute("ts1")
			.read(H5::PredType::NATIVE_DOUBLE, &attr_value);
		set_attr("ts1",attr_value);
		savegroup.openAttribute("hk1")
			.read(H5::PredType::NATIVE_DOUBLE, &attr_value);
		set_attr("hk1",attr_value);
		savegroup.openAttribute("hs1")
			.read(H5::PredType::NATIVE_DOUBLE, &attr_value);
		set_attr("hs1",attr_value);
		// load mesh geometry
		resize(xdim,ydim);
		load_dataset_2d(savegroup,"x",m_x);
		load_dataset_2d(savegroup,"y",m_y);
		// load essential mesh functions
		if (!defined("phi")) {
		 	add_function("phi");
		}
		load_dataset_2d(savegroup,"phi",mf["phi"]);
		if (!defined("psi")) {
		 	add_function("psi");
		}
		load_dataset_2d(savegroup,"psi",mf["psi"]);
		if (!defined("c")) {
		 	add_function("c");
		}
		load_dataset_2d(savegroup,"c",mf["c"]);
		if (!defined("vx")) {
		 	add_function("vx");
		}
		load_dataset_2d(savegroup,"vx",mf["vx"]);
		if (!defined("vy")) {
		 	add_function("vy");
		}
		load_dataset_2d(savegroup,"vy",mf["vy"]);
		if (!defined("phi_t")) {
		 	add_function("phi_t");
		}
		load_dataset_2d(savegroup,"phi_t",mf["phi_t"]);
		if (!defined("phi_h")) {
		 	add_function("phi_h");
		}
		load_dataset_2d(savegroup,"phi_h",mf["phi_h"]);
		if (!defined("phi_growth")) {
		 	add_function("phi_growth");
		}
		load_dataset_2d(savegroup,"phi_growth",mf["phi_growth"]);
	}
	catch (H5::FileIException filerr) {
		ostringstream ss;
		ss << "DMesh::load: cannot open file `"
			<< filename << "', HDF5 details: ";
		ss << filerr.getDetailMsg();
		throw MeshFileException(ss.str());
	}
	catch (H5::Exception h5err) {
		ostringstream ss;
		ss << "DMesh::load: HDF5 exception: " << h5err.getDetailMsg();
		throw MeshException(ss.str());
	}
}
#endif // ifdef HAVE_LIBHDF5

double
DMesh::vecxvec(const pair<double,double>& v1,
	const pair<double,double>& v2)
const {
	return v1.second*v2.first-v1.first*v2.second;
}

bool
DMesh::is_in_triangle (const DMesh::Point& p, const DMesh::Triangle& t)
const {
	// precision
	double err=1e-20;
	// triangle vertices
	double x1 = x(t.a.first,t.a.second);
	double x2 = x(t.b.first,t.b.second);
	double x3 = x(t.c.first,t.c.second);
	double y1 = y(t.a.first,t.a.second);
	double y2 = y(t.b.first,t.b.second);
	double y3 = y(t.c.first,t.c.second);
	// test point
	double x = p.first;
	double y = p.second;
	// triangle sides as vectors
	pair<double,double> ab(x2-x1,y2-y1);
	pair<double,double> bc(x3-x2,y3-y2);
	pair<double,double> ca(x1-x3,y1-y3);
	// vectors from all the vertices to the test point
	pair<double,double> ap(x-x1,y-y1);
	pair<double,double> bp(x-x2,y-y2);
	pair<double,double> cp(x-x3,y-y3);
	if ((vecxvec(ap,ab)*vecxvec(bp,bc)>=-err) &&
		(vecxvec(bp,bc)*vecxvec(cp,ca)>=-err)) {
		return true; // p inside t
	} else {
		return false;
	}
}

double
DMesh::interpolate_in_triangle(string const fid, const DMesh::Point& p,
	const DMesh::Triangle& t)
const throw(MeshException) {
	/*
	 * looking for a f(x,y)=a*x+b*y+c, solving this system first
	 *
	 * a*x1+b*y1+c = f1   |        | x1  y1  1 |        | f1 |
	 * a*x2+b*y2+c = f2   | => A = | x2  y2  1 |, rhs = | f2 |
	 * a*x3+b*y3+c = f3   |        | x3  y3  1 |        | f3 |
	 *
	 * a, b, c are unknowns
	 */
	double x1 = x(t.a.first,t.a.second);
	double x2 = x(t.b.first,t.b.second);
	double x3 = x(t.c.first,t.c.second);
	double y1 = y(t.a.first,t.a.second);
	double y2 = y(t.b.first,t.b.second);
	double y3 = y(t.c.first,t.c.second);
	double f1 = get(fid,t.a.first,t.a.second);
	double f2 = get(fid,t.b.first,t.b.second);
	double f3 = get(fid,t.c.first,t.c.second);
	// det(A)
	double D = x1*(y2-y3) - x2*(y1-y3) + x3*(y1-y2);
	if (fabs(D) < 1.0e-20) { // degenerate triangle
		ostringstream ss;
		ss << "DMesh::interpolate_in_triangle: "
			<< "det(A)=" << D << " vertices:"
			<< t.a.first << "," << t.a.second << " "
			<< t.b.first << "," << t.b.second << " "
			<< t.c.first << "," << t.c.second;
		throw MeshException(ss.str());
	}
	// det(A:s/x/f/))
	double D1 = f1*(y2-y3) - f2*(y1-y3) + f3*(y1-y2);
	// det(A:s/y/f/))
	double D2 = x1*(f2-f3) - x2*(f1-f3) + x3*(f1-f2);
	// det(A:s/\<1\>/f./))
	double D3 = x1*(y2*f3-y3*f2) - x2*(y1*f3-y3*f1) + x3*(y1*f2-y2*f1);
	double a=D1/D;
	double b=D2/D;
	double c=D3/D;
	double value=a*p.first+b*p.second+c;
	return value;
}

double
DMesh::interpolate(string const fid, const double x, const double y)
const throw(MeshException) {
	// for each quandrangle of the grid do
	for (int i=0; i<(xdim-1); i++) {
		for (int j=0; j<(ydim-1); j++) {
			/*   Triangulation first:
			 *
			 *   B (i,j+1)  ----- C (i+1,j+1)
			 *      |                 |
			 *      |                 |
			 *   A (i,j) -------- D (i+1,j)
			 *
			 *   Possible triangulations: ABC+ACD and ABD+BCD.
			 *
			 *   If A and C are in different semiplanes with
			 *   respect to BD, ABD+BCD is valid triangulation.
			 */
			pair<double,double>
				BA(this->x(i,j)-this->x(i,j+1),
					this->y(i,j)-this->y(i,j+1));
			pair<double,double>
				BD(this->x(i+1,j)-this->x(i,j+1),
					this->y(i+1,j)-this->y(i,j+1));
			pair<double,double>
				BC(this->x(i+1,j+1)-this->x(i,j+1),
					this->y(i+1,j+1)-this->y(i,j+1));
			Triangle t1, t2;
			// TODO: check that triangulation respects inner-outer
			if (vecxvec(BA,BD)*vecxvec(BC,BD)<0) {
				// A, C in different semiplanes, using ABD+BCD
				Triangle thet1(MP(i,j),MP(i,j+1),MP(i+1,j));
				Triangle thet2(MP(i,j+1),MP(i+1,j+1),MP(i+1,j));
				t1=thet1;
				t2=thet2;
			} else { // ABC+ACD
				Triangle thet1(MP(i,j),MP(i,j+1),MP(i+1,j+1));
				Triangle thet2(MP(i,j),MP(i+1,j+1),MP(i+1,j));
				t1=thet1;
				t2=thet2;
			}
			// if (x,y) in triangle t1 or triangle t2, interpolate
			Point p(x,y);
			if (is_in_triangle(p,t1)) {
				return interpolate_in_triangle(fid,p,t1);
			} else if (is_in_triangle(p,t2)) {
				return interpolate_in_triangle(fid,p,t2);
			}
			// else means not in this quad, proceed further
		}
	}
	return 0.0; // not defined, really, for outer points; but no error
}


