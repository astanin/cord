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

template<class fid_t>
string id2str(fid_t& id) {
	if (typeid(id) == typeid(int)) {
		switch (id) {
		case PHI:
			return "phi";
			break;
		case PHI1:
			return "phi1";
			break;
		case PHI2:
			return "phi2";
			break;
		case O2:
			return "c";
			break;
		case GLC:
			return "c_g";
			break;
		case PSI:
			return "psi";
			break;
		case VX:
			return "vx";
			break;
		case VY:
			return "vy";
			break;
		case PHI_T:
			return "phi_t";
			break;
		case PHI_H:
			return "phi_h";
			break;
		case PHI_GROWTH:
			return "phi_growth";
			break;
		case TMP1:
			return "tmp1";
			break;
		case TMP2:
			return "tmp2";
			break;
		case VASC:
			return "vasc";
			break;
		default:
			std::ostringstream ss;
			ss << "fid_" << id;
			return ss.str();
			break;
		}
	} else {
		std::ostringstream ss;
		ss << id;
		return ss.str();
	}
}

template<class fid_t>
DMesh<fid_t>::DMesh
(int _xdim, int _ydim, double _xmin, double _xmax, double _ymin, double _ymax) :
	time(0), xdim(_xdim), ydim(_ydim), xmin(_xmin),
	xmax(_xmax), ymin(_ymin), ymax(_ymax), m_x(_xdim,_ydim),
	m_y(_xdim, _ydim), mf(), attr() {
		init_mesh();
}

template<class fid_t>
DMesh<fid_t>::DMesh(const DMesh<fid_t>& om) :
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
	typename pile_of_arrays::const_iterator i = om.mf.begin();
	for (; i != om.mf.end(); ++i) {
		fid_t mfkey=i->first;
		array2d mfc;
		mfc.resize(i->second.shape());
		mfc=i->second.copy();
		mf.insert(make_pair(mfkey,mfc));
	}
	attr=om.attr;
}

template<class fid_t>
DMesh<fid_t>&
DMesh<fid_t>::operator=(const DMesh<fid_t>& rhs) {
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
		typename pile_of_arrays::const_iterator i=rhs.mf.begin();
		for (; i != rhs.mf.end(); ++i) {
			fid_t mfkey=i->first;
			array2d mfc;
			mfc.resize(i->second.shape());
			mfc=i->second.copy();
			mf.insert(make_pair(mfkey,mfc));
		}
		attr=rhs.attr;
		return *this;
	}
}

template<class fid_t>
double
DMesh<fid_t>::x(const int i, const int j)
const throw(MeshException) {
	if (!this->is_out_of_domain(i,j)) {
		return m_x(i,j);
	} else {
		ostringstream ss;
		ss << "DMesh::x: indices out of range: " << i << ", " << j ;
		throw MeshException(ss.str());
	}
}

template<class fid_t>
double
DMesh<fid_t>::y(const int i, const int j)
const throw(MeshException) {
	if (!this->is_out_of_domain(i,j)) {
		return m_y(i,j);
	} else {
		ostringstream ss;
		ss << "DMesh::y: indices out of range: " << i << ", " << j ;
		throw MeshException(ss.str());
	}
}

template<class fid_t>
void
DMesh<fid_t>::set_x(const int i, const int j, double const value)
throw(MeshException) {
	if (!this->is_out_of_domain(i,j)) {
		m_x(i,j)=value;
	} else {
		ostringstream ss;
		ss << "DMesh::set_x: indices out of range: " << i << ", " << j ;
		throw MeshException(ss.str());
	}
}

template<class fid_t>
void
DMesh<fid_t>::set_y(const int i, const int j, double const value)
throw(MeshException) {
	if (!this->is_out_of_domain(i,j)) {
		m_y(i,j)=value;
	} else {
		ostringstream ss;
		ss << "DMesh::set_y: indices out of range: " << i << ", " << j ;
		throw MeshException(ss.str());
	}
}

template<class fid_t>
void
DMesh<fid_t>::init_mesh(void) {
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

template<class fid_t>
void
DMesh<fid_t>::resize(int const _xdim, int const _ydim) {
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

template<class fid_t>
void
DMesh<fid_t>::add_function(fid_t const fid, double const value)
throw() {
	array2d a;
	mf[fid]=a;
	mf[fid].resize(xdim,ydim);
	mf[fid]=value;
}

template<class fid_t>
void
DMesh<fid_t>::add_function_ifndef(fid_t const fid, double const value)
throw() {
	if (!defined(fid)) {
		add_function(fid,value);
	}
}

template<class fid_t>
void
DMesh<fid_t>::remove_function(fid_t const fid)
throw() {
	// TODO: do we have to clean up here?
	mf.erase(fid);
}

template<class fid_t>
void
DMesh<fid_t>::remove_function_ifdef(fid_t const fid)
throw() {
	if (defined(fid)) {
		remove_function(fid);
	}
}


template<class fid_t>
bool
DMesh<fid_t>::defined(fid_t fid)
const throw() {
	if (mf.find(fid) == mf.end()) {
		return false;
	} else {
		return true;
	}
}

template<class fid_t>
vector<fid_t>
DMesh<fid_t>::get_fids(void)
const throw() {
	vector<fid_t> fids;
	typename pile_of_arrays::const_iterator i = mf.begin();
	for (; i != mf.end(); ++i) {
		fids.push_back(i->first);
	}
	return fids;
}

template<class fid_t>
double
DMesh<fid_t>::get(fid_t const fid, const int i, const int j)
const throw(MeshException) {
	if (!defined(fid)) { // there is no such function
		string msg="DMesh::get: no such mesh function: ";
		msg+=id2str(fid);
		throw MeshException(msg);
	}
	if (this->is_out_of_domain(i,j)) {
		ostringstream ss;
		ss << "DMesh::get: indices out of range: " << i << ", " << j;
		ss << " (fid=" << id2str(fid) << ")";
		throw MeshException(ss.str());
	}
	return ((pile_of_arrays)mf)[fid](i,j);
}

template<class fid_t>
void
DMesh<fid_t>::set(fid_t const fid, const int i, const int j, const double value)
throw(MeshException) {
	if (!defined(fid)) { // there is no such function
		string msg="DMesh::set: no such mesh function: ";
		msg+=id2str(fid);
		throw MeshException(msg);
	}
	if (this->is_out_of_domain(i,j)) {
		ostringstream ss;
		ss << "DMesh::set: indices out of range: " << i << ", " << j;
		ss << " (fid=" << id2str(fid) << ")";
		throw MeshException(ss.str());
	}
	mf[fid](i,j)=value;
}

template<class fid_t>
const blitz::Array<double,2>
DMesh<fid_t>::operator[](fid_t const fid)
const throw(MeshException) {
	if (!defined(fid)) { // there is no such function
		string msg="DMesh::operator[]: no such mesh function: ";
		msg+=id2str(fid);
		throw MeshException(msg);
	}
	return ((pile_of_arrays)mf)[fid];
}

template<class fid_t>
double
DMesh<fid_t>::triangle_area(double const x1, double const y1, double const x2,
	double const y2, double const x3, double const y3)
const {
	pair<double,double> ab(x2-x1,y2-y1);
	pair<double,double> ac(x3-x1,y3-y1);
	double s=0.5*fabs(vecxvec(ab,ac)); // half the vector product
	return s;
}

template<class fid_t>
bool
DMesh<fid_t>::is_corner(const int i, const int j)
const throw() {
	// if corner point
	if ( ((i==0) || (i==(get_xdim()-1)))
			&& ((j==0) || (j==(get_ydim()-1))) ) {
		return true;
	} else {
		return false;
	}
}

template<class fid_t>
double
DMesh<fid_t>::area(const int i, const int j)
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

template<class fid_t>
vector<string>
DMesh<fid_t>::get_attrs(void) const throw() {
	vector<string> names;
	typename map<string,double>::const_iterator i = attr.begin();
	for (; i != attr.end(); ++i) {
		names.push_back(i->first);
	}
	return names;
}

template<class fid_t>
void
DMesh<fid_t>::set_attr(string const attribute, double const value)
throw(MeshException) {
	attr[attribute]=value;
}

template<class fid_t>
double
DMesh<fid_t>::get_attr(string const attribute)
const throw(MeshException) {
	typename map<string,double>::const_iterator i;
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

template<class fid_t>
bool
DMesh<fid_t>::attr_defined(string const attribute)
const {
	if (attr.find(attribute) == attr.end()) {
		return false;
	} else {
		return true;
	}
}

#ifdef HAVE_LIBHDF5
template<class fid_t>
void
DMesh<fid_t>::save(string const filename)
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
		typename pile_of_arrays::const_iterator i = mf.begin();
		for (; i != mf.end(); ++i) {
			dataset = savegroup.createDataSet(id2str(i->first),
					doubletype, ds2d);
			dataset.write(i->second.data(),
					H5::PredType::NATIVE_DOUBLE);
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
		typename map<string,double>::const_iterator ai = attr.begin();
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

template<class fid_t>
void
DMesh<fid_t>::load_dataset_2d(H5::Group const g, string const& setname,
	blitz::Array<double,2>& array) {
	try {
	hsize_t dims[2];
	H5::DataSet set;
	try {
		set=g.openDataSet(setname);
	} catch (H5::Exception h5err) { // cannot open data set
		array=0.0;
		return;
	}
	H5::DataSpace space=set.getSpace();
	if (space.getSimpleExtentNdims() != 2) {
		string errmsg="DMesh::load_dataset_2d: dataset `FOO' is not 2D";
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
	return;
	} catch (H5::Exception h5err) {
		ostringstream ss;
		ss << "DMesh::load_dataset_2d: "
			<< "setname=" << setname << ": HDF5 exception: "
			<< h5err.getDetailMsg();
		throw MeshException(ss.str());
	}
}

template<class fid_t>
void
DMesh<fid_t>::load_dataset_1d(H5::Group const g, string const& setname,
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

template<class fid_t>
void
DMesh<fid_t>::load(string const filename)
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
		savegroup.openAttribute("growth_rate")
			.read(H5::PredType::NATIVE_DOUBLE, &attr_value);
		set_attr("growth_rate",attr_value);
		savegroup.openAttribute("o2_uptake")
			.read(H5::PredType::NATIVE_DOUBLE, &attr_value);
		set_attr("o2_uptake",attr_value);
		savegroup.openAttribute("permability")
			.read(H5::PredType::NATIVE_DOUBLE, &attr_value);
		set_attr("permability",attr_value);
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
		savegroup.openAttribute("host_activity")
			.read(H5::PredType::NATIVE_DOUBLE, &attr_value);
		set_attr("host_activity",attr_value);
		savegroup.openAttribute("conversion_rate")
			.read(H5::PredType::NATIVE_DOUBLE, &attr_value);
		set_attr("conversion_rate",attr_value);
		savegroup.openAttribute("anaerobic_rate")
			.read(H5::PredType::NATIVE_DOUBLE, &attr_value);
		set_attr("anaerobic_rate",attr_value);
		savegroup.openAttribute("glc_diffusion")
			.read(H5::PredType::NATIVE_DOUBLE, &attr_value);
		set_attr("glc_diffusion",attr_value);
		// load mesh geometry
		resize(xdim,ydim);
		load_dataset_2d(savegroup,"x",m_x);
		load_dataset_2d(savegroup,"y",m_y);
		// load essential mesh functions
		if (!defined(PHI)) {
		 	add_function(PHI);
		}
		load_dataset_2d(savegroup,"phi",mf[PHI]);
		if (!defined(PHI1)) {
		 	add_function(PHI1);
		}
		load_dataset_2d(savegroup,"phi1",mf[PHI1]);
		if (!defined(PHI2)) {
		 	add_function(PHI2);
		}
		load_dataset_2d(savegroup,"phi2",mf[PHI2]);
		if (!defined(PSI)) {
		 	add_function(PSI);
		}
		load_dataset_2d(savegroup,"psi",mf[PSI]);
		if (!defined(O2)) {
		 	add_function(O2);
		}
		load_dataset_2d(savegroup,"c",mf[O2]);
		if (!defined(GLC)) {
		 	add_function(GLC);
		}
		load_dataset_2d(savegroup,"c_g",mf[GLC]);
		if (!defined(VX)) {
		 	add_function(VX);
		}
		load_dataset_2d(savegroup,"vx",mf[VX]);
		if (!defined(VY)) {
		 	add_function(VY);
		}
		load_dataset_2d(savegroup,"vy",mf[VY]);
		if (!defined(PHI_T)) {
		 	add_function(PHI_T);
		}
		load_dataset_2d(savegroup,"phi_t",mf[PHI_T]);
		if (!defined(PHI_H)) {
		 	add_function(PHI_H);
		}
		load_dataset_2d(savegroup,"phi_h",mf[PHI_H]);
//		if (!defined(PHI_GROWTH)) {
//		 	add_function(PHI_GROWTH);
//		}
//		load_dataset_2d(savegroup,"phi_growth",mf[PHI_GROWTH]);
		if (!defined(VASC)) {
		 	add_function(VASC);
		}
		load_dataset_2d(savegroup,"vasc",mf[VASC]);
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

template<class fid_t>
double
DMesh<fid_t>::vecxvec(const pair<double,double>& v1,
	const pair<double,double>& v2)
const {
	return v1.second*v2.first-v1.first*v2.second;
}

template<class fid_t>
bool DMesh<fid_t>::is_in_triangle
(const DMesh<fid_t>::Point& p, const DMesh<fid_t>::Triangle& t)
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

template<class fid_t>
double
DMesh<fid_t>::interpolate_in_triangle
(fid_t const fid, const DMesh<fid_t>::Point& p, const DMesh<fid_t>::Triangle& t)
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

template<class fid_t>
double
DMesh<fid_t>::interpolate(fid_t const fid, const double x, const double y)
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

// templates
template string id2str<int>(int& id);
template class DMesh<int>;
