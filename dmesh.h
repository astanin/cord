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

#ifndef DMESH_H
#define DMESH_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>

using std::string;
using std::pair;
using std::map;
using std::vector;
using std::ofstream;

#include "amesh2d.h"
#include "dmeshops.h"

#ifdef HAVE_LIBHDF5
// saving to HDF5 file
#include <H5Cpp.h>
#endif

#include "function.h"

template<class fid_t> class DMesh;

template<class fid_t>
void
dump2dx_scalar_field (DMesh<fid_t> const& m, fid_t const& fid, ofstream& fs);

template<class fid_t>
void
dump2dx_vector_field(DMesh<fid_t> const& m, fid_t const& vxfid,
	fid_t const& vyfid, string const& fieldname, ofstream& fs);

template<class fid_t>
void
dump2dx(DMesh<fid_t> const& m, string const& filename);

template<class fid_t>
string id2str(fid_t& id);

/** @brief 2D mesh with rectangular-like connections and named functions */
template<class fid_t>
class DMesh : public AMesh2D<fid_t>
{
	typedef map<fid_t, array2d > pile_of_arrays;

private:

	double time; ///< value of time variable corresponding to mesh state
	int xdim; ///< mesh dimensionality along x-axis (the number of vertices)
	int ydim; ///< mesh dimensionality along y-axis (the number of vertices)
	double xmin; ///< lower x bound of the outbound rectangle
	double xmax; ///< upper x bound of the outbound rectangle
	double ymin; ///< lower y bound of the outbound rectangle
	double ymax; ///< upper y bound of the outbound rectangle
	array2d m_x; ///< physical x-coordinates of vertex positions (i,j)
	array2d m_y; ///< physical y-coordinates of vertex positions (i,j)

	/** set of named functions defined on the mesh (mesh functions) */
	pile_of_arrays mf;

// interpolation infrastructure

	/// represents meshpoint's coords in index-space
	typedef pair<int,int> MP;

	/// represents physical coords (x,y) of a point
	typedef pair<double,double> Point;

	/// represents triangle vertices
	class Triangle {
	public:
		MP a;
		MP b;
		MP c;
		Triangle(): a(), b(), c() {}
		Triangle(const MP& _a, const MP& _b, const MP& _c) :
			a(_a), b(_b), c(_c) {}
	};

	/// returns (signed) vector product of vectors @c v1 and @c v2
	double
	vecxvec(const pair<double,double>& v1,
			const pair<double,double>& v2) const;

	/// returns true if @c p is in in the triangle @t
	bool is_in_triangle(const Point& p, const Triangle& t) const;

	/// interpolation procedure within a triangle
	double
	interpolate_in_triangle(fid_t const fid,
		const Point& p, const Triangle& t) const throw(MeshException);

	/// triangle area
	double triangle_area(double const x1, double const y1, double const x2,
		double const y2, double const x3, double const y3) const;

	/** @brief return non-zero if (i,j) is mesh corner */
	bool is_corner(const int i, const int j) const throw();

protected:

	/** uniformly map mesh vertices to [xmin;xmax]*[ymin;ymax] rectangle */
	void
	init_mesh(void);

	/** resize mesh according to new dimensions, reset all mesh functions */
	void
	resize(int const _xdim, int const _ydim);

#ifdef HAVE_LIBHDF5
	/** load 2D dataset into blitz::Array<double,2> */
	void
	load_dataset_2d(H5::Group const g, string const& setname,
		blitz::Array<double,2>& array);

	/** load 1D dataset into blitz::Array<double,1> */
	void
	load_dataset_1d(H5::Group const g, string const& setname,
		blitz::Array<int,1>& array);
#endif

	/** set of named mesh attributes and parameters */
	map<string, double> attr;

public:

	virtual ~DMesh() {}

	virtual DMesh<fid_t>* clone(void) const {
		return new DMesh<fid_t>(*this);
	}

	/** constructs nx*ny mesh and maps it to [_xmin;_xmax]*[_ymin;_ymax]
	 * rectangle */
	DMesh(int _xdim=20, int _ydim=20, double _xmin=0.0, double _xmax=1.0,
		double _ymin=0.0, double _ymax=1.0);

	/** copy constructor */
	DMesh(const DMesh<fid_t>& om);

	/** assignment operator */
	DMesh& operator=(const DMesh<fid_t>& rhs);

	/// delete all mesh functions and mesh attributes
	virtual void
	clear(void) {
		mf.clear();
		attr.clear();
	}

// mesh geometry

	virtual double
		x(const int i, const int j) const throw(MeshException);

	virtual double
		y(const int i, const int j) const throw(MeshException);

	virtual void
		set_x(const int i, const int j, double const value)
		throw(MeshException);

	virtual void
		set_y(const int i, const int j, double const value)
		throw(MeshException);

	virtual int
		get_xdim(void) const {return xdim;}

	virtual int
		get_ydim(void) const {return ydim;}

	virtual double
		get_xmin(void) const {return xmin;}

	virtual double
		get_xmax(void) const {return xmax;}

	virtual double
		get_ymin(void) const {return ymin;}

	virtual double
		get_ymax(void) const {return ymax;}


	/** @brief area corresponding to mesh point (i,j) */
	virtual double area(const int i, const int j) const throw();

	// mesh evolution variable (time)

	virtual double get_time(void) const {return time;}

	virtual void set_time(double const _time) throw() { time=_time; }

	virtual void inc_time(double const _dt) throw() { time+=_dt; }

	// mesh functions

	virtual void
	add_function(fid_t const fid, double const value=0.0)
	throw();

	virtual void
	add_function_ifndef(fid_t const fid, double const value=0.0)
	throw();

	virtual void
	remove_function(fid_t fid) throw();

	virtual void
	remove_function_ifdef(fid_t const fid) throw();

	virtual vector<fid_t>
	get_fids(void) const throw();

	virtual bool
	defined(fid_t fid) const throw();

	virtual double
	get(fid_t const fid, const int i, const int j)
	const throw(MeshException);

	virtual void
	set(fid_t const fid, const int i, const int j,
		const double value) throw(MeshException);

	virtual const blitz::Array<double,2>
	operator[](fid_t const fid) const throw(MeshException);


// mesh attributes
	virtual vector<string>
	get_attrs(void) const throw();

	virtual void
	set_attr(string const attribute,
		double const value) throw(MeshException);

	virtual double
	get_attr(string const attribute) const throw(MeshException);

	virtual bool
	attr_defined(string const attribute) const;

// regridding and interpolation

	/// interpolate mesh function in a given point
	virtual double
	interpolate (fid_t const fid, const double x, const double y)
	const throw(MeshException);

// file operations

#ifdef HAVE_LIBHDF5
	/** @brief save mesh state to the file */
	virtual void
	save(string const filename) const throw(MeshException);

	/** @brief restore mesh state from the file */
	virtual void
	load(string const filename) throw(MeshFileException,MeshException);
#endif

// friends

	friend void
	dump2dx_scalar_field<>
	(DMesh<fid_t> const& m, fid_t const& fid, ofstream& fs);

	friend void
	dump2dx_vector_field<>
	(DMesh<fid_t> const& m, fid_t const& vxfid, fid_t const& vyfid,
		string const& fieldname, ofstream& fs);

	friend void
	dump2dx<>(DMesh<fid_t> const& m, string const& filename);

};

#endif /* DMESH_H */

