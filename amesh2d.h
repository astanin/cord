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

#ifndef AMESH2D_H
#define AMESH2D_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <vector>

#include <blitz/array.h>
using std::string;
using std::vector;
using std::exception;

/** @brief exception class that describes the error of AMesh2D classes */
class MeshException : public exception {
protected:
	string msg; ///< exception message
public:
	MeshException(string const errormsg="unspecified error")
	throw() : msg(errormsg) {}
	~MeshException() throw() {};
	virtual const char* what() const throw() { return msg.c_str(); }
};

/** @brief exception class for file operation exceptions */
class MeshFileException : public MeshException {
public:
	MeshFileException(string const errormsg="unspecified file error")
	throw() : MeshException(errormsg) {}
	~MeshFileException() throw() {};
};

/// Meshes whose state depends on time should implement this interface
class AEvolvableInterface {
public:
	virtual ~AEvolvableInterface() {}
	virtual double get_time(void) const = 0;
	virtual void set_time(double const _time) = 0;
	// increment time variable
	virtual void inc_time(double const _dt) = 0;
};

/// Meshes which may be saved to disk or loaded from disk should implement
/// APermanentInterface
class APermanentInterface {
public:
	virtual ~APermanentInterface() {}
	virtual void
	save(string const filename)
	const throw(MeshException) = 0;
	virtual void
	load(string const filename)
	throw(MeshFileException,MeshException) = 0;
};

typedef blitz::Array<double,2> array2d; ///< 2D array of floats type

/** @brief abstract 2D rectangular mesh class
 * Mesh is a mapping of (i,j) points in index space [0:xdim-1],[0:ydim-1] to
 * their (x,y) physical coordinates.
 * Rectangularity means that for a given mesh point (i0,j0) their neighbours
 * are any mesh points (i,j) for which we have abs(i0-i)+abs(j0-j)==1;
 * Additionally any number of mesh functions may be defined on the mesh points.
 * Mesh function are double-valued and are referenced with std::string names.
 * Mesh functions are distinguished by function ids (e.g. their literal names).
 * */
template<class fid_t>
class AMesh2D : public AEvolvableInterface
#ifdef HAVE_LIBHDF5
	, public APermanentInterface
#endif
{
public:
	virtual ~AMesh2D() {}
	virtual AMesh2D<fid_t>* clone() const = 0;

	/// return physical x for mesh point (i,j)
	virtual double x(const int i, const int j) const = 0;
	/// return physical y for mesh point (i,j)
	virtual double y(const int i, const int j) const = 0;
	/// set physical x for mesh point (i,j)
	virtual void set_x(const int i, const int j, double const value) = 0;
	/// set physical y for mesh point (i,j)
	virtual void set_y(const int i, const int j, double const value) = 0;
	/// the number of vertices along x-axis
	virtual int get_xdim(void) const = 0;
	/// the number of vertices along y-axis
	virtual int get_ydim(void) const = 0;
	/// lower x bound of the outbound rectangle
	virtual double get_xmin(void) const = 0;
	/// upper x bound of the outbound rectangle
	virtual double get_xmax(void) const = 0;
	/// lower y bound of the outbound rectangle
	virtual double get_ymin(void) const = 0;
	/// upper y bound of the outbound rectangle
	virtual double get_ymax(void) const = 0;
	/// average x step; WARNING: assuming uniform rectangular grid
	virtual double get_dx(void) const {
		return (x(get_xdim()-1,0)-x(0,0))/(get_xdim()-1);
	}
	/// average y step; WARNING: assuming uniform rectangular grid
	virtual double get_dy(void) const {
		return (y(0,get_ydim()-1)-y(0,0))/(get_ydim()-1);
	}
	/// area corresponding to mesh point (i,j)
	virtual double area(int const i, int const j) const = 0;


	/// delete all mesh functions and attributes associated with the mesh
	virtual void clear(void) = 0;

	/** @brief add double-valued mesh function
	 * @param fid[in]	name (ID) of the function
	 * @param value[in]	initial value */
	virtual void
	add_function(fid_t const fid, double const value=0.0) = 0;

	/** @brief add double-valued mesh function if does not exist yet */
	virtual void
	add_function_ifndef(fid_t const fid, double const value=0.0) = 0;

	/** @brief remove mesh function from the mesh
	 * @param fid[in]	name (ID) of the function */
	virtual void remove_function(fid_t const fid) = 0;

	/** @brief remove mesh function from the mesh if one exists */
	virtual void
	remove_function_ifdef(fid_t const fid) = 0;

	/** @brief return vector of names of attached mesh functions */
	virtual vector<fid_t> get_fids(void) const = 0;

	/** @brief returns true if mesh function @c fid is defined, false
	 * otherwise */
	virtual bool defined(fid_t const fid) const = 0;

	/** @brief get mesh function value
	 * @param fid[in]	name (id) of the function
	 * @param i[in]		x-axis index
	 * @param j[in]		y-axis index
	 * @return		value of @c fid(i,j) */
	virtual double
	get(fid_t const fid, const int i, const int j) const = 0;

	/** @brief set mesh function value
	 * @param fid[in]	name (id) of the function
	 * @param i[in]		x-axis index
	 * @param j[in]		y-axis index
	 * @param value[in]	new value of the @c fid(i,j) */
	virtual void
	set(fid_t const fid, const int i, const int j, const double value)
	= 0;

	/** @brief return reference to array corresponding to mesh function
	 * @c fid
	 *
	 * operator[] provides more math-like access to mesh functions and
	 * allows to read and set them with @c mesh["fid"](i,j) syntax,
	 * that may be handy in numerical code.
	 *
	 * Warning: this way to access mesh functions breaks incapsulation,
	 * but increased readability of numerical code looks more important.
	 *
	 * @param fid[in]	name (id) of the function
	 * @return		array2d with function values */
	virtual array2d const operator[](fid_t const fid) const = 0;

	/// interpolate mesh function in a given point
	virtual double
	interpolate (fid_t const fid, double const x, double const y)
	const throw(MeshException) = 0;

	/** @brief return vector of names of attached attributes */
	virtual vector<string> get_attrs(void) const = 0;

	/** @brief set attribute associated with mesh */
	virtual void
	set_attr(string const attribute, double const value)
	throw(MeshException) = 0;

	/** @brief get value of attributes associated with mesh */
	virtual double
	get_attr(string const attribute) const throw(MeshException) = 0;

	/** @brief check if attribute @c attribute is associated with mesh */
	virtual bool
	attr_defined(string const attribute) const = 0;

	/** @brief return true if either index i or j is out of range */
	virtual int is_out_of_domain(const int i, const int j) const {
		if ((i >= get_xdim()) || (i < 0)
			|| (j >= get_ydim()) || (j < 0)) {
			return true;
		} else {
			return false;
		}
	}
	/** @brief check if (@c i,@c j) belongs to domain (is either inner or
	 * boundary point) */
	virtual bool is_inner(const int i, const int j) const {
		return (!is_out_of_domain(i,j));
	}
	/** @brief check if (@c i,@c j) lies on the domain boundary */
	virtual bool is_border(const int i, const int j) const {
		if ((i == 0) || (i == (get_xdim()-1))
			|| (j == 0) || (j == (get_ydim()-1))) {
			return true;
		} else {
			return false;
		}
	}
};

#endif /* AMESH2D_H */
