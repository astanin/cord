/* Tumour cord growth model */

/*
 * Copyright (C) 2006 Sergey Astanin, Luigi Preziosi, Andrea Tosin
 * * This program is free software; you can redistribute it and/or modify
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

#ifndef NEWTON_H
#define NEWTON_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <iostream>
#include "function.h"

using std::cerr;

class ARootFinder {
public:
	virtual ~ARootFinder() {};
	virtual void iterate(void) throw(FunctionException) = 0;
	/// @return current root estimation
	virtual double get_root(void) const = 0; 
	/// stop condition for iterative process
	virtual bool is_found(void) const = 0; 
};

class NewtonRootFinder : public ARootFinder {
private:
	ADoubleFunction const& f; ///< f(x), we look for 
	ADoubleFunction const& df; ///< derivative of f(x)
	double x0; ///< current root estimation
	double absresidual; ///< absolute residual |f(x0)| to be achieved
public:
	NewtonRootFinder(ADoubleFunction const& _f, ADoubleFunction const& _df,
		double const _x0, double const _residual) :
		f(_f), df(_df), x0(_x0), absresidual(_residual) {}
	virtual void iterate(void) throw(FunctionException) {
		double vdf=df(x0);
		if (fabs(vdf)<1e-99) {
			throw FunctionException("NewtonRootFinder: derivative "
			"is almost zero, bad initial root approximation");
		}
		double vf=f(x0);
		double x1=x0-vf/vdf;
		x0=x1;
	}
	virtual double get_root(void) const {
		return x0;
	}
	virtual bool is_found(void) const {
		double residual=fabs(f(x0));
		if (residual < absresidual) {
			return true;
		} else {
			return false;
		}
	}
};

#endif



