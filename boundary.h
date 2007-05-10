/* Tumour cord growth model */

/*
 * Copyright (C) 2005-2007 Sergey Astanin, Luigi Preziosi, Andrea Tosin
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but

 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef BOUNDARY_H
#define BOUNDARY_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>

/* Boundary conditions classes */

/**
 * Boundary condition is represented as 
 * \[
 * 	a u + b \frac{\partial u}{\partial \mathbf{n}} = c
 * \]
 * where $\mathbf{n}$ is an _external_ normal to the boundary
 */
class BoundaryCondition {
private:
	double _a;
	double _b;
	double _c;
public:
	typedef enum {
		UNKNOWN_BC = 0,
		DIRICHLET_BC = 1,
		NEUMANN_BC = 2,
		ROBIN_BC = 3
	} Type;

	BoundaryCondition() : _a(0.0), _b(1.0), _c(0.0) {}
	BoundaryCondition(double const a, double const b, double const c) :
		_a(a), _b(b), _c(c) {}

	static BoundaryCondition
	createDirichletBC(double const val=0.0) {
		BoundaryCondition bc(1.0,0.0,val);
		return bc;
	}

	static BoundaryCondition
	createNeumannBC(double const val=0.0) {
		BoundaryCondition bc(0.0,1.0,val);
		return bc;
	}

	Type get_type() const {
		if ((fabs(_a)+fabs(_b)) < 1e-99) {
			return UNKNOWN_BC;
		} else if (fabs(_a) < 1e-99) {
			return NEUMANN_BC;
		} else if (fabs(_b) < 1e-99) {
			return DIRICHLET_BC;
		} else {
			return ROBIN_BC;
		}
	}
	void set_a(double const a) { this->_a = a; }
	void set_b(double const b) { this->_b = b; }
	void set_c(double const c) { this->_c = c; }
	double a() const { return _a; }
	double b() const { return _b; }
	double c() const { return _c; }
};

#endif


