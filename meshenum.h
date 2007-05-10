/* Tumour cord growth model */

/*
 * Copyright (C) 2005-2006 Sergey Astanin, Luigi Preziosi
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

#include "config.h"

#ifndef MESHENUM_H
#define MESHENUM_H

class EnumeratorException {
public:
	EnumeratorException(int i=-1, int j=-1, int k=-1) : i(i), j(j), k(k) {}
	int i;
	int j;
	int k;
};

/**
 * MeshEnumerator class provides an interface to enumerate all or some of
 * the mesh points
 */
class MeshEnumerator {
public:
	virtual ~MeshEnumerator() {}
	virtual int size(void) const = 0;
	virtual int operator()(int const i, int const j) const = 0;
	virtual int k(int const i, int const j) const {
		return this->operator()(i,j);
	}
	virtual int i(int const k) const = 0;
	virtual int j(int const k) const = 0;
};

#include "amesh2d.h"

class FullEnumerator : public MeshEnumerator {
	int xdim;
	int ydim;
public:
	FullEnumerator(AMesh2D const& m) :
		xdim(m.get_xdim()), ydim(m.get_ydim()) {}
	virtual ~FullEnumerator() {}
	virtual int size(void) const {
		return xdim*ydim;
	}
	virtual int operator()(int const i, int const j) const {
		return (i+(xdim*j));
	}
	virtual int i(int const k) const {
		return k%xdim;
	}
	virtual int j(int const k) const {
		return k/xdim;
	}
};

#include "params.h"
#include <iostream>

/**
 * SkipDirichletEnumerator enumerates only those mesh points where
 * Dirichlet boundary condition is not defined.
 */
class SkipDirichletEnumerator : public MeshEnumerator {
private:
	int xdim;
	int ydim;
	int xmin;
	int xmax;
	int ymin;
	int ymax;
public:
	SkipDirichletEnumerator(AMesh2D const& m, BCSet const& bcs)
		: xdim(m.get_xdim()), ydim(m.get_ydim()) {
		xmin=0;
		xmax=xdim-1;
		ymin=0;
		ymax=ydim-1;
		if (bcs.get_north().get_type() == BC::DIRICHLET_BC) {
			ymax--;
			ydim--;
		}
		if (bcs.get_south().get_type() == BC::DIRICHLET_BC) {
			ymin++;
			ydim--;
		}
		if (bcs.get_west().get_type() == BC::DIRICHLET_BC) {
			xmin++;
			xdim--;
		}
		if (bcs.get_east().get_type() == BC::DIRICHLET_BC) {
			xmax--;
			xdim--;
		}
	}
	virtual ~SkipDirichletEnumerator() {}
	virtual int size(void) const {
		return xdim*ydim;
	}
	virtual int operator()(int const i, int const j) const {
		if ((i>xmax) || (i<xmin) || (j>ymax) || (j<ymin)) {
			return -1;
		}
		return ((i-xmin)+(xdim*(j-ymin)));
	}
	virtual int i(int const k) const {
		if ((k<0) || (k>size())) {
			throw EnumeratorException(-1,-1,k);
		}
		return xmin+k%xdim;
	}
	virtual int j(int const k) const {
		if ((k<0) || (k>size())) {
			throw EnumeratorException(-1,-1,k);
		}
		return ymin+k/xdim;
	}
};

#endif


