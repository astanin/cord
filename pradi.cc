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

#include <stdlib.h>

#include "pradi.h"
#include <memory>

using std::auto_ptr;

string dbg_stamp(double t);

template<class fid_t>
AMesh2D<fid_t>*
step_peaceman_rachford_adi_x(const BCSet& bcs,
	double dt, const AMesh2D<fid_t>& m1, fid_t const var,
	fid_t const D_coef_var, fid_t const reaction_term_var)
throw(MeshException) {
	try {
	auto_ptr<AMesh2D<fid_t> > m2(m1.clone());
	int xdim=m1.get_xdim();
	int ydim=m1.get_ydim();
	// WARNING: valid for uniform grid only
	double dx=m1.get_dx();
	double dy=m1.get_dy();
	array2d u=m1[var];
	array2d D=m1[D_coef_var];
	array2d R=m1[reaction_term_var];
	array2d m2_u=(*m2)[var];
	for (int j=1; j<(ydim-1); ++j) {
		// solution vector
		vector<double> u2(xdim);
#ifdef HAVE_LIBLAPACK
		TridiagMatrix A(xdim);
#else
#ifdef HAVE_LIBUMFPACK
		UMFPACKMatrix A(xdim,xdim);
#else
		LSolverMatrix A(xdim, u2, 1e-5 /* accuracy */,
				1000 /* max iters */, LSolverMatrix::BICG);
#endif
#endif
		vector<double> rhs(xdim);
		// fill tridiagonal matrix and right-hand-side
		// inner points
		for (int i=1; i<(xdim-1); ++i) {
			double D_0, D_ip, D_im, D_jp, D_jm;
			D_0=D(i,j);
			D_ip=D(i+1,j);
			D_im=D(i-1,j);
			D_jp=D(i,j+1);
			D_jm=D(i,j-1);
			A.set(i,i,  2*dx*dx+dt*(D_ip+D_im+2*D_0));
			A.set(i,i+1,-dt*(D_ip+D_0));
			A.set(i,i-1,-dt*(D_im+D_0));
			double u_0=u(i,j);
			double u_jp=u(i,j+1);
			double u_jm=u(i,j-1);
			double f=-R(i,j);
			rhs.at(i)= 2*dx*dx*u_0 +
				(dt*dx*dx/(dy*dy))*(
					+ (D_jp+D_0)*u_jp
					- (D_jp+2*D_0+D_jm)*u_0
					+ (D_0+D_jm)*u_jm
				) + (2*dt*dx*dx)*f;
		}
		// boundary conditions
		BoundaryCondition bc;
		// west boundary
		bc=bcs.get_west();
		A.set(0,0,bc.a()*dx+bc.b());
		A.set(0,1,-bc.b());
		rhs.at(0)=bc.c()*dx;
		// east boundary
		bc=bcs.get_east();
		A.set(xdim-1,xdim-1,bc.a()*dx+bc.b());
		A.set(xdim-1,xdim-2,-bc.b());
		rhs.at(xdim-1)=bc.c()*dx;
		// solution
		u2=A.solve(rhs);
		for (int i=0; i<xdim; ++i) {
			m2_u(i,j)=u2[i];
		}
	}
	// satisfy other boundary conditions (along y-axis)
	for (int i=0; i<xdim; ++i) {
		BoundaryCondition bc;
		// north boundary
		bc=bcs.get_north();
		m2_u(i,ydim-1)=(bc.c()*dy+bc.b()*m2_u(i,ydim-2))/
					(bc.a()*dy+bc.b());
		// south boundary
		bc=bcs.get_south();
		m2_u(i,0)=(bc.c()*dy+bc.b()*m2_u(i,1))/(bc.a()*dy+bc.b());
	}
	return m2.release();
	} catch(MeshException& e) {
		ostringstream ss;
		ss << "step_peaceman_rachford_adi_x: " << e.what();
		throw MeshException(ss.str());
	}
}

template<class fid_t>
AMesh2D<fid_t>*
step_peaceman_rachford_adi_y(const BCSet& bcs,
	double dt, const AMesh2D<fid_t>& m1, fid_t const var,
	fid_t const D_coef_var, fid_t const reaction_term_var)
throw(MeshException) {
	try {
	auto_ptr<AMesh2D<fid_t> > m2(m1.clone());
	int xdim=m1.get_xdim();
	int ydim=m1.get_ydim();
	// WARNING: valid for uniform grid only
	double dx=m1.get_dx();
	double dy=m1.get_dy();
	array2d u=m1[var];
	array2d D=m1[D_coef_var];
	array2d R=m1[reaction_term_var];
	array2d m2_u=(*m2)[var];
	for (int i=1; i<(xdim-1); ++i) {
		// solution vector
		vector<double> u2(ydim);
#ifdef HAVE_LIBLAPACK
		TridiagMatrix A(ydim);
#else
#ifdef HAVE_LIBUMFPACK
		UMFPACKMatrix A(ydim,ydim);
#else
		LSolverMatrix A(ydim, u2, 1e-5 /* accuracy */,
				1000 /* max iters */, LSolverMatrix::BICG);
#endif
#endif
		vector<double> rhs(ydim);
		// fill tridiagonal matrix and right-hand-side
		// inner points
		for (int j=1; j<(ydim-1); ++j) {
			double D_0, D_ip, D_im, D_jp, D_jm;
			D_0=D(i,j);
			D_ip=D(i+1,j);
			D_im=D(i-1,j);
			D_jp=D(i,j+1);
			D_jm=D(i,j-1);
			A.set(j,j, 2*dy*dy+dt*(D_jp+D_jm+2*D_0));
			A.set(j,j+1,-dt*(D_jp+D_0));
			A.set(j,j-1,-dt*(D_jm+D_0));
			double u_0=u(i,j);
			double u_ip=u(i+1,j);
			double u_im=u(i-1,j);
			double f=-R(i,j);
			rhs.at(j)=2*dy*dy*u_0 +
				(dt*dy*dy/(dx*dx))*(
					+ (D_ip+D_0)*u_ip
					- (D_ip+2*D_0+D_im)*u_0
					+ (D_0+D_im)*u_im
				) + (2*dt*dy*dy)*f;
		}
		// boundary conditions
		BoundaryCondition bc;
		// south boundary
		bc=bcs.get_south();
		A.set(0,0,bc.a()*dy+bc.b());
		A.set(0,1,-bc.b());
		rhs.at(0)=bc.c()*dy;
		// north boundary
		bc=bcs.get_north();
		A.set(ydim-1,ydim-1,bc.a()*dy+bc.b());
		A.set(ydim-1,ydim-2,-bc.b());
		rhs.at(ydim-1)=bc.c()*dy;
		// solution
		u2=A.solve(rhs);
		for (int j=0; j<ydim; ++j) {
			m2_u(i,j)=u2[j];
		}
	}
	// satisfy other boundary conditions (along x-axis)
	for (int j=0; j<ydim; ++j) {
		BoundaryCondition bc;
		// west
		bc=bcs.get_west();
		m2_u(0,j)=(bc.c()*dx+bc.b()*m2_u(1,j))/(bc.a()*dx+bc.b());
		// east
		bc=bcs.get_east();
		m2_u(xdim-1,j)=(bc.c()*dx+bc.b()*m2_u(xdim-2,j))/
					(bc.a()*dx+bc.b());
	}
	return m2.release();
	} catch(MeshException& e) {
		ostringstream ss;
		ss << "step_peaceman_rachford_adi_y: " << e.what();
		throw MeshException(ss.str());
	}
}

// templates
template
AMesh2D<int>*
step_peaceman_rachford_adi_x<int>(const BCSet& bcs,
	double dt, const AMesh2D<int>& m1, int const var,
	int const D_coef_var, int const reaction_term_var);
template
AMesh2D<int>*
step_peaceman_rachford_adi_y<int>(const BCSet& bcs,
	double dt, const AMesh2D<int>& m1, int const var,
	int const D_coef_var, int const reaction_term_var);

