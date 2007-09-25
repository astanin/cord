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
 *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>

#include "slesolve.h"
#include "rdsolve.h"
#include "pradi.h"
#include "meshenum.h"

#include <memory>
using std::auto_ptr;

#ifdef ENABLE_ADI_ISO_FIX
#include <stdlib.h>
#endif

#include <iostream>
using std::cerr;

double
div_term(array2d const& var, array2d const& Dvar, int i, int j,
	double dx, double dy) {
	double fx=0.5/(dx*dx);
	double fy=0.5/(dy*dy);
	return fx*((Dvar(i+1,j)+Dvar(i,j))*(var(i+1,j)-var(i,j))
		  -(Dvar(i,j)+Dvar(i-1,j))*(var(i,j)-var(i-1,j)))
		+
	       fy*((Dvar(i,j+1)+Dvar(i,j))*(var(i,j+1)-var(i,j))
	          -(Dvar(i,j)+Dvar(i,j-1))*(var(i,j)-var(i,j-1)));
}

template<class fid_t>
AMesh2D<fid_t>*
rd_step_euler(BCSet const& bcs, double dt, const AMesh2D<fid_t>& m1,
	fid_t const var, fid_t const Dvar, fid_t const Rvar)
throw(MeshException) {
	try {
	auto_ptr<AMesh2D<fid_t> > m2(m1.clone());
	int xdim=m2->get_xdim();
	int ydim=m2->get_ydim();
	// for inner points
	array2d m1_var=m1[var];
	array2d m1_R=m1[Rvar];
	array2d m1_D=m1[Dvar];
	array2d m2_var=(*m2)[var];
	double dx=m1.get_dx();
	double dy=m1.get_dy();
	for (int i=1; i<(xdim-1); ++i) {
		for (int j=1; j<(ydim-1); ++j) {
			double var2=m1_var(i,j)
				+dt*div_term(m1_var,m1_D,i,j,dx,dy)
				+dt*m1_R(i,j);
			m2_var(i,j)=var2;
		}
	}
	// boundary points
	for (int j=1; j < (ydim-1); ++j) {
		BoundaryCondition bc;
		// east boundary
		bc=bcs.get_east();
		m2_var(xdim-1,j)=(bc.c()*dx+bc.b()*m2_var(xdim-2,j))/
					(bc.a()*dx+bc.b());
		// west boundary
		bc=bcs.get_west();
		m2_var(0,j)=(bc.c()*dx+bc.b()*m2_var(1,j))/(bc.a()*dx+bc.b());
	}
	for (int i=0; i < xdim; ++i) {
		BoundaryCondition bc;
		// north boundary
		bc=bcs.get_north();
		m2_var(i,ydim-1)=(bc.c()*dy+bc.b()*m2_var(i,ydim-2))/
					(bc.a()*dy+bc.b());
		// south boundary
		bc=bcs.get_south();
		m2_var(i,0)=(bc.c()*dy+bc.b()*m2_var(i,1))/(bc.a()*dy+bc.b());
	}
	return m2.release();
	} catch(MeshException& e) {
		ostringstream ss;
		ss << "rd_step_euler: " << e.what();
		throw MeshException(ss.str());
	}
}

void
build_boundary_point_eq(ASparseMatrix& A, vector<double>& rhs,
	int const k_boundary, int const k_inner,
	BoundaryCondition const& bc, double const dx) {
	int kb=k_boundary;
	int ki=k_inner;
	double a_ki_kb;
	switch (bc.get_type()) {
	case BoundaryCondition::NEUMANN_BC:
		a_ki_kb=A.get(ki,kb);
		A.set(kb,kb ,-a_ki_kb);
		A.set(kb,ki, a_ki_kb);
		rhs.at(kb)=bc.c()*dx*a_ki_kb/bc.b();
		break;
	case BoundaryCondition::DIRICHLET_BC: // do nothing
		break;
	default:
		ostringstream ss;
		ss << "build_boundary_point_eq: "
			<< "unsupported boundary condition ("
			<< bc.get_type()<<") in "<< ki<<"-th equation";
		throw MeshException(ss.str());
		break;
	}
}

template<class fid_t>
void
rd_step_implicit_fill_matrix(const AMesh2D<fid_t>& m, const BCSet& bcs,
	ASparseMatrix& A, vector<double>& rhs, MeshEnumerator const& k,
	double const dt, fid_t const var, fid_t const Dvar, fid_t const Rvar)
throw(MeshException) {
	int xdim=m.get_xdim();
	int ydim=m.get_ydim();
	double dx=m.get_dx();
	double dy=m.get_dy();
	double fx=0.5*dt/(dx*dx);
	double fy=0.5*dt/(dy*dy);
	fid_t D=Dvar;
	// reset right hand side
	int ksize=k.size();
	for (int i=0; i<ksize; ++i) {
		rhs.at(i)=0.0;
	}
	// for inner points
	array2d m_D=m[D];
	array2d m_var=m[var];
	array2d m_R=m[Rvar];
	for (int i=1; i<(xdim-1); ++i) {
		for (int j=1; j<(ydim-1); ++j) {
			int k0 =k(i  ,j);
			int k0p=k(i  ,j+1);
			int kp0=k(i+1,j);
			int k0m=k(i  ,j-1);
			int km0=k(i-1,j);
			// approximate diffusion term
			A.set(k0,k0,1.0+fx*(m_D(i+1,j))+
				2*(fx+fy)*(m_D(i,j))+fx*(m_D(i-1,j))
				+fy*(m_D(i,j+1))+fy*(m_D(i,j-1)));
			// the following is valid for stationary boundary
			// conditions only, using var from the previous time
			// step instead of bc.c()/bc.a()
			if (k0p >= 0) {
				A.set(k0,k0p,-fy*(m_D(i,j)+m_D(i,j+1)));
			} else {
				rhs.at(k0)+=fy*(m_D(i,j)+m_D(i,j+1))
						*m_var(i,j+1);
			}
			if (k0m >= 0) {
				A.set(k0,k0m,-fy*(m_D(i,j)+m_D(i,j-1)));
			} else {
				rhs.at(k0)+=fy*(m_D(i,j)+m_D(i,j-1))
						*m_var(i,j-1);
			}
			if (kp0 >= 0) {
				A.set(k0,kp0,-fx*(m_D(i,j)+m_D(i+1,j)));
			} else {
				rhs.at(k0)+=fx*(m_D(i,j)+m_D(i+1,j))
						*m_var(i+1,j);
			}
			if (km0 >= 0) {
				A.set(k0,km0,-fx*(m_D(i,j)+m_D(i-1,j)));
			} else {
				rhs.at(k0)+=fx*(m_D(i,j)+m_D(i-1,j))
						*m_var(i-1,j);
			}
			// approximate reaction term
			rhs.at(k0)+=m_var(i,j)+dt*m_R(i,j);
		}
	}
	// construct equations for boundary points
	for (int i=1; i<(xdim-1); ++i) {
		int k0, k0m;
		BoundaryCondition bc;
		// north
		bc=bcs.get_north();
		k0 =k(i,ydim-1);
		k0m=k(i,ydim-2);
		build_boundary_point_eq(A,rhs,k0,k0m,bc,dy);
		// south
		bc=bcs.get_south();
		k0 =k(i,0);
		k0m=k(i,1);
		build_boundary_point_eq(A,rhs,k0,k0m,bc,dy);
	}
	for (int j=1; j<(ydim-1); ++j) {
		int k0, k0m;
		BoundaryCondition bc;
		// west
		bc=bcs.get_west();
		k0 =k(0,j);
		k0m=k(1,j);
		build_boundary_point_eq(A,rhs,k0,k0m,bc,dx);
		// east
		bc=bcs.get_east();
		k0 =k(xdim-1,j);
		k0m=k(xdim-2,j);
		build_boundary_point_eq(A,rhs,k0,k0m,bc,dx);
	}
	return;
}

template<class fid_t>
AMesh2D<fid_t>*
rd_step_implicit(BCSet const& bcs, double dt, const AMesh2D<fid_t>& m1,
	fid_t const var, fid_t const Dvar, fid_t const Rvar)
throw(MeshException) {
	try {
	SkipDirichletEnumerator<fid_t> kenum(m1,bcs);
	auto_ptr<AMesh2D<fid_t> > m2(m1.clone());
	// initial estimate (previous solution)
	vector<double> u0(kenum.size());
	for (int k=0; k<kenum.size(); ++k) {
		u0.at(k)=m2->get(var,kenum.i(k),kenum.j(k));
	}
	// SLE matrix
	auto_ptr<ASparseMatrix> pA(build_sle_solver_matrix(kenum.size(), u0));
	vector<double> rhs(kenum.size(),0.0);
	rd_step_implicit_fill_matrix(*m2,bcs,*pA,rhs,kenum,dt,var,Dvar,Rvar);
	// solve SLE
	try {
		u0=pA->solve(rhs);
	} catch (SparseMatrixException& e) {
		ostringstream ss;
		ss << "SparseMatrixException in solve(): " << e.what();
		throw MeshException(ss.str());
	}
	// update m2
	for (int k=0; k<(int)u0.size(); ++k) {
		int i=kenum.i(k);
		int j=kenum.j(k);
		m2->set(var,i,j,u0.at(k));
	}
	// update Dirichlet boundary points (excluded from SLE)
	update_dirichlet_points(*m2,bcs,var);
	return m2.release();
	} catch(MeshException& e) {
		ostringstream ss;
		ss << "rd_step_implicit: " << e.what();
		throw MeshException(ss.str());
	}
}

template<class fid_t>
AMesh2D<fid_t>*
rd_step_adi
(BCSet const& bcs, double dt, const AMesh2D<fid_t>& m1, fid_t const var,
	fid_t const Dvar, fid_t const Rvar)
throw(MeshException) {
	try {
	auto_ptr<AMesh2D<fid_t> > m2(m1.clone());
	fid_t D=Dvar;
	fid_t R=Rvar;
	if (Dvar == NONE) {
		m2->add_function_ifndef(D_ONE);
		array2d d=(*m2)[D_ONE];
		d=1.0;
		D=D_ONE;
	}
	if (Rvar == NONE) {
		m2->add_function_ifndef(R_ZERO);
		array2d r=(*m2)[R_ZERO];
		r=0.0;
		R=R_ZERO;
	}
	// ADI
#ifdef ENABLE_ADI_ISO_FIX
	if (rand()%2) { // choose sequence of directions: 50% x,y - 50% y,x
#endif
		m2.reset(step_peaceman_rachford_adi_x<fid_t>
			(bcs, 0.5*dt, *m2, var, D, R));
		m2.reset(step_peaceman_rachford_adi_y<fid_t>
			(bcs, 0.5*dt, *m2, var, D, R));
#ifdef ENABLE_ADI_ISO_FIX
	} else {
		m2.reset(step_peaceman_rachford_adi_y<fid_t>
			(bcs, 0.5*dt, *m2, var, D, R));
		m2.reset(step_peaceman_rachford_adi_x<fid_t>
			(bcs, 0.5*dt, *m2, var, D, R));
	}
#endif
	if (Dvar == NONE) {
		m2->remove_function_ifdef(D_ONE);
	}
	if (Rvar == NONE) {
		m2->remove_function_ifdef(R_ZERO);
	}
	return m2.release();
	} catch(MeshException& e) {
		ostringstream ss;
		ss << "rd_step_adi: " << e.what();
		throw MeshException(ss.str());
	}
}


template<class fid_t>
AMesh2D<fid_t>*
reaction_diffusion_step(BCSet const& bcs, double dt,
	AMesh2D<fid_t> const& m1, fid_t const var,
	fid_t const Dvar, fid_t const Rvar,
	MP::rd_solver_t solver)
throw(MeshException) {
	auto_ptr<AMesh2D<fid_t> > m2(m1.clone());
	// default coefficients if not given
	fid_t D;
	if (Dvar == NONE) { // assume D==1.0
		D=D_ONE;
		m2->remove_function_ifdef(D);
		m2->add_function(D,1.0);
	} else {
		D=Dvar;
	}
	fid_t R;
	if (Rvar == NONE) { // assume reaction term R==0.0
		R=R_ZERO;
		m2->remove_function_ifdef(R);
		m2->add_function(R,0.0);
	} else {
		R=Rvar;
	}
	// solve
	switch(solver) {
	case MP::RDS_EXPLICIT:
		m2.reset(rd_step_euler(bcs,dt,*m2,var,D,R));
		break;
	case MP::RDS_ADI:
		m2.reset(rd_step_adi(bcs,dt,*m2,var,D,R));
		break;
	case MP::RDS_IMPLICIT:
		m2.reset(rd_step_implicit(bcs,dt,*m2,var,D,R));
		break;
	default:
		throw MeshException("reaction_diffusion_step: unknown solver");
		break;
	}
	// remove added functions
	if (Dvar == NONE) {
		m2->remove_function_ifdef(D);
	}
	if (Rvar == NONE) {
		m2->remove_function_ifdef(R);
	}
	return m2.release();
}

template<class fid_t>
AMesh2D<fid_t>*
reaction_diffusion_step(BCSet const& bcs, double dt,
	AMesh2D<fid_t> const& m1, fid_t const var,
	double const D, double const R,
	MP::rd_solver_t solver)
throw(MeshException) {
	auto_ptr<AMesh2D<fid_t> > m2(m1.clone());
	m2->remove_function_ifdef(D_TMP);
	m2->remove_function_ifdef(Q_TMP);
	m2->add_function(D_TMP,D);
	m2->add_function(Q_TMP,R);
	m2.reset(reaction_diffusion_step(bcs,dt,m1,var,
				D_TMP,Q_TMP,solver));
	m2->remove_function_ifdef(D_TMP);
	m2->remove_function_ifdef(Q_TMP);
	return m2.release();
}

// templates

template
AMesh2D<int>*
reaction_diffusion_step<int>(BCSet const& bcs, double dt,
	AMesh2D<int> const& m1, int const var,
	int const Dvar, int const Rvar,
	MP::rd_solver_t solver);

