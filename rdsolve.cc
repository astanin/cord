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
div_term(const AMesh2D& m, string var, string Dvar, int i, int j) {
	double fx=0.5/(m.get_dx()*m.get_dx());
	double fy=0.5/(m.get_dy()*m.get_dy());
	return fx*((m[Dvar](i+1,j)+m[Dvar](i,j))*(m[var](i+1,j)-m[var](i,j))
		  -(m[Dvar](i,j)+m[Dvar](i-1,j))*(m[var](i,j)-m[var](i-1,j)))
		+
	       fy*((m[Dvar](i,j+1)+m[Dvar](i,j))*(m[var](i,j+1)-m[var](i,j))
	          -(m[Dvar](i,j)+m[Dvar](i,j-1))*(m[var](i,j)-m[var](i,j-1)));
}

AMesh2D*
rd_step_euler(BCSet const& bcs, double dt, const AMesh2D& m1, string const var,
	string const Dvar, string const Rvar)
throw(MeshException) {
	try {
	auto_ptr<AMesh2D> m2(m1.clone());
	int xdim=m2->get_xdim();
	int ydim=m2->get_ydim();
	// for inner points
	for (int i=1; i<(xdim-1); ++i) {
		for (int j=1; j<(ydim-1); ++j) {
			double var2=m1.get(var,i,j)
				+dt*div_term(m1,var,Dvar,i,j)
				+dt*m1.get(Rvar,i,j);
			m2->set(var,i,j,var2);
		}
	}
	// boundary points
	for (int j=1; j < (ydim-1); ++j) {
		double dx=m2->get_dx();
		BoundaryCondition bc;
		// east boundary
		bc=bcs.get_east();
		m2->set(var,xdim-1,j,(bc.c()*dx+bc.b()*m2->get(var,xdim-2,j))/
					(bc.a()*dx+bc.b()));
		// west boundary
		bc=bcs.get_west();
		m2->set(var,0,j,(bc.c()*dx+bc.b()*m2->get(var,1,j))/
					(bc.a()*dx+bc.b()));
	}
	for (int i=0; i < xdim; ++i) {
		double dy=m2->get_dy();
		BoundaryCondition bc;
		// north boundary
		bc=bcs.get_north();
		m2->set(var,i,ydim-1,(bc.c()*dy+bc.b()*m2->get(var,i,ydim-2))/
					(bc.a()*dy+bc.b()));
		// south boundary
		bc=bcs.get_south();
		m2->set(var,i,0,(bc.c()*dy+bc.b()*m2->get(var,i,1))/
					(bc.a()*dy+bc.b()));
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

void
rd_step_implicit_fill_matrix(const AMesh2D& m, const BCSet& bcs,
	ASparseMatrix& A, vector<double>& rhs, MeshEnumerator const& k,
	double const dt, string const var, string const Dvar, string const Rvar)
throw(MeshException) {
	int xdim=m.get_xdim();
	int ydim=m.get_ydim();
	double dx=m.get_dx();
	double dy=m.get_dy();
	double fx=0.5*dt/(dx*dx);
	double fy=0.5*dt/(dy*dy);
	string D=Dvar;
	// reset right hand side
	for (int i=0; i<k.size(); ++i) {
		rhs.at(i)=0.0;
	}
	// for inner points
	for (int i=1; i<(xdim-1); ++i) {
		for (int j=1; j<(ydim-1); ++j) {
			int k0 =k(i  ,j);
			int k0p=k(i  ,j+1);
			int kp0=k(i+1,j);
			int k0m=k(i  ,j-1);
			int km0=k(i-1,j);
			// approximate diffusion term
			A.set(k0,k0,1.0+fx*(m[D](i+1,j))+
				2*(fx+fy)*(m[D](i,j))+fx*(m[D](i-1,j))
				+fy*(m[D](i,j+1))+fy*(m[D](i,j-1)));
			// the following is valid for stationary boundary
			// conditions only, using var from the previous time
			// step instead of bc.c()/bc.a()
			if (k0p >= 0) {
				A.set(k0,k0p,-fy*(m[D](i,j)+m[D](i,j+1)));
			} else {
				rhs.at(k0)+=fy*(m[D](i,j)+m[D](i,j+1))
						*m[var](i,j+1);
			}
			if (k0m >= 0) {
				A.set(k0,k0m,-fy*(m[D](i,j)+m[D](i,j-1)));
			} else {
				rhs.at(k0)+=fy*(m[D](i,j)+m[D](i,j-1))
						*m[var](i,j-1);
			}
			if (kp0 >= 0) {
				A.set(k0,kp0,-fx*(m[D](i,j)+m[D](i+1,j)));
			} else {
				rhs.at(k0)+=fx*(m[D](i,j)+m[D](i+1,j))
						*m[var](i+1,j);
			}
			if (km0 >= 0) {
				A.set(k0,km0,-fx*(m[D](i,j)+m[D](i-1,j)));
			} else {
				rhs.at(k0)+=fx*(m[D](i,j)+m[D](i-1,j))
						*m[var](i-1,j);
			}
			// approximate reaction term
			rhs.at(k0)+=m[var](i,j)+dt*m[Rvar](i,j);
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

AMesh2D*
rd_step_implicit(BCSet const& bcs, double dt, const AMesh2D& m1,
	string const var, string const Dvar, string const Rvar)
throw(MeshException) {
	try {
	SkipDirichletEnumerator kenum(m1,bcs);
	auto_ptr<AMesh2D> m2(m1.clone());
	int xdim=m2->get_xdim();
	int ydim=m2->get_ydim();
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

AMesh2D*
rd_step_adi(BCSet const& bcs, double dt, const AMesh2D& m1, string const var,
	string const Dvar, string const Rvar)
throw(MeshException) {
	try {
	auto_ptr<AMesh2D> m2(m1.clone());
	// ADI
#ifdef ENABLE_ADI_ISO_FIX
	if (rand()%2) { // choose sequence of directions: 50% x,y - 50% y,x
#endif
		m2.reset(step_peaceman_rachford_adi_x
			(bcs, 0.5*dt, *m2, var, Dvar, Rvar));
		m2.reset(step_peaceman_rachford_adi_y
			(bcs, 0.5*dt, *m2, var, Dvar, Rvar));
#ifdef ENABLE_ADI_ISO_FIX
	} else {
		m2.reset(step_peaceman_rachford_adi_y
			(bcs, 0.5*dt, *m2, var, Dvar, Rvar));
		m2.reset(step_peaceman_rachford_adi_x
			(bcs, 0.5*dt, *m2, var, Dvar, Rvar));
	}
#endif
	return m2.release();
	} catch(MeshException& e) {
		ostringstream ss;
		ss << "rd_step_adi: " << e.what();
		throw MeshException(ss.str());
	}
}


AMesh2D*
reaction_diffusion_step(BCSet const& bcs, double dt,
	AMesh2D const& m1, string const var,
	string const Dvar, string const Rvar,
	MP::rd_solver_t solver)
throw(MeshException) {
	auto_ptr<AMesh2D> m2(m1.clone());
	// default coefficients if not given
	string D;
	if (Dvar.empty()) { // assume D==1.0
		D="D_ONE";
		m2->remove_function_ifdef(D);
		m2->add_function(D,1.0);
	} else {
		D=Dvar;
	}
	string R;
	if (Rvar.empty()) { // assume reaction term R==0.0
		R="R_ZERO";
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
	if (Dvar.empty()) {
		m2->remove_function_ifdef(D);
	}
	if (Rvar.empty()) {
		m2->remove_function_ifdef(R);
	}
	return m2.release();
}

AMesh2D*
reaction_diffusion_step(BCSet const& bcs, double dt,
	AMesh2D const& m1, string const var,
	double const D, double const R,
	MP::rd_solver_t solver)
throw(MeshException) {
	auto_ptr<AMesh2D> m2(m1.clone());
	m2->remove_function_ifdef("Dvar_tmp");
	m2->remove_function_ifdef("Rvar_tmp");
	m2->add_function("Dvar_tmp",D);
	m2->add_function("Rvar_tmp",R);
	m2.reset(reaction_diffusion_step(bcs,dt,m1,var,
				"Dvar_tmp","Rvar_tmp",solver));
	m2->remove_function_ifdef("Dvar_tmp");
	m2->remove_function_ifdef("Rvar_tmp");
	return m2.release();
}

