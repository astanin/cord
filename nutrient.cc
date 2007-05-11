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

#include <iostream>
#include "nutrient.h"
#include "global.h"
#include "pradi.h"
#include "growth.h"
#include "meshenum.h"

#ifndef NUTRIENT_DUMB_EXPLICIT
#include "spmat.h"
#endif

#include <memory>

using std::cerr;
using std::string;
using std::vector;
using std::ostringstream;
using std::ios_base;
using std::setiosflags;
using std::setw;
using std::auto_ptr;

extern int verbose;

string dbg_stamp(double t);

/// norm-I: max(abs(f))
double norm_1(AMesh2D const& m, string const fid) {
	using namespace blitz;
	return max(fabs(m[fid]));
}

/// norm-II: sum(abs(f))
double norm_2(AMesh2D const& m, string const fid) {
	using namespace blitz;
	double norm=sum(fabs(m[fid]));
	return norm;
}

double
best_dt(AMesh2D const& m) throw(MeshException) {
	// WARNING: assuming uniform rectangular grid
	double dx=m.get_dx();
	double dy=m.get_dy();
	double h=(dx<dy?dx:dy);
	double dt=0.1*h*h/1.0; // D == 1.0 in non-dimensional model
	return dt;
}

double
consumption_rate(AMesh2D const& m, int const i, int const j)
throw(MeshException) {
	double consume=0.0;
	if (m.get("psi",i,j) > 0) { // tumour point
		consume=m.get("phi",i,j)*m.get_attr("consumption_c");
		double g=growth_term(m,i,j);
		consume+=m.get_attr("gconsumption_c")*0.5*(g+fabs(g));
	}
	return consume;
}

double
consumption_term(AMesh2D const& m, int const i, int const j)
throw(MeshException) {
	return m.get("c",i,j)*consumption_rate(m,i,j);
}

/** @brief fill @c A and @c rhs to construct equation for boundary point
 * @c k_boundary, which satisfies boundary condition @c bc and uses given inner
 * point @c k_inner and point distance @c dx
 */
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
eval_nutrient_build_sle_matrix(const AMesh2D& m, const BCSet& bcs,
	ASparseMatrix& A, vector<double>& rhs, MeshEnumerator const& k)
throw(MeshException) {
	int xdim=m.get_xdim();
	int ydim=m.get_ydim();
	int i, j;
	try {
		for (i=1; i<(xdim-1); ++i) {
			for (j=1; j<(ydim-1); ++j) {
				/* div(grad(c)) stencil:
				 *
				 *    (  )---[0p]---(  )
				 *      |     |      |
				 *      |     |      |
				 *    [m0]---[0 ]---[p0]
				 *      |     |      |
				 *      |     |      |
				 *    (  )---[0m]---(  )
				 */
				int k0 =k(i  ,j);
				int k0p=k(i  ,j+1);
				int kp0=k(i+1,j);
				int k0m=k(i  ,j-1);
				int km0=k(i-1,j);
				// WARNING: assuming rectangular uniform grid
				double dx2f=1.0/(m.get_dx()*m.get_dx());
				double dy2f=1.0/(m.get_dy()*m.get_dy());
				// fill in laplacian matrix
				A.set(k0,k0, -2.0*(dx2f+dy2f));
	// WARNING/TODO: this is valid only for stationary boundary conditions,
	// should use bc.c()/bc.a() for appropriate bc instead of m.get("c")
	// from the previous time step
				if (k0p >= 0) {
					A.set(k0,k0p, dy2f);
				} else {
					rhs.at(k0)=-dy2f*m.get("c",i,j+1);
				}
				if (k0m >= 0) {
					A.set(k0,k0m, dy2f);
				} else {
					rhs.at(k0)=-dy2f*m.get("c",i,j-1);
				}
				if (kp0 >= 0) {
					A.set(k0,kp0, dx2f);
				} else {
					rhs.at(k0)=-dx2f*m.get("c",i+1,j);
				}
				if (km0 >= 0) {
					A.set(k0,km0, dx2f);
				} else {
					rhs.at(k0)=-dx2f*m.get("c",i-1,j);
				}
				// add consumption term
				A.set(k0,k0,A.get(k0,k0)
					-consumption_rate(m,i,j));
			}
		}
	} catch (SparseMatrixException& e) {
		ostringstream ss;
		ss << "eval_nutrient_build_sle_matrix: "
			<< "failed to construct inner points equations: "
			<< " i=" << i << " j=" << j << ", reason: "
			<< e.what();
		throw MeshException(ss.str());
	}
	try {
		// construct equations for boundary points
		for (i=0; i<xdim; ++i) {
			int k0, k0m;
			double dy=m.get_dy();
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
		for (j=1; j<(ydim-1); ++j) {
			int k0, k0m;
			double dx=m.get_dx();
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
	} catch (SparseMatrixException& e) {
		ostringstream ss;
		ss << "eval_nutrient_build_sle_matrix: "
			<< "failed to construct boundary points equations: "
			<< e.what();
		throw MeshException(ss.str());
	}
}

void
update_dirichlet_points(AMesh2D& m, const BCSet& bcs, string const var) {
	int xdim=m.get_xdim();
	int ydim=m.get_ydim();
	for (int i=0; i<xdim; ++i) {
		BoundaryCondition bc;
		// north
		bc=bcs.get_north();
		if (bc.get_type() == BoundaryCondition::DIRICHLET_BC) {
			m.set(var,i,ydim-1,bc.c()/bc.a());
		}
		// south
		bc=bcs.get_south();
		if (bc.get_type() == BoundaryCondition::DIRICHLET_BC) {
			m.set(var,i,0,bc.c()/bc.a());
		}
	}
	for (int j=1; j<(ydim-1); ++j) {
		BoundaryCondition bc;
		// west
		bc=bcs.get_west();
		if (bc.get_type() == BoundaryCondition::DIRICHLET_BC) {
			m.set(var,0,j,bc.c()/bc.a());
		}
		// east
		bc=bcs.get_east();
		if (bc.get_type() == BoundaryCondition::DIRICHLET_BC) {
			m.set(var,xdim-1,j,bc.c()/bc.a());
		}
	}
}

AMesh2D*
eval_nutrient_directly(const Params& p, const AMesh2D& m1)
throw(MeshException) {
	try {
		SkipDirichletEnumerator kenum(m1,p.c_bc);
		auto_ptr<AMesh2D> m2(m1.clone());
		// initial estimate (previous solution)
		vector<double> c(kenum.size());
		for (int k=0; k<kenum.size(); ++k) {
			c.at(k)=m2->get("c",kenum.i(k),kenum.j(k));
		}
		auto_ptr<ASparseMatrix> pA;
#ifdef HAVE_LIBUMFPACK
		if (poisson_solver == SOLVER_UMFPACK) {
			pA.reset(new UMFPACKMatrix(kenum.size(),kenum.size()));
		} else
#endif
		if (poisson_solver == SOLVER_CG) {
			pA.reset(new LSolverMatrix(kenum.size(), c,
				poisson_solver_accuracy,
				poisson_solver_max_iterations,
				LSolverMatrix::CG));
		} else if (poisson_solver == SOLVER_BICG) {
			pA.reset(new LSolverMatrix(kenum.size(), c,
				poisson_solver_accuracy,
				poisson_solver_max_iterations,
				LSolverMatrix::BICG));
		} else if (poisson_solver == SOLVER_BICGSTAB) {
			pA.reset(new LSolverMatrix(kenum.size(), c,
				poisson_solver_accuracy,
				poisson_solver_max_iterations,
				LSolverMatrix::BICGstab));
		} else if (poisson_solver == SOLVER_GMRES) {
			pA.reset(new LSolverMatrix(kenum.size(), c,
				poisson_solver_accuracy,
				poisson_solver_max_iterations,
				LSolverMatrix::GMRES, gmres_restart_after));
		} else {
			throw MeshException("unknown poisson_solver");
		}
		vector<double> rhs(kenum.size());// right hand side vector
		eval_nutrient_build_sle_matrix(*m2,p.c_bc,*pA,rhs,kenum);
		// solve SLE
		try {
			c=pA->solve(rhs);
		} catch (SparseMatrixException& e) {
			ostringstream ss;
			ss << "SparseMatrixException exception in solve(): "
				<< e.what() << "\n";
			throw MeshException(ss.str());
		}
		// update m2
		for (int k=0; k<(int)c.size(); ++k) {
			int i=kenum.i(k);
			int j=kenum.j(k);
			m2->set("c",i,j,c.at(k));
		}
		// update Dirichlet boundary points (excluded from SLE)
		update_dirichlet_points(*m2,p.c_bc,"c");
		return m2.release();
	} catch (MeshException& e) {
		ostringstream ss;
		ss << "eval_nutrient_directly: " << e.what();
		throw MeshException(ss.str());
	}
}

double
eval_nutrient_residual(const AMesh2D& m)
throw(MeshException) {
	double global_res=0.0;
	double dx=m.get_dx();
	double dy=m.get_dy();
	const auto_ptr<AMesh2D> m2(m.clone());
	m2->remove_function_ifdef("c_residual");
	m2->add_function("c_residual");
	// WARNING: ignoring boundary points
	// inner points
	for (int i=1; (i<m2->get_xdim()-1) ; ++i) {
		for (int j=1; (j<m2->get_ydim()-1) ; ++j) {
			double loc_res=0.0;
			loc_res=(m2->get("c",i+1,j)-2*m2->get("c",i,j)
					+m2->get("c",i-1,j))/(dx*dx)+
				(m2->get("c",i,j+1)-2*m2->get("c",i,j)
				 	+m2->get("c",i,j-1))/(dy*dy)-
				consumption_term(*m2,i,j);
			loc_res=fabs(loc_res);
			m2->set("c_residual",i,j,loc_res);
		}
	}
	global_res=norm_1(*m2,"c_residual");
	return global_res;
}

void
eval_nutrient_set_boundaries(const BCSet& bcs, AMesh2D& m)
throw(MeshException) {
	int xdim=m.get_xdim();
	int ydim=m.get_ydim();
	double dx=m.get_dx();
	double dy=m.get_dy();
	BoundaryCondition bc;
	for (int i=0; i<xdim; ++i) {
		// north boundary
		bc=bcs.get_north();
		m.set("c",i,ydim-1,(bc.c()*dy+bc.b()*m.get("c",i,ydim-2))/
					(bc.a()*dy+bc.b()));
		// south boundary
		bc=bcs.get_south();
		m.set("c",i,0,(bc.c()*dy+bc.b()*m.get("c",i,1))/
					(bc.a()*dy+bc.b()));
	}
	for (int j=0; j<ydim; ++j) {
		// west boundary
		bc=bcs.get_west();
		m.set("c",0,j,(bc.c()*dx+bc.b()*m.get("c",1,j))/
					(bc.a()*dx+bc.b()));
		// east boundary
		bc=bcs.get_east();
		m.set("c",xdim-1,j,(bc.c()*dx+bc.b()*m.get("c",xdim-2,j))/
					(bc.a()*dx+bc.b()));
	}
}

AMesh2D*
eval_nutrient_iteratively_step_explicit(const Params& p, double const dt,
	const AMesh2D& m)
throw(MeshException) {
	double dx=m.get_dx();
	double dy=m.get_dy();
	auto_ptr<AMesh2D> m2(m.clone());
	// inner points
	for (int i=1; (i<m.get_xdim()-1) ; ++i) {
		for (int j=1; (j<m.get_ydim()-1) ; ++j) {
			double loc_res=
				(m.get("c",i+1,j)-2*m.get("c",i,j)
					+m.get("c",i-1,j))/(dx*dx)+
				(m.get("c",i,j+1)-2*m.get("c",i,j)
				 	+m.get("c",i,j-1))/(dy*dy)-
				consumption_term(m,i,j);
			m2->set("c",i,j,m.get("c",i,j)+dt*loc_res);
		}
	}
	// boundary points
	eval_nutrient_set_boundaries(p.c_bc,*m2);
	return m2.release();
}

AMesh2D*
eval_nutrient_iteratively_step_implicit(const Params& p, double const dt,
	const AMesh2D& m)
throw(MeshException) {
	auto_ptr<AMesh2D> m2(m.clone());
	m2->add_function_ifndef("consumption");
	// initialize consumption variable
	int xdim=m.get_xdim();
	int ydim=m.get_ydim();
	for (int i=0; i<xdim ; ++i) {
		for (int j=0; j<ydim ; ++j) {
			m2->set("consumption",i,j,consumption_term(m,i,j));
		}
	}
	// ADI
	m2.reset(step_peaceman_rachford_adi_x
		(p.c_bc, 0.5*dt, *m2, "c", "", "consumption"));
	m2.reset(step_peaceman_rachford_adi_y
		(p.c_bc, 0.5*dt, *m2, "c", "", "consumption"));
	// finishing
	m2->remove_function_ifdef("consumption");
	return m2.release();
}

AMesh2D*
eval_nutrient_iteratively(const Params& p, const AMesh2D& m,
		double const epsilon)
throw(MeshException) {
	try {
	auto_ptr<AMesh2D> m2(m.clone());
	double dt;
	if (poisson_solver_iteration_step < 0) { // automatic time step
		if (poisson_solver == SOLVER_ITERATIVE_EXPLICIT) {
			dt=best_dt(*m2);
		} else { // implicit
			double dx=m.get_dx();
			double dy=m.get_dy();
			dt=10.0*(dx*dx+dy*dy); // heuristic
		}
	} else {
		dt=poisson_solver_iteration_step;
	}
	if (verbose > 1) { // extra verbose
		cerr << dbg_stamp(m.get_time())
			<< "eval_nutrient_iteratively: dt=" << dt;
		if (poisson_solver == SOLVER_ITERATIVE_EXPLICIT) {
			cerr << " explicit method" << "\n";
		} else {
			cerr << " implicit method" << "\n";
		}
	}
	// iterative method here
	double res=0.0;
	res=eval_nutrient_residual(*m2);
	int iter=0;
	while ((res > epsilon) && (iter <= poisson_solver_max_iterations)) {
		if (poisson_solver == SOLVER_ITERATIVE_EXPLICIT) {
			m2.reset(
			eval_nutrient_iteratively_step_explicit(p,dt,*m2));
		} else {
			m2.reset(
			eval_nutrient_iteratively_step_implicit(p,dt,*m2));
		}
		res=eval_nutrient_residual(*m2);
		++iter;
		if ((verbose > 1) && (!(iter%50))) {
			cerr << dbg_stamp(m2->get_time())
				<< "eval_nutrient_iteratively: iters=" << iter
					<< " residual=" << res << "\n";
		}
	};
	if (verbose > 1) {
		cerr << dbg_stamp(m2->get_time())
			<< "eval_nutrient_iteratively: iters=" << iter
				<< " residual=" << res << "\n";
	}
	return m2.release();
	} catch (MeshException& e) {
		ostringstream ss;
		ss << "eval_nutrient_iteratively: " << e.what();
		throw MeshException(ss.str());
	}
}

AMesh2D*
eval_nutrient(const Params& p, const AMesh2D& m1, double const epsilon)
throw(MeshException) {
	if (!m1.defined("c")) {
		throw MeshException("eval_nutrient: c not defined");
	}
	if (!m1.defined("phi")) {
		throw MeshException("eval_nutrient: phi not defined");
	}
	if (!m1.defined("psi")) {
		throw MeshException("eval_nutrient: psi not defined");
	}
	if (!m1.attr_defined("consumption_c")) {
		throw MeshException("eval_nutrient: "
			"basic nutrient consumption rate not defined");
	}
	if (!m1.attr_defined("gconsumption_c")) {
		throw MeshException("eval_nutrient: "
			"growth-related nutrient consumption rate not defined");
	}
	try {
		if (p.c_equation != Params::EQ_POISSON) {
			throw MeshException(
				"unsupported nutrient equation type");
		}
		AMesh2D *m2;
		if (
#ifdef HAVE_LIBUMFPACK
			poisson_solver == SOLVER_UMFPACK ||
#endif
			poisson_solver == SOLVER_CG ||
			poisson_solver == SOLVER_GMRES ||
			poisson_solver == SOLVER_BICG ||
			poisson_solver == SOLVER_BICGSTAB) {
			m2=eval_nutrient_directly(p, m1);
		} else if (poisson_solver == SOLVER_ITERATIVE_IMPLICIT ||
			poisson_solver == SOLVER_ITERATIVE_EXPLICIT) {
			m2=eval_nutrient_iteratively(p, m1, epsilon);
		} else {
			throw MeshException("eval_nutrient: unknown method");
		}
		return m2;
	} catch (MeshException& e) {
		ostringstream ss;
		ss << "eval_nutrient: " << e.what();
		throw MeshException(ss.str());
	}
}

