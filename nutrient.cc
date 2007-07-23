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
#include "utils.h"
#include "rdsolve.h"
#include "growth.h"
#include "meshenum.h"
#include "slesolve.h"

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

template<class fid_t>
double
best_dt(AMesh2D<fid_t> const& m) throw(MeshException) {
	// WARNING: assuming uniform rectangular grid
	double dx=m.get_dx();
	double dy=m.get_dy();
	double h=(dx<dy?dx:dy);
	double dt=0.1*h*h/1.0; // D == 1.0 in non-dimensional model
	return dt;
}

template<class fid_t>
double
uptake_per_cell(AMesh2D<fid_t> const& m, array2d const& phi, array2d const& psi,
			int const i, int const j) {
	double consume=0.0;
	if (psi(i,j) > 0) { // tumour point
		double p=phi(i,j);
		consume=p*f_atp_per_cell(p)*m.get_attr("o2_uptake");
	}
	return consume;
}

template<class fid_t>
double
uptake_per_cell(AMesh2D<fid_t> const& m, int const i, int const j)
throw(MeshException) {
	return uptake_per_cell(m,m[PHI],m[PSI],i,j);
}

template<class fid_t>
double
uptake_term(AMesh2D<fid_t> const& m, int const i, int const j)
throw(MeshException) {
	double psi=m[PSI](i,j);
	double host_activity=m.get_attr("host_activity");
	double theta=m.get_attr("upkeep_per_cell");
	double alpha=m.get_attr("o2_uptake");
	double host_uptake=H(-psi)*host_activity
				*theta*alpha*N_ATP_PER_GLUCOSE/6;
	return m[CO2](i,j)*uptake_per_cell(m,i,j) + host_uptake;
}

template<class fid_t>
void
eval_nutrient_fill_sle_matrix(const AMesh2D<fid_t>& m, const BCSet& bcs,
	ASparseMatrix& A, vector<double>& rhs, MeshEnumerator const& k)
throw(MeshException) {
	int xdim=m.get_xdim();
	int ydim=m.get_ydim();
	int i, j;
	try {
		array2d c_arr=m[CO2];
		array2d phi_arr=m[PHI];
		array2d psi_arr=m[PSI];
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
	// should use bc.c()/bc.a() for appropriate bc instead of m.get(CO2)
	// from the previous time step
				if (k0p >= 0) {
					A.set(k0,k0p, dy2f);
				} else {
					rhs.at(k0)=-dy2f*c_arr(i,j+1);
				}
				if (k0m >= 0) {
					A.set(k0,k0m, dy2f);
				} else {
					rhs.at(k0)=-dy2f*c_arr(i,j-1);
				}
				if (kp0 >= 0) {
					A.set(k0,kp0, dx2f);
				} else {
					rhs.at(k0)=-dx2f*c_arr(i+1,j);
				}
				if (km0 >= 0) {
					A.set(k0,km0, dx2f);
				} else {
					rhs.at(k0)=-dx2f*c_arr(i-1,j);
				}
				// add consumption term
				A.set(k0,k0,A.get(k0,k0)
				-uptake_per_cell(m,phi_arr,psi_arr,i,j));
			}
		}
	} catch (SparseMatrixException& e) {
		ostringstream ss;
		ss << "eval_nutrient_fill_sle_matrix: "
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
		ss << "eval_nutrient_fill_sle_matrix: "
			<< "failed to construct boundary points equations: "
			<< e.what();
		throw MeshException(ss.str());
	}
}

template<class fid_t>
AMesh2D<fid_t>*
eval_nutrient_directly(const Params& p, const AMesh2D<fid_t>& m1)
throw(MeshException) {
	try {
		SkipDirichletEnumerator<fid_t> kenum(m1,p.c_bc);
		auto_ptr<AMesh2D<fid_t> > m2(m1.clone());
		// initial estimate (previous solution)
		vector<double> c(kenum.size());
		array2d c_arr=(*m2)[CO2];
		for (int k=0; k<kenum.size(); ++k) {
			c.at(k)=c_arr(kenum.i(k),kenum.j(k));
		}
		auto_ptr<ASparseMatrix> pA(build_sle_solver_matrix(
			kenum.size(), c, Method::it().p_solver_accuracy));
		vector<double> rhs(kenum.size());// right hand side vector
		eval_nutrient_fill_sle_matrix(*m2,p.c_bc,*pA,rhs,kenum);
		// solve SLE
		try {
			c=pA->solve(rhs);
		} catch (SparseMatrixException& e) {
			ostringstream ss;
			ss << "SparseMatrixException solve(): " << e.what();
			throw MeshException(ss.str());
		}
		// update m2
		for (int k=0; k<(int)c.size(); ++k) {
			int i=kenum.i(k);
			int j=kenum.j(k);
			m2->set(CO2,i,j,c.at(k));
		}
		// update Dirichlet boundary points (excluded from SLE)
		update_dirichlet_points<fid_t>(*m2,p.c_bc,CO2);
		return m2.release();
	} catch (MeshException& e) {
		ostringstream ss;
		ss << "eval_nutrient_directly: " << e.what();
		throw MeshException(ss.str());
	}
}

template<class fid_t>
double
eval_nutrient_residual(const AMesh2D<fid_t>& m)
throw(MeshException) {
	double global_res=0.0;
	double dx=m.get_dx();
	double dy=m.get_dy();
	const auto_ptr<AMesh2D<fid_t> > m2(m.clone());
	m2->remove_function_ifdef(C_RESIDUAL);
	m2->add_function(C_RESIDUAL);
	// WARNING: ignoring boundary points
	// inner points
	array2d c=(*m2)[CO2];
	for (int i=1; (i<m2->get_xdim()-1) ; ++i) {
		for (int j=1; (j<m2->get_ydim()-1) ; ++j) {
			double loc_res=0.0;
			loc_res=(c(i+1,j)-2*c(i,j)+c(i-1,j))/(dx*dx)+
				(c(i,j+1)-2*c(i,j)+c(i,j-1))/(dy*dy)-
				uptake_term(*m2,i,j);
			loc_res=fabs(loc_res);
			m2->set(C_RESIDUAL,i,j,loc_res);
		}
	}
	global_res=norm_1<fid_t>(*m2,C_RESIDUAL);
	return global_res;
}

template<class fid_t>
void
eval_nutrient_init_consumption(AMesh2D<fid_t>& m, fid_t const consumption) {
	m.add_function_ifndef(consumption);
	// initialize consumption variable
	int xdim=m.get_xdim();
	int ydim=m.get_ydim();
	for (int i=0; i<xdim ; ++i) {
		for (int j=0; j<ydim ; ++j) {
			m.set(consumption,i,j,uptake_term(m,i,j));
		}
	}
}

template<class fid_t>
AMesh2D<fid_t>*
eval_nutrient_iteratively(const Params& p, const AMesh2D<fid_t>& m,
		double const epsilon)
throw(MeshException) {
	try {
	auto_ptr< AMesh2D<fid_t> > m2(m.clone());
	double dt;
	if (Method::it().p_solver_relax_step < 0) { // automatic time step
		if (Method::it().rd_solver == MP::RDS_EXPLICIT) {
			dt=best_dt(*m2);
		} else { // implicit
			double dx=m.get_dx();
			double dy=m.get_dy();
			dt=10.0*(dx*dx+dy*dy); // heuristic
		}
	} else {
		dt=Method::it().p_solver_relax_step;
	}
	if (verbose > 1) { // extra verbose
		cerr << dbg_stamp(m.get_time())
			<< "eval_nutrient_iteratively: dt=" << dt << "\n";
	}
	// iterative method here
	double res=0.0;
	res=eval_nutrient_residual<fid_t>(*m2);
	int iter=0;
	while ((res>epsilon) && (iter<=Method::it().p_solver_relax_max_iters)) {
		eval_nutrient_init_consumption<fid_t>(*m2,Q_TMP);
		switch (Method::it().rd_solver) {
		case MP::RDS_EXPLICIT:
			m2.reset(reaction_diffusion_step<fid_t>(p.c_bc,dt,*m2,
				CO2,NONE,Q_TMP,MP::RDS_EXPLICIT));
			break;
		case MP::RDS_IMPLICIT:
			m2.reset(reaction_diffusion_step<fid_t>(p.c_bc,dt,*m2,
				CO2,NONE,Q_TMP,MP::RDS_IMPLICIT));
			break;
		case MP::RDS_ADI:
			m2.reset(reaction_diffusion_step<fid_t>(p.c_bc,dt,*m2,
				CO2,NONE,Q_TMP,MP::RDS_ADI));
			break;
		default:
			throw MeshException("eval_nutrient_iteratively: "
				"unknown method");
			break;
		}
		m2->remove_function_ifdef(Q_TMP);
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

template<class fid_t>
AMesh2D<fid_t>*
eval_nutrient_diffusion
(Params const& p, double const dt, AMesh2D<fid_t> const& m1)
throw(MeshException) {
	auto_ptr< AMesh2D<fid_t> > m2(m1.clone());
	eval_nutrient_init_consumption<fid_t>(*m2,Q_TMP);
	m2->remove_function_ifdef(D_C_TMP);
	m2->add_function(D_C_TMP,1.0);
	m2.reset(reaction_diffusion_step<fid_t>(p.c_bc,dt,*m2,CO2,D_C_TMP,Q_TMP));
	m2->remove_function_ifdef(Q_TMP);
	m2->remove_function_ifdef(D_C_TMP);
	return m2.release();
}

template<class fid_t>
AMesh2D<fid_t>*
eval_nutrient_poisson
(Params const& p, AMesh2D<fid_t> const& m1, double const epsilon)
throw(MeshException) {
	AMesh2D<fid_t> *m2;
	if (Method::it().p_solver == MP::PS_SLE) {
		m2=eval_nutrient_directly<fid_t>(p, m1);
	} else if (Method::it().p_solver == MP::PS_RELAX) {
		m2=eval_nutrient_iteratively<fid_t>(p, m1, epsilon);
	} else {
		throw MeshException("eval_nutrient: unknown method");
	}
	return m2;
}

template<class fid_t>
AMesh2D<fid_t>*
eval_nutrient(const Params& p, const AMesh2D<fid_t>& m1,
	double const epsilon, double const dt)
throw(MeshException) {
	if (!m1.defined(CO2)) {
		throw MeshException("eval_nutrient: CO2 not defined");
	}
	if (!m1.defined(PHI)) {
		throw MeshException("eval_nutrient: PHI not defined");
	}
	if (!m1.defined(PSI)) {
		throw MeshException("eval_nutrient: PSI not defined");
	}
	if (!m1.attr_defined("o2_uptake")) {
		throw MeshException("eval_nutrient: "
			"oxygen consumption rate not defined");
	}
	if (!m1.attr_defined("host_activity")) {
		throw MeshException("eval_nutrient: "
			"host behaviour (host_activity) not defined");
	}
	try {
		if (p.c_equation == Params::EQ_POISSON) {
			return eval_nutrient_poisson<fid_t>(p,m1,epsilon);
		} else if (p.c_equation == Params::EQ_DIFFUSION) {
			return eval_nutrient_diffusion<fid_t>(p,dt,m1);
		} else {
			throw MeshException("eval_nutrient: "
				"unsupported nutrient equation");
		}
	} catch (MeshException& e) {
		ostringstream ss;
		ss << "eval_nutrient: " << e.what();
		throw MeshException(ss.str());
	}
}

// templates
template
AMesh2D<int>*
eval_nutrient<int>(const Params& p, const AMesh2D<int>& m1,
	double const epsilon, double const dt);


