/*
 * Copyright (C) 2006-2007 Sergey Astanin, Luigi Preziosi, Andrea Tosin
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

#ifndef SPMAT_H
#define SPMAT_H

#include "config.h"

#include <vector>
#include <iostream>
#include <string>
#include <iostream>
#include <sstream>
#include <exception>

extern "C" {
#if HAVE_UFSPARSE_UMFPACK_H
#include <ufsparse/umfpack.h>
#else
#if HAVE_UMFPACK_UMFPACK_H
#include <umfpack/umfpack.h>
#endif
#endif
}

class SparseMatrixException : public std::exception {
	std::string _what;
	int _code;
public:
	SparseMatrixException(std::string const& msg) :
		_what(msg), _code(-1) { }
	SparseMatrixException(std::string const& msg, int const code) :
		_code(code) {
		std::ostringstream ss;
		ss << msg << "; error " << code;
		_what=ss.str();
	}
	~SparseMatrixException() throw() {}
	const char* what() const throw() { return _what.c_str(); }
	int code() const throw() { return _code; }
};

class ASparseMatrix {
public:
	virtual ~ASparseMatrix() {};
	// "virtual" constructors
	virtual ASparseMatrix* clone(void) const = 0;
	// matrix data manipulation
	virtual int get_rows(void) const = 0;
	virtual int get_cols(void) const = 0;
	virtual void set(int const i, int const j, double const value) = 0;
	virtual double get(int const i, int const j) const = 0;
	/// solver A*x=rhs, where A is the matrix; returns x
	virtual std::vector<double>
	solve(std::vector<double>& rhs) throw(SparseMatrixException) = 0;
};

#ifdef HAVE_LIBUMFPACK
class UMFPACKMatrix : public ASparseMatrix {
private:
	/// unidimensional index bijected on (i,j) index-space
	int k(int const i, int const j) const {
		return (i+j*n_cols);
	}

	int n_rows;
	int n_cols;
	std::vector<int> Ti; ///< triplets' first index
	std::vector<int> Tj; ///< triplets' second index
	std::vector<double> Tx; ///< triplets' value;
	/** n_rows*n_cols vector, describes matrix pattern, -1 -- empty cell,
	    otherwise, index in triplet form */
	std::vector<int> pattern;

public:
	UMFPACKMatrix(int const _rows, int const _cols) :
		n_rows(_rows), n_cols(_cols), Ti(0), Tj(0), Tx(0),
		pattern(_rows*_cols) {
		std::vector<int>::iterator i;
		for (i=pattern.begin(); i!=pattern.end(); ++i) {
			*i=-1; // empty cell;
		}
	}

	virtual UMFPACKMatrix* clone(void) const {
		return new UMFPACKMatrix(*this);
	}

	~UMFPACKMatrix() {}

	virtual int get_rows(void) const {
		return n_rows;
	}

	virtual int get_cols(void) const {
		return n_cols;
	}

	virtual void set(int const i, int const j, double const value) {
		if ((i < 0) || (i >= n_rows) || (j < 0) || (j >= n_cols)) {
			std::ostringstream ss;
			ss << "UMFPACKMatrix::set: indices i=" << i
				<< ", j=" << j << " are out of range ("
				<< n_rows << "," << n_cols << ")";
			throw SparseMatrixException(ss.str(), -1);
		}
		// do change triplet, if it is already there
		if (pattern.at(k(i,j)) != -1) {
			int Tx_index=pattern.at(k(i,j));
			Tx.at(Tx_index)=value;
		} else { // add new triplet
			Ti.push_back(i);
			Tj.push_back(j);
			Tx.push_back(value);
			// pattern contains indices of triplet arrays
			pattern.at(k(i,j))=Tx.size()-1;
		}
	}

	virtual double get(int const i, int const j) const {
		int Tx_index=pattern.at(k(i,j));
		if (Tx_index == -1) {
			return 0.0; // empty cell
		} else {
			return Tx.at(Tx_index);
		}
	}

	virtual std::vector<double>
	solve(std::vector<double>& rhs) throw(SparseMatrixException) {
		if (rhs.size() != (size_t) n_rows) {
			throw SparseMatrixException("solve: inappropriate "
				"length of right hand side");
		}
		// convert to compressed-column form
		int nz=Tx.size();
		int *Ap=new int[n_cols+1];
		int *Ai=new int[nz];
		double *Ax=new double[nz];
		if ((Ap == (int*)0) || (Ai == (int*)0) || (Ax == (double*)0)) {
			if (Ap) {
				delete[] Ap; Ap=0;
			}
			if (Ai) {
				delete[] Ai; Ai=0;
			}
			if (Ax) {
				delete[] Ax; Ax=0;
			}
			throw SparseMatrixException("solve: cannot alloc "
				"memory for compressed-column form");
		}
		int status=umfpack_di_triplet_to_col(n_rows, n_cols, nz,
				&*Ti.begin(), &*Tj.begin(), &*Tx.begin(),
				Ap, Ai, Ax, (int*)0);
		if (status) {
			delete[] Ap; Ap=0;
			delete[] Ai; Ai=0;
			delete[] Ax; Ax=0;
			throw SparseMatrixException("solve: triplet_to_col "
				"failed", status);
		}
		// start solving
		double *px=new double[n_rows];
		if (px == (double*)0) {
			delete[] Ap; Ap=0;
			delete[] Ai; Ai=0;
			delete[] Ax; Ax=0;
			throw SparseMatrixException("solve: cannot alloc "
				"solution vector", status);
		}
		void *symbolic, *numeric;
		status=umfpack_di_symbolic(n_rows, n_cols, Ap, Ai, Ax,
				&symbolic, (double*)0, (double*)0);
		if (status) {
			delete[] Ap; Ap=0;
			delete[] Ai; Ai=0;
			delete[] Ax; Ax=0;
			delete[] px; px=0;
			throw SparseMatrixException("solve: symbolic failed",
					status);
		}
		status=umfpack_di_numeric(Ap, Ai, Ax, symbolic, &numeric,
				(double*)0, (double*)0);
		if (status) {
			delete[] Ap; Ap=0;
			delete[] Ai; Ai=0;
			delete[] Ax; Ax=0;
			delete[] px; px=0;
			throw SparseMatrixException("solve: numeric failed",
					status);
		}
		umfpack_di_free_symbolic(&symbolic);
		status=umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, px,
				&*rhs.begin(), numeric, (double*)0, (double*)0);
		umfpack_di_free_numeric(&numeric);
		if (status) {
			delete[] Ap; Ap=0;
			delete[] Ai; Ai=0;
			delete[] Ax; Ax=0;
			delete[] px; px=0;
			throw SparseMatrixException("solve: solve failed",
					status);
		}
		std::vector<double> x(px,px+n_rows);
		// clean up
		delete[] Ap; Ap=0;
		delete[] Ai; Ai=0;
		delete[] Ax; Ax=0;
		delete[] px; px=0;
		return x;
	}

};
#endif // if UMFPACK is present

#include "lsolver/matrix.h"
#include "lsolver/cghs.h"
#include "lsolver/bicgsq.h"
#include "lsolver/bicgstab.h"
#include "lsolver/gmres.h"

class LSolverMatrix : public ASparseMatrix, public ALSolverMatrix {
public:
	typedef enum {
		CG=0,
		BICG=1,
		BICGstab=2,
		GMRES=3
	} lsolver_method;
private:
	int k(int const i, int const j) const { return i*n+j; };
	int n;
	std::vector<double> x;
	double accuracy;
	int max_iterations;
	std::vector<int> Ti; ///< triplets' first index
	std::vector<int> Tj; ///< triplets' second index
	std::vector<double> Ta; ///< triplets' value;
	lsolver_method method;
	int gmres_restart;
	int find_entry(int const i, int const j) const {
		int n=Ti.size();
		for (int k=0; k<n; ++k) {
			if ((Ti[k]==i) && (Tj[k]==j)) {
				return k;
			}
		}
		return -1;
	}
public:
	LSolverMatrix(int const rows, const std::vector<double>& x0,
		double const accuracy, int const max_iterations,
		lsolver_method const method=BICGstab, int gmres_restart=-1)
		: n(rows), x(x0), accuracy(accuracy),
		max_iterations(max_iterations), method(method),
		gmres_restart(gmres_restart) {
	}
	virtual ~LSolverMatrix() {
	}
	virtual LSolverMatrix* clone(void) const {
		return new LSolverMatrix(*this);
	}
	virtual int get_rows() const { return n; }
	virtual int get_cols() const { return n; }
	virtual void set(int const i, int const j, double const value) {
		if ((i < 0) || (i >= n) || (j < 0) || (j >= n)) {
			std::ostringstream ss;
			ss << "LSolverMatrix::set: indices i=" << i
				<< ", j=" << j << " are out of range ("
				<< n << "," << n << ")";
			throw SparseMatrixException(ss.str(), -1);
		}
		// WARNING: DOES NOT CHECK IF THE TRIPLET EXISTS
	//	// do change triplet, if it is already there
	//	int k=find_entry(i,j);
	//	if (k != -1) {
	//		Ta.at(k)=value;
	//	} else { // add new triplet
			Ti.push_back(i);
			Tj.push_back(j);
			Ta.push_back(value);
	//	}
	}
	virtual double get(int const i, int const j) const {
		int Ta_index=find_entry(i,j);
		if (Ta_index == -1) {
			return 0.0; // empty cell
		} else {
			return Ta.at(Ta_index);
		}
	}
	virtual std::vector<double>
	solve(std::vector<double>& rhs) throw(SparseMatrixException) {
		if (rhs.size() != (size_t) n) {
			throw SparseMatrixException("LSolverMatrix::solve: "
				"inappropriate length of right hand side");
		}
		// WARNING: vector<double> x encapsulation is broken here
		// x is modified within iterative solver
		int iters=-1;
		switch (method) {
		case CG:
			iters=cghs(n,*this,rhs.empty()?NULL:&rhs[0],
				x.empty()?NULL:&x[0],accuracy);
			break;
		case BICG:
			iters=bicgsq(n,*this,rhs.empty()?NULL:&rhs[0],
				x.empty()?NULL:&x[0],accuracy);
			break;
		case BICGstab:
			iters=bicgstab(n,*this,rhs.empty()?NULL:&rhs[0],
				x.empty()?NULL:&x[0],accuracy);
			break;
		case GMRES:
			iters=gmres(gmres_restart>0?gmres_restart:(1+n/10),
				n,*this, rhs.empty()?NULL:&rhs[0],
				x.empty()?NULL:&x[0],accuracy);
			break;
		default:
			std::ostringstream ss;
			ss << "LSolverMatrix::solve: unknown method ("
				<< static_cast<int>(method) << ")";
			throw SparseMatrixException(ss.str());
			break;
		}
		return x;
	}
	// LSolver specific methods
	virtual void
	mult(const double* v, double *w) const {
		// reset output
		for (int i=0; i<n; ++i) {
			*(w+i)=0.0;
		}
		int size=Ta.size();
		for (int k=0; k<size; ++k) {
			int i=Ti[k];
			int j=Tj[k];
			double aij=Ta[k];
			*(w+i)=*(w+i)+aij*v[j];
		}
	}
	virtual void
	set_gmres_restart(int const m) {
		gmres_restart=m;
	}
};

#ifdef HAVE_LIBLAPACK
int
dgtsv(int N, int NRHS, double *DL, double *D, double *DU, double *B, int ldb);

/// class for general tridiagonal matrix (with diagonal domination)
class TridiagMatrix : public ASparseMatrix {
	int n;
	std::vector<double> u;
	std::vector<double> d;
	std::vector<double> l;
public:
	TridiagMatrix(int const n) : n(n), u(n-1), d(n), l(n-1) {}
	virtual ~TridiagMatrix() {};
	virtual TridiagMatrix* clone(void) const {
		return new TridiagMatrix(*this);
	}
	virtual int get_rows(void) const { return n; }
	virtual int get_cols(void) const { return n; }
	virtual void set(int const i, int const j, double const value) {
		if ((i < 0) || (i >= n) || (j < 0) || (j >= n)) {
			throw SparseMatrixException("TridiagMatrix::set: "
					"index out of range");
		}
		int offset=j-i;
		if (abs(offset) > 1) {
			throw SparseMatrixException("TridiagMatrix::set: "
					"access out of tridiagonal band");
		}
		switch(offset) {
		case 0:
			d.at(i)=value; break;
		case 1:
			u.at(i)=value; break;
		case -1:
			l.at(i-1)=value; break;
		default:
			throw SparseMatrixException("TridiagMatrix::set: "
				"invalid offset");
			break;
		}
	}
	virtual double get(int const i, int const j) const {
		if ((i < 0) || (i >= n) || (j < 0) || (j >= n)) {
			throw SparseMatrixException("TridiagMatrix::get: "
					"index out of range");
		}
		int offset=j-i;
		if ((offset>1)||(offset<-1)) { // off the three main diagonals
			return 0.0;
		} else {
			switch(offset) {
			case 0:
				return d.at(i); break;
			case 1:
				return u.at(i); break;
			case -1:
				return l.at(i-1); break;
			default:
				throw SparseMatrixException("TridiagMatrix::get"
					": invalid offset");
				break;
			}
		}
	}
	virtual std::vector<double>
	solve(std::vector<double>& rhs) throw(SparseMatrixException) {
		if (rhs.size() != static_cast<size_t>(n)) {
			throw SparseMatrixException("TridiagMatrix::solve: "
					"wrong dimension of rhs");
		}
		std::vector<double> x(n);
		x=rhs;
		int ret;
		ret=dgtsv(n,1,l.empty()?NULL:&l[0],d.empty()?NULL:&d[0],
			u.empty()?NULL:&u[0],x.empty()?NULL:&x[0],n);
		if (ret) {
			throw SparseMatrixException("TridiagMatrix::solve: "
				"error in dgtsv (LAPACK) routing");
		}
		return x;
	}
};
#endif

bool
is_symmetric(ASparseMatrix const& A);

#endif /* SPMAT_H */
