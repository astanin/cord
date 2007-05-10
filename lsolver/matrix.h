#ifndef ALSOLVERMATRIX_H
#define ALSOLVERMATRIX_H

/// example matrix class interface to be used with LSolver routines
class ALSolverMatrix {
public:
	virtual ~ALSolverMatrix() {}

	/// matrix-vector multiplication (product with v, save result in w)
	virtual void
	mult(const double *v, double *w) const = 0;

	/// set restart iteration count of the GMRES method
	virtual void
	set_gmres_restart(int const gmres_restart) = 0;
};

#endif

