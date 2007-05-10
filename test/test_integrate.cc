#include <iostream>
#include <gsl/gsl_integration.h>
#include <cmath>

// link like this:
// g++ -o test_integrate -L/usr/lib/sse2 -lcblas -lgsl -lm test_integrate.cc

using namespace std;

double piece_linear(const double arg, void *params) {
	if (arg < 1.0) {
		return 0.0;
	} else {
		return (arg-1.0);
	}
}

int main(int argc, char *argv[]) {
	double result=-1.0;
	double error=-1.0;
	size_t limit=10;
	gsl_integration_workspace *w=gsl_integration_workspace_alloc(limit);
	gsl_function F;
	F.function=&piece_linear;
	F.params=0;
	gsl_integration_qag(&F, 0.0, 2.0, 1e-6, 1e-3, limit, GSL_INTEG_GAUSS15,
		w, &result, &error); 
	gsl_integration_workspace_free(w);
	cerr << "\\int_0^2f(x)dx= " << result << " error= " << error << "\n";
	if (fabs(result-0.5)>1e-6) {
		cerr << "FAILURE\n";
		return 1;
	} else {
		cerr << "OK\n";
		return 0;
	}
}

