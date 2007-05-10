#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "function.h"
#include "newton.h"

// link like this
// g++ -Wall -pedantic -ansi -I../ -lm -o test_newton test_newton.cc

using namespace std;

int main(int argc, char *argv[]) {
	int max_iter=100;
	double e=1e-10;
	srand(time(0)%1000);
	double exact=10.0*(rand()/(RAND_MAX+1.0))-5.0; // in [-5.0;+5.0)
	double x0=10.0*(rand()/(RAND_MAX+1.0))-5.0; // in [-5.0;+5.0)
	SigmaFunction f(exact,1.0,1.0);
	SigmaPrimeFunction df(exact,1.0,1.0);
	NewtonRootFinder solver((ADoubleFunction&)f,(ADoubleFunction&)df,x0,e);
	int i=0;
	try {
		while ((i<max_iter) && (!solver.is_found())) {
			solver.iterate();
			++i;
		}
		if (solver.is_found()) {
			if (fabs(solver.get_root()-exact)>=e) {
				cerr << "Root finding error exceeds " << e
					<< "\nFAILURE\n";
				return 1;
			} else {
				cerr << "root=" << exact
					<< " found=" << solver.get_root()<<"\n";
				cerr << "OK\n";
				return 0;
			}
		} else {
			cerr << "Root not found in " << i << "iterations\n";
			cerr << "FAILURE\n";
			return 1;
		}
	} catch (FunctionException& e) {
		cerr << "Exception raised: " << e.what() << "\n";
		cerr << "FAILURE\n";
		return 1;
	}
}
