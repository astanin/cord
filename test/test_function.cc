#include "function.h"
#include <iostream>
#include <math.h>
#include <time.h>

using namespace std;

// link like this:
// g++ -I../ -lm -o test_function test_function.cc

double c_twox(double const arg) {
	return 2.0*arg;
}

int main(int argc, char* argv[]) {
	double arg0=5.0;
	double k1=1.0;
	double s1=0.0;
	SigmaFunction sigma(arg0,k1,s1);
	SigmaPrimeFunction sigmaprime(arg0,k1,s1);
	double x1=1e-10+10.0*(rand()/(RAND_MAX+1.0)); // x1 in (0;10.0)
	double x2=1e-10+10.0*(rand()/(RAND_MAX+1.0)); // x2 in (0;10.0)
	if (x2 < x1) { // swap
		double tmp=x1;
		x1=x2;
		x2=tmp;
	}
	if (sigma(x1) >= sigma(x2)) {
		cerr << "sigma("<<x1<<")="<<sigma(x1)<<" >= "
			<< "sigma("<<x2<<")="<<sigma(x2)<<"\n";
		cerr << "Sigma is not monotonous\n";
		return 1;
	} else {
		cerr << "Sigma appears to monotonous\n";
	}
	if (( fabs((sigma(arg0+x1)-sigma(arg0))/(x1)-
		(sigma(arg0+x2)-sigma(arg0))/(x2)) > 1e-99 ) ||
		( fabs((sigma(arg0-x1)-sigma(arg0))/(-x1)-
		(sigma(arg0-x2)-sigma(arg0))/(-x2)) > 1e-99 )) {
		cerr << "x1=" << x1 <<  "(sigma(arg0+x1)-sigma(arg0))/(x1)="
			<< (sigma(arg0+x1)-sigma(arg0))/(x1) << "\n";
		cerr << "x2=" << x2 << "(sigma(arg0+x2)-sigma(arg0))/(x2)="
			<< (sigma(arg0+x2)-sigma(arg0))/(x2) << "\n";
		cerr << "x1=" << x1 <<  "(sigma(arg0-x1)-sigma(arg0))/(-x1)="
			<< (sigma(arg0-x1)-sigma(arg0))/(-x1) << "\n";
		cerr << "x2=" << x2 << "(sigma(arg0-x2)-sigma(arg0))/(-x2)="
			<< (sigma(arg0-x2)-sigma(arg0))/(-x2) << "\n";
		cerr << "Sigma is not piece linear\n";
		return 1;
	} else {
		cerr << "Sigma appears to be piece linear\n";
	}
	if ( (fabs(sigmaprime(arg0+x1)-sigmaprime(arg0+x2))+
		fabs(sigmaprime(arg0-x1)-sigmaprime(arg0-x2))) > 1e-99 ) {
		cerr << "SigmaPrime is not piece constant\n";
		return 1;
	} else {
		cerr << "SigmaPrime appears to be piece constant\n";
	}
	SumFunction sum((ADoubleFunction*)&sigmaprime,(ADoubleFunction*)&sigma);
	if (fabs(sum(x1)-sigmaprime(x1)-sigma(x1)) > 1e-99) {
		cerr << "SumFunction(ptrs)!=sigmaprime+sigma\n";
		return 1;
	} else {
		cerr << "SumFunction(ptrs) returns the sum\n";
	}
	SumFunction sumr((ADoubleFunction&)sigmaprime,(ADoubleFunction&)sigma);
	if (fabs(sumr(x1)-sigmaprime(x1)-sigma(x1)) > 1e-99) {
		cerr << "SumFunction(refs)!=sigmaprime+sigma\n";
		return 1;
	} else {
		cerr << "SumFunction(refs) returns the sum\n";
	}
	ProductFunction p((ADoubleFunction*)&sigmaprime,
				(ADoubleFunction*)&sigma);
	if (fabs(p(x1)-sigmaprime(x1)*sigma(x1)) > 1e-99) {
		cerr << "ProductFunction(sigmaprime,sigma)!=sigmaprime*sigma\n";
		return 1;
	} else {
		cerr << "ProductFunction returns the product\n";
	}
	ConstFunction cf(arg0);
	if (fabs(cf(x1)-cf(x2)) > 1e-99) {
		cerr << "ConstFunction is not constant\n";
		return 1;
	} else {
		cerr << "ConstFunction appears to be constant\n";
	}
	cerr << "OK\n";
	ConstFunction two(2.0);
	LinearFunction x(1.0);
	ProductFunction twox((ADoubleFunction*)&two,(ADoubleFunction*)&x);
	int i;
	double r=0.0;
	int N=100000000;
	clock_t start, end;
	start=clock();
	for (i=0; i<N; ++i) {
		r=twox(3.1416);
	}
	end=clock();
	cerr << "ProductFunction: " << N << " calls made in " << (end-start)*1.0/CLOCKS_PER_SEC << " sec\n";
	start=clock();
	for (i=0; i<N; ++i) {
		r=c_twox(3.1416);
	}
	end=clock();
	cerr << "C function calls " << N << " calls made in " << (end-start)*1.0/CLOCKS_PER_SEC << " sec\n";
	return 0;
}

