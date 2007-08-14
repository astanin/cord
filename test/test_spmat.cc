#include "spmat.h"
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <limits>

#include <cstdlib>
#include <time.h>

// build line:
// g++ -o test_spmat -DHAVE_LIBUMFPACK -DHAVE_LIBLAPACK -I.. test_spmat.cc ../spmat.cc -lumfpack -lamd -lgslcblas -llapack -lblas

#ifdef DHAVE_CONFIG_H
#include "config.h"
#endif

using std::floor;
using std::fabs;

using namespace std;

bool is_solution(ASparseMatrix &A, vector<double>& x, vector<double>& rhs) {
	int N=x.size();
	cout.setf(ios::right | ios::fixed);
	cout.precision(2);
	double err;
	double sum;
	double xnorm=numeric_limits<double>::epsilon();
	double enorm=numeric_limits<double>::epsilon();
	for (int i=0; i<N; ++i) {
		sum=0.0;
//		cout << "|";
		for (int j=0; j<N; ++j) {
//			cout.width(6);
//			cout << A.get(i,j) << " ";
			sum+=A.get(i,j)*x[j];
		}
//		cout << "|";
//		cout << " |";
//		cout.width(6);
//		cout << x[i];
//		cout << "|";
//		if (i==N/2) {
//			cout << " = ";
//		} else {
//			cout << "   ";
//		}
//		cout << "|";
//		cout.width(6);
//		cout << rhs[i];
//		cout << "|";
		sum-=rhs[i];
//		cout << "\n";
		xnorm+=fabs(x[i]);
		enorm+=fabs(sum);
	}
	cout << "error: "<<enorm<<"; rel.error: "<<enorm/xnorm << "; ";
	if (enorm < 0.01*xnorm) {
		cout << "OK\n";
		return true;
	} else {
		cout << "FAILED\n";
		return false;
	}
}

int test_mult(const int N) {
	std::vector<double> x(N);
	std::vector<double> y(N);
	LSolverMatrix A(N,x,1e-5,1000,LSolverMatrix::CG);
	// identity matrix multiplication
	for (int i=0; i<N; ++i) {
		x[i]=1.0;
		A.set(i,i,2.0);
	}
	A.mult(x.empty()?NULL:&x[0],y.empty()?NULL:&y[0]);
	for (int i=0; i<N; ++i) {
		if (fabs(y[i]-2.0) > numeric_limits<double>::epsilon()) {
			cout << "2 I[][] * ones[]  != twos[], i=" << i << "\n";
			cout << "FAILED\n";
			return 1;
		}
	}
	// tridiagonal matrix multiplication now
	A.set(0,1,-2.0);
	A.set(N-1,N-2,-2.0);
	for (int i=1; i<(N-1); ++i) {
		A.set(i,i-1,-1.0);
		A.set(i,i+1,-1.0);
	}
	A.mult(x.empty()?NULL:&x[0],y.empty()?NULL:&y[0]);
	for (int i=0; i<N; ++i) {
		if (fabs(y[i]) > numeric_limits<double>::epsilon()) {
			cout << "D[][] * ones[]  != zeros[], "
				"y[" << i << "]=" << y[i] << "\n";
			cout << "FAILED\n";
			return 1;
		}
	}
	cout << "OK\n";
	return 0;
}

int test_set_get(ASparseMatrix &A,const int N) {
	for (int i=0; i<N; ++i) {
		for (int j=0; j<N; ++j) {
			A.set(i,j,(i+j));
			double g=A.get(i,j);
			if (fabs(g-(i+j)*1.0) >
				numeric_limits<double>::epsilon()) {
				cout << "set at (" << i << "," << j << "): "
					<< (i+j) << "; found: " << g << "\n";
				cout << "FAILED\n";
				return 1;
			}
		}
	}
	// once again
	for (int i=0; i<N; ++i) {
		for (int j=0; j<N; ++j) {
			double g=A.get(i,j);
			if (fabs(g-(i+j)*1.0) >
				numeric_limits<double>::epsilon()) {
				cout << "at (" << i << "," << j << "): "
					<< "found: " << g << " instead of "
					<< (i+j) << "\n";
				cout << "FAILED\n";
				return 1;
			}
		}
	}
	cout << "OK\n";
	return 0;
}

int test(ASparseMatrix &A,const int N, const double DIAG) {
	vector<double> rhs(N);
	A.set(0,0,1.0);
	rhs[0]=0.0;
	A.set(N-1,N-1,1.0);
	rhs[N-1]=0.0;
	for (int i=1; i<(N-1); ++i) {
		A.set(i,i,DIAG);
		rhs[i]=i*(N-1-i)/(N*N*0.25);
		for (int j=0; j<N; ++j) {
			if ((abs(i-j)<=2)&&(i!=j)) {
				double a=rand()%((int)floor(DIAG));
				A.set(i,j,a);
			}
		}
	}
	vector<double> x=A.solve(rhs);
	bool ok=is_solution(A,x,rhs);
	if (ok) {
		return 0;
	} else {
		return 1;
	}
}

#ifdef HAVE_LIBLAPACK
int td_test(TridiagMatrix& T, int const N) {
	double lb=1.0;
	double rb=N*1.0;
	vector<double> rhs(N);
	rhs[0]=lb;
	rhs[N-1]=rb;
	T.set(0,0,1.0);
	T.set(N-1,N-1,1.0);
	for (int i=1; i<(N-1); ++i) {
		T.set(i,i-1,-1);
		T.set(i,i,2);
		T.set(i,i+1,-1);
	}
	vector<double> x=rhs;
	x=T.clone()->solve(x);
	bool ok=is_solution(T,x,rhs);
	if (ok) {
		return 0;
	} else {
		return 1;
	}
}
#endif

int is_symm_test(int const N) {
	int ret=0;
	std::vector<double> x0(N);
	LSolverMatrix Ta(N,x0, 0.0, 0);
	Ta.set(0,0,1.0);
	Ta.set(N-1,N-1,1.0);
	for (int i=1; i<(N-1); ++i) {
		Ta.set(i,i-1,-1);
		Ta.set(i,i,2);
		Ta.set(i,i+1,-1);
	}
	cout << "Asymmetry test: ";
	if (is_symmetric(Ta)) {
		cout << "FAILED\n";
		ret++;
	} else {
		cout << "OK\n";
	}
	LSolverMatrix Tb(N, x0, 0.0, 0);
	Tb.set(0,0, 1);
	Tb.set(0,1,-1);
	Tb.set(N-1,N-1, 1);
	Tb.set(N-1,N-2,-1);
	for (int i=1; i<(N-1); ++i) {
		Tb.set(i,i-1,-1);
		Tb.set(i,i,2);
		Tb.set(i,i+1,-1);
	}
	cout << "Symmetry test: ";
	if (!is_symmetric(Tb)) {
		cout << "FAILED\n";
		ret++;
	} else {
		cout << "OK\n";
	}
	return ret;
}

int main(int argc, char *argv[]) {
	int ret=0;
	const int N=100;
	const double DIAG=10;
	srand(time(0));
	std::vector<double> x(N);
	{
		LSolverMatrix A(N,x,1e-5,1000,LSolverMatrix::GMRES,N/2);
		cout << "LSolverMatrix set-get test: ";
		ret+=test_set_get(A,N);
	}
	{
		cout << "LSolverMatrix mult test: ";
		ret+=test_mult(N);
	}
	{
		LSolverMatrix A(N,x,1e-5,1000,LSolverMatrix::BICG);
		cout << "LSolverMatrix (BiCG) solve test: ";
		ret+=test(A,N,DIAG);
	}
	{
		LSolverMatrix A(N,x,1e-5,1000,LSolverMatrix::GMRES,N);
		cout << "LSolverMatrix (GMRES) solve test: ";
		ret+=test(A,N,DIAG);
	}
#ifdef HAVE_LIBUMFPACK
	{
		cout << "UMFPACKMatrix set-get test: ";
		UMFPACKMatrix C(N,N);
		ret+=test_set_get(C,N);
	}
	{
		cout << "UMFPACKMatrix solve test: ";
		UMFPACKMatrix C(N,N);
		ret+=test(C,N,DIAG);
	}
#endif
#ifdef HAVE_LIBLAPACK
	{
		cout << "TridiagMatrix solve test: ";
		TridiagMatrix T(N);
		ret+=td_test(T,N);
	}
#endif
	{
		ret+=is_symm_test(N);
	}
	return ret;
}

