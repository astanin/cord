#include <ctime>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <blitz/array.h>

using std::time;
using std::cout;
using std::vector;
using std::string;
using std::map;
typedef blitz::Array<double, 1> array1d;

int main(int argc, char*argv[]) {
	int N=10000;
	int T=0;
	// adjust to CPU
	int tstart=time(0);
	do {
		T++;
	} while (time(0) <= (tstart+1));
	T*=100;
	cout << "iterations to do: << " << T << "\n";
	double x;
	double y;
	// direct access
	{
		array1d A(N);
		int t0=time(0);
		for (int i=0; i<T; ++i) {
			x=A(i%N);
			y+=x;
		}
		int t1=time(0);
		cout << "direct array access time: "<<t1-t0
				<<", op/s="<<T*1.0/(t1-t0)<<"\n";
	}
	// array access
	{
		int NVAR=10;
		std::vector<array1d> vA;
		for (int i=0; i<NVAR; ++i) {
			array1d *a=new array1d(N);
			vA.push_back(*a);
		}
		int t0=time(0);
		for (int i=0; i<T; ++i) {
			x=vA[i%NVAR](i%N);
			y+=x;
		}
		int t1=time(0);
		cout << "vector<array> access time: "<<t1-t0
				<<", op/s="<<T*1.0/(t1-t0)<<"\n";
	}
	// map<int,> access
	{
		int NVAR=5;
		vector<int> keys;
		keys.push_back(1);
		keys.push_back(3);
		keys.push_back(5);
		keys.push_back(7);
		keys.push_back(1001);
		map<int,array1d> mA;
		for (int i=0; i<NVAR; ++i) {
			int key=keys[i%NVAR];
			array1d *a=new array1d(N);
			mA.insert(std::make_pair(key,*a));
		}
		int t0=time(0);
		for (int i=0; i<T; ++i) {
			x=mA[7](i%N);
			y+=x;
		}
		int t1=time(0);
		cout << "map<int,array> access time: "<<t1-t0
				<<", op/s="<<T*1.0/(t1-t0)<<"\n";
	}
	// map<string,> access
	{
		int NVAR=5;
		vector<string> keys;
		keys.push_back("aaa");
		keys.push_back("bbb");
		keys.push_back("ccc");
		keys.push_back("ddd");
		keys.push_back("fff");
		map<string,array1d> mA;
		for (int i=0; i<NVAR; ++i) {
			string key=keys[i%NVAR];
			array1d *a=new array1d(N);
			mA.insert(make_pair(key,*a));
		}
		int t0=time(0);
		for (int i=0; i<T; ++i) {
			x=mA["ddd"](i%N);
			y+=x;
		}
		int t1=time(0);
		cout << "map<string,array> access time: "<<t1-t0
				<<", op/s="<<T*1.0/(t1-t0)<<"\n";
	}
}

