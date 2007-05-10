#include "meshenum.h"
#include "dmesh.h"

#include <memory>
#include <cstdio>

using namespace std;

int main(int argc, char *argv[]) {
	int WIDTH=6;
	int HEIGHT=5;
	int ret=0;
	std::auto_ptr<DMesh> m(new DMesh(WIDTH,HEIGHT));
	BCSet bcs(BC::createDirichletBC(),
		BC::createNeumannBC(),
		BC::createDirichletBC(),
		BC::createNeumannBC());
	SkipDirichletEnumerator k(*m,bcs);
	try {
		printf("ksize=%-d: ",k.size());
		if (k.size() != (WIDTH*(HEIGHT-2))) {
			printf("FAILED\n");
			ret--;
		} else {
			printf("OK\n");
		}
		for (int j=0; j<(HEIGHT); ++j) {
			for (int i=0; i<(WIDTH); ++i) {
				printf("%4d ",k(i,j));
				if ((j!=0) && (j!=(HEIGHT-1))) {
					if (k(i,j) < 0) {
						printf("\nFAILED\n");
						ret--;
					}
				} else {
					if (k(i,j) > 0) {
						printf("\nFAILED\n");
						ret--;
					}
				}
			}
			printf("\n");
		}
		printf("(i,j)->(k) conversion: OK\n");
		for (int j=1; j<(HEIGHT-1); ++j) {
			for (int i=0; i<(WIDTH); ++i) {
			if ((k.i(k(i,j)) != i)) {
				printf("(k)->(i) conversion at (%d,%d): "
					"FAILED\n",i,j);
				ret--;
			}
			if ((k.j(k(i,j)) != j)) {
				printf("(k)->(j) conversion at (%d,%d): "
					"FAILED\n",i,j);
				ret--;
			}
			}
		}
		printf("(k)->(i,j) conversion: OK\n");
	} catch (EnumeratorException& e) {
		printf("Enumeration exception: i=%d, j=%d, k=%d\n",
			e.i, e.j, e.k);
		printf("FAILED\n");
		ret--;
	}
	try {
		printf("i(k=-4)=%d\n",k.i(-4));
		printf("No exception, FAILED\n");
		ret--;
	} catch (EnumeratorException& e) {
		printf("Enumeration exception: i=%d, j=%d, k=%d\n",
			e.i, e.j, e.k);
		printf("invalid (k)->(i,j) conversion: OK\n");
	}
	return ret;
}

