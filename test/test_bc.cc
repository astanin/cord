#include "../boundary.h"
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main(int argc, char *argv[]) {
	string names[4] = { "unknown", "Dirichlet", "Neumann", "Robin" };
	BoundaryCondition dirichlet(1,0,0);
	BoundaryCondition neumann(0,1,0);
	int type;
	cout << "bc type: " << names[(type=dirichlet.get_type())] << "\n";
	if (type == BoundaryCondition::DIRICHLET_BC) {
		cout << "OK\n";
	} else {
		cout << "Failure\n";
		return -1;
	}
	cout << "bc type: " << names[(type=neumann.get_type())] << "\n";
	if (type == BoundaryCondition::NEUMANN_BC) {
		cout << "OK\n";
	} else {
		cout << "Failure\n";
		return -1;
	}
	BoundaryCondition d2=BoundaryCondition::createDirichletBC();
	cout << "bc type: " << names[(type=dirichlet.get_type())] << "\n";
	if (type == BoundaryCondition::DIRICHLET_BC) {
		cout << "OK\n";
	} else {
		cout << "Failure\n";
		return -1;
	}
	BoundaryCondition n2=BoundaryCondition::createNeumannBC();
	cout << "bc type: " << names[(type=neumann.get_type())] << "\n";
	if (type == BoundaryCondition::NEUMANN_BC) {
		cout << "OK\n";
	} else {
		cout << "Failure\n";
		return -1;
	}
}

