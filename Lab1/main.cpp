#include "fem.h"
#include "common.h"
#include "portrait.h"
#include "solver.h"

int main() {

	FEM test("area.txt", "cond_1.txt");
	Portrait t(test.nx, test.ny, test.nz);
	t.build(test.A);
	test.build();
	Solver slae_solver(test.A, test.b, test.q);
	test.export_q();

	//Внутренее ребро
	Point p0(7.0, 0.0, 3.0);
	Point p(7.0, 10.0, 3.0);
	cout << test.valueAtEdge(p0, p);
	//test github 
	int a = 1;
	return a;
}