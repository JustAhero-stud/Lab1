#include "solver.h"

Solver::Solver(Matrix& A_glob, std::vector<double>& b_glob, std::vector<double>& result)
{
	this->A = A_glob;
	this->b = b_glob;

	r.resize(A.N);
	z.resize(A.N);
	p.resize(A.N);
	temp1.resize(A.N);
	temp2.resize(A.N);

	LOS(result);
}

void Solver::mult(std::vector<double>& MV, std::vector<double> vec) {
	for (int i = 0; i < A.N; i++) {
		int k0 = A.IA[i];
		int k1 = A.IA[i + 1];
		MV[i] = A.DI[i] * vec[i];
		for (int k = k0; k < k1; k++) {
			int j = A.JA[k];
			MV[i] += vec[j] * A.AL[k];
			MV[j] += vec[i] * A.AL[k];
		}
	}
}

double Solver::skal_mult(std::vector<double> vec1, std::vector<double> vec2) {
	double s = 0;
	for (int i = 0; i < A.N; i++) {
		s += vec1[i] * vec2[i];
	}
	return s;
}

void Solver::LOS(std::vector<double>& result)
{
	// инициализация
	mult(temp1, result);
	for (int i = 0; i < A.N; i++) {
		temp2[i] = b[i] - temp1[i];
	}
	r = temp2;
	z = temp2;
	mult(p, r);

	//iteration
	double nev = skal_mult(r, r);
	int k = 0;
	for (; k < max_it && nev > err; k++) {
		skal1 = skal_mult(p, r);
		skal2 = skal_mult(p, p);

		double alfa = skal1 / skal2;
		for (int i = 0; i < A.N; i++) {
			result[i] += alfa * z[i];
			r[i] -= alfa * p[i];
		}

		mult(temp1, r);
		skal1 = skal_mult(p, temp1);

		double beta = -skal1 / skal2;

		for (int i = 0; i < A.N; i++) {
			z[i] = r[i] + beta * z[i];
			p[i] = temp1[i] + beta * p[i];
		}

		nev = skal_mult(r, r);
	}
	std::cout << "iter:" << k << " nev = " << nev << std::endl;
}