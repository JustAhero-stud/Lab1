#pragma once
#ifndef SOLVER_H
#define SOLVER_H
#include "common.h"

class Solver {
public:
	int max_it = 20000;
	double err = 1e-20;
	double skal1, skal2;
	std::vector<double> r, p, z, temp1, temp2;

	Matrix A;
	std::vector<double> b;

	Solver() {};
	Solver(Matrix& A_glob, std::vector<double>& b_glob, std::vector<double>& result);

	void mult(std::vector<double>& MV, std::vector<double> vec);
	double skal_mult(std::vector<double> vec1, std::vector<double> vec2);
	void LOS(std::vector<double>& result);
};
#endif 