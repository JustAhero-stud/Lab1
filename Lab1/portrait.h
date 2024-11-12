#pragma once
#ifndef PORTRAIT_H
#define PORTRAIT_H
#include "common.h"

class Portrait
{
public:
	int nx, ny, nz, numEdges;
	std::vector<std::set<int>> connections;
	Matrix A;

	Portrait(int nx, int ny, int nz);

	void build_connections();
	void build(Matrix& A);
};
#endif 
