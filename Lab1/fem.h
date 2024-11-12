#pragma once
#ifndef FEM_H
#define FEM_H
#include "common.h"

class FEM
{
public:
	int nx, ny, nz, n_ke, basisSize = 12;

	double mu = 1.0;
	vector<double> x_coor, y_coor, z_coor;
	vector<Point> centerEdges;
	int numEdges;
	vector<int> ind_funcA;
	FEM(string area, string cond);
	void build();
	vector<double> calc_f_local(int p, int s, int q);
	void create_mesh(string area);
	void resize_matrix();
	Matrix A;
	vector<double> b, q;
	vector<double> b_local;
	vector<vector<double>> A_local, G, M,
		G1 = { {2., 1., -2., -1.},
			   {1., 2., -1., -2.},
			   {-2., -1., 2., 1.},
			   {-1., -2., 1., 2.} },
		G2 = { {2., -2., 1., -1.},
			   {-2., 2., -1., 1.},
			   {1., -1., 2., -2.},
			   {-1., 1., -2., 2.} },
		G3 = { {-2., 2., -1., 1.},
			   {-1., 1., -2., 2.},
			   {2., -2., 1., -1.},
			   {1., -1., 2., -2.} },
		D =  { {4., 2., 2., 1.},
			   {2., 4., 1., 2.},
			   {2., 1., 4., 2.},
			   {1., 2., 2., 4.} };
	void makeA_Local_And_b_Local(int p, int s, int q);
	void addLocalToGlobal(int p, int s, int q);
	void makeG(double hx, double hy, double hz);
	void makeM(double hx, double hy, double hz);
	void addFirstCondition();
	double valueAtEdge(Point edge0, Point edge);
	double getTau(Point value, Point edge0, Point edge);
	Point findElem(Point edge0, Point edge);
	vector<double> getSolution(Point elem);
	Point getPsi(int iFunc, Point point);
	set<int> cond_1;
	void initFirstCondition(string cond);
	void calc_centerEdges();
	Point calc_localCenter(int p, int s, int q, int i);
	void export_q();
};

#endif // !FEM_H

