#include "fem.h"

//#define Fx(x, y, z) (1.0)
//#define Fy(x, y, z) (x * x * x * x - 12.0 * x * x)
//#define Fz(x, y, z) (1.0)
//#define Ax(x, y, z) (1.0)
//#define Ay(x, y, z) (x * x * x * x)
//#define Az(x, y, z) (1.0)

//#define Fx(x, y, z) (1.0)
//#define Fy(x, y, z) (x * x * x - 6.0 * x)
//#define Fz(x, y, z) (1.0)
//#define Ax(x, y, z) (1.0)
//#define Ay(x, y, z) (x * x * x)
//#define Az(x, y, z) (1.0)

//#define Fx(x, y, z) (1.0)
//#define Fy(x, y, z) (x * x - 2.0)
//#define Fz(x, y, z) (1.0)
//#define Ax(x, y, z) (1.0)
//#define Ay(x, y, z) (x * x)
//#define Az(x, y, z) (1.0)

#define Fx(x, y, z) (1.0)
#define Fy(x, y, z) (z)
#define Fz(x, y, z) (1.0)
#define Ax(x, y, z) (1.0)
#define Ay(x, y, z) (z)
#define Az(x, y, z) (1.0)

//#define Fx(x, y, z) (1.0)
//#define Fy(x, y, z) (5.0)
//#define Fz(x, y, z) (10.0)
//#define Ax(x, y, z) (1.0)
//#define Ay(x, y, z) (5.0)
//#define Az(x, y, z) (10.0)


FEM::FEM(string area, string cond) {
	create_mesh(area);
	resize_matrix();
	initFirstCondition(cond);
	calc_centerEdges();
	b.resize(numEdges);
	q.resize(numEdges);
}

void FEM::build()
{
	for (int q = 0; q < nz; q++) {
		for (int s = 0; s < ny; s++) {
			for (int p = 0; p < nx; p++) {
				makeA_Local_And_b_Local(p, s, q);
				addLocalToGlobal(p, s, q);
			}
		}
	}
	addFirstCondition();
}

Point FEM::calc_localCenter(int p, int s, int q, int i) {
	double x = (x_coor[p + 1] + x_coor[p]) / 2.0;
	double y = (y_coor[s + 1] + y_coor[s]) / 2.0;
	double z = (z_coor[q + 1] + z_coor[q]) / 2.0;

	switch (i)
	{
	case 0: return Point(x, y_coor[s], z_coor[q]);
	case 1: return Point(x, y_coor[s + 1], z_coor[q]);
	case 2: return Point(x, y_coor[s], z_coor[q + 1]);
	case 3: return Point(x, y_coor[s + 1], z_coor[q + 1]);
	case 4: return Point(x_coor[p], y, z_coor[q]);
	case 5: return Point(x_coor[p + 1], y, z_coor[q]);
	case 6: return Point(x_coor[p], y, z_coor[q + 1]);
	case 7: return Point(x_coor[p + 1], y, z_coor[q + 1]);
	case 8: return Point(x_coor[p], y_coor[s], z);
	case 9: return Point(x_coor[p + 1], y_coor[s], z);
	case 10: return Point(x_coor[p], y_coor[s + 1], z);
	case 11: return Point(x_coor[p + 1], y_coor[s + 1], z);
	default: return Point();
	}
}

void FEM::calc_centerEdges() {
	centerEdges.resize(numEdges);
	ind_funcA.resize(numEdges);
	for (int q = 0; q < nz; q++) {
		for (int s = 0; s < ny; s++) {
			for (int p = 0; p < nx; p++) {
				for (int i = 0; i < 12; i++)
				{
					int num = getGlobalNumber(nx, ny, p, s, q, i);
					/*if (num == 25) {
						cout <<i << ": " << p << " " << s << " " << q << endl;
					}*/
					centerEdges[num] = calc_localCenter(p, s, q, i);
					if (i >= 0 && i < 4) ind_funcA[num] = 0;
					if (i >= 4 && i < 8) ind_funcA[num] = 1;
					if (i >= 8 && i < 12) ind_funcA[num] = 2;
				}
			}
		}
	}
}

void FEM::initFirstCondition(string cond) {
	ifstream in(cond);
	int n;
	in >> n;

	for (int i = 0; i < n; i++) {
		int num_first;
		in >> num_first;
		switch (num_first) {
		case 0: {//нижн€€
			for (int p = 0; p < nx * (ny + 1) + ny * (nx + 1); p++) 
				cond_1.insert(p);
			break;
		}
		case 1: { //верхн€€
			int start = getGlobalNumber(nx, ny, 0, 0, nz-1, 2);
			for (int p = start; p < numEdges; p++)
				cond_1.insert(p);

			break;
		}
		case 2: { //лева€
			for (int q = 0; q < nz; q++) {
				int start = getGlobalNumber(nx, ny, 0, 0, q, 8);
				for (int p = 0; p < ny + 1; p++)
					cond_1.insert(start + p * (nx + 1));
			}
			for (int p = 0; p < ny; p++) {
				for (int q = 0; q < nz; q++) {
					cond_1.insert(getGlobalNumber(nx, ny, 0, p, q, 4));
				}
				cond_1.insert(getGlobalNumber(nx, ny, 0, p, nz - 1, 6));
			}
			break;
		}
		case 3: { //права€
			for (int q = 0; q < nz; q++) {
				int start = getGlobalNumber(nx, ny, nx - 1, 0, q, 9);
				for (int p = 0; p < ny + 1; p++)
					cond_1.insert(start + p * (nx + 1));
			}
			for (int p = 0; p < ny; p++) {
				for (int q = 0; q < nz; q++) {
					cond_1.insert(getGlobalNumber(nx, ny, nx - 1, p, q, 5));
				}
				cond_1.insert(getGlobalNumber(nx, ny, nx - 1, p, nz - 1, 7));
			}
			break;
		}
		case 4: { //передн€€
			for (int q = 0; q < nz; q++) {
				int start = getGlobalNumber(nx, ny, 0, 0, q, 8);
				for (int p = 0; p < nx + 1; p++)
					cond_1.insert(start + p);
			}
			for (int s = 0; s < nx; s++) {
				for (int q = 0; q < nz; q++) {
					cond_1.insert(getGlobalNumber(nx, ny, s, 0, q, 0));
				}
				cond_1.insert(getGlobalNumber(nx, ny, s, 0, nz - 1, 2));
			}
			break;
		}
		case 5: { //задн€€
			for (int q = 0; q < nz; q++) {
				int start = getGlobalNumber(nx, ny, 0, ny - 1, q, 10);
				for (int p = 0; p < nx + 1; p++)
					cond_1.insert(start + p);
			}
			for (int s = 0; s < nx; s++) {
				for (int q = 0; q < nz; q++) {
					cond_1.insert(getGlobalNumber(nx, ny, s, ny - 1, q, 1));
				}
				cond_1.insert(getGlobalNumber(nx, ny, s, ny - 1, nz - 1, 3));
			}
			break;	
		}
		}
	}
}

void FEM::addFirstCondition() {
	double B = 10e+15;
	for (int cond : cond_1){
		int index = cond;
		if (ind_funcA[index] == 0) b[index] = 
			Ax(centerEdges[index].x, centerEdges[index].y, centerEdges[index].z) * B;
		if (ind_funcA[index] == 1) b[index] = 
			Ay(centerEdges[index].x, centerEdges[index].y, centerEdges[index].z) * B;
		if (ind_funcA[index] == 2) b[index] = 
			Az(centerEdges[index].x, centerEdges[index].y, centerEdges[index].z) * B;
		A.DI[index] = B;
	}
}

vector<double> FEM::calc_f_local(int p, int s, int q) {
	vector<double> f_loc(basisSize);
	double x = (x_coor[p + 1] + x_coor[p]) / 2.0;
	double y = (y_coor[s + 1] + y_coor[s]) / 2.0;
	double z = (z_coor[q + 1] + z_coor[q]) / 2.0;

	f_loc[0] = Fx(x, y_coor[s], z_coor[q]);
	f_loc[1] = Fx(x, y_coor[s+1], z_coor[q]);
	f_loc[2] = Fx(x, y_coor[s], z_coor[q+1]);
	f_loc[3] = Fx(x, y_coor[s+1], z_coor[q+1]);

	f_loc[4] = Fy(x_coor[p], y, z_coor[q]);
	f_loc[5] = Fy(x_coor[p+1], y, z_coor[q]);
	f_loc[6] = Fy(x_coor[p], y, z_coor[q+1]);
	f_loc[7] = Fy(x_coor[p+1], y, z_coor[q+1]);

	f_loc[8] = Fz(x_coor[p], y_coor[s], z);
	f_loc[9] = Fz(x_coor[p+1], y_coor[s], z);
	f_loc[10] = Fz(x_coor[p], y_coor[s+1], z);
	f_loc[11] = Fz(x_coor[p+1], y_coor[s+1], z);

	//f_loc[0] = Fx(x_coor[p], y, z_coor[q]);
	//f_loc[1] = Fx(x_coor[p + 1], y, z_coor[q]);
	//f_loc[2] = Fx(x_coor[p], y, z_coor[q + 1]);
	//f_loc[3] = Fx(x_coor[p + 1], y, z_coor[q + 1]);

	//f_loc[4] = Fy(x, y_coor[s], z_coor[q]);
	//f_loc[5] = Fy(x, y_coor[s + 1], z_coor[q]);
	//f_loc[6] = Fy(x, y_coor[s], z_coor[q + 1]);
	//f_loc[7] = Fy(x, y_coor[s + 1], z_coor[q + 1]);

	//f_loc[8] = Fz(x_coor[p], y_coor[s], z);
	//f_loc[9] = Fz(x_coor[p + 1], y_coor[s], z);
	//f_loc[10] = Fz(x_coor[p], y_coor[s + 1], z);
	//f_loc[11] = Fz(x_coor[p + 1], y_coor[s + 1], z);

	return f_loc;
}

void FEM::makeA_Local_And_b_Local(int p, int s, int q) {
	double hx = x_coor[p + 1] - x_coor[p];
	double hy = y_coor[s + 1] - y_coor[s];
	double hz = z_coor[q + 1] - z_coor[q];
	makeM(hx, hy, hz);
	makeG(hx, hy, hz);
	
	vector<double> f_local = calc_f_local(p, s, q);

	b_local.clear();
	b_local.resize(basisSize);

	for (int i = 0; i < basisSize; i++) {
		for (int j = 0; j < basisSize; j++) {
			A_local[i][j] = G[i][j] + M[i][j];
			b_local[i] += M[i][j] * f_local[j];
		}
	}

}

void FEM::addLocalToGlobal(int p, int s, int q) {
	vector<int> globalNum(basisSize);
	for (int i = 0; i < basisSize; i++)
		globalNum[i] = getGlobalNumber(nx, ny, p, s, q, i);

	for (int i = 0; i < basisSize; i++)
	{
		A.DI[globalNum[i]] += A_local[i][i];
		b[globalNum[i]] += b_local[i];

		for (int j = 0; j < i; j++)
		{
			auto currentNum = globalNum[i];
			auto nextNum = globalNum[j];
			if (currentNum < nextNum)
				std::swap(currentNum, nextNum);

			auto begin = A.JA.begin() + A.IA[currentNum];
			if (A.IA[currentNum + 1] > A.IA[currentNum])
			{
				auto end = A.JA.begin() + A.IA[currentNum + 1] - 1;
				auto iter = lower_bound(begin, end, nextNum);
				auto index = iter - A.JA.begin();
				A.AL[index] += A_local[i][j];
			}
		}
	}
}


void FEM::makeG(double hx, double hy, double hz) {
	double koef_1 = (hx * hy) / (6.0 * hz);
	double koef_2 = (hx * hz) / (6.0 * hy);
	double koef_3 = -1.0*hz / 6.0;
	double koef_4 = hy / 6;
	double koef_5 = (hy * hz) / (6 * hx);
	double koef_6 = -1.0*hx / 6.0;
	double inv_mu = 1.0 / mu;
	//диагональные блоки
	for (int k = 0; k < basisSize; k += 4) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				if (k == 0) G[k + i][k + j] = inv_mu * (koef_1 * G1[i][j] + koef_2 * G2[i][j]);
				if (k == 4) G[k + i][k + j] = inv_mu * (koef_1 * G1[i][j] + koef_5 * G2[i][j]);
				if (k == 8) G[k + i][k + j] = inv_mu * (koef_2 * G1[i][j] + koef_5 * G2[i][j]);
			}
		}
	}
	//остальные блоки
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			G[4 + i][j] = G[i][4 + j] = inv_mu * koef_3 * G2[i][j];
			G[8 + i][4 + j] = G[4 + i][8 + j] = inv_mu * koef_6 * G1[i][j];
			G[i][8 + j] = inv_mu * koef_4 * G3[i][j];
			G[8 + i][j] = inv_mu * koef_4 * G3[j][i];
		}
	}	
}

void FEM::makeM(double hx, double hy, double hz) {
	for (int k = 0; k < basisSize; k += 4) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				M[k + i][k + j] = hx * hy * hz / 36.0 * D[i][j];
			}
		}
	}
}

vector<double> split_area(double min, double max, int n) {
	vector<double> coor;
	coor.resize(n+1);
	double h = (max - min) / n;
	for (int i = 0; i < n; i++) {
		coor[i] = min;
		min += h;
	}
	coor[n] = max;
	return coor;
}

void FEM::create_mesh(string area) {
	ifstream fin(area);
	fin >> nx >> ny >> nz;
	double min, max;
	fin >> min >> max;
	x_coor = split_area(min, max, nx);
	fin >> min >> max;
	y_coor = split_area(min, max, ny);
	fin >> min >> max;
	z_coor = split_area(min, max, nz);
	n_ke = nx * ny * nz;
	numEdges = (nz + 1) * (nx * (ny + 1) + ny * (nx + 1)) + (nx + 1) * (ny + 1) * nz;
}

int getGlobalNumber(int xNum, int yNum, int p, int s, int q, int i) {
	int jump_xy = p + s * (2 * xNum);
	int jump_xyz = 3 * xNum * yNum + 2*(xNum + yNum) + 1;
	int jump = (xNum * (yNum + 1) + yNum * (xNum + 1))* 
		(q + 1) + (xNum + 1) * (yNum + 1) * q + p;
	switch (i)
	{
	case 0: return jump_xy + jump_xyz * q + s;
	case 1: return jump_xy + jump_xyz * q + 2 * xNum + 1 + s;
	case 2: return jump_xy + jump_xyz * (q + 1) + s;
	case 3: return jump_xy + jump_xyz * (q + 1) + 2 * xNum + 1 + s;
	case 4: return jump_xy + jump_xyz * q + xNum + s;
	case 5: return jump_xy + jump_xyz * q + xNum + 1 + s;
	case 6: return jump_xy + jump_xyz * (q + 1) + xNum + s;
	case 7: return jump_xy + jump_xyz * (q + 1) + xNum + 1 + s;
	case 8: return jump + s * (xNum + 1);
	case 9: return jump + 1 + s * (xNum + 1);
	case 10: return jump + xNum + 1 + s * (xNum + 1);
	case 11: return jump + xNum + 2 + s * (xNum + 1);
	default: return 0;
	}
}

Point getMiddleEdge(Point edge0, Point edge) {
	Point middle;
	middle.x = (edge.x + edge0.x) / 2.0;
	middle.y =  (edge.y + edge0.y) / 2.0;
	middle.z =  (edge.z + edge0.z) / 2.0;
	return middle;
}

void FEM::resize_matrix() {
	M.resize(basisSize);
	A_local.resize(basisSize);
	G.resize(basisSize);

	for (int i = 0; i < basisSize; i++) {
		M[i].resize(basisSize);
		A_local[i].resize(basisSize);
		G[i].resize(basisSize);
	}
}

void FEM::export_q()
{
	std::ofstream out_q("q.txt");

	for (int i = 0; i < numEdges; i++)
		out_q << "(" << centerEdges[i].x << ", " << centerEdges[i].y << ", " << centerEdges[i].z << ") \t " << q[i] << std::endl;
	out_q.close();
}



Point FEM::findElem(Point edge0, Point edge) {
	Point Elem;
	Point middle = getMiddleEdge(edge0, edge);
	
	for (int p = 0; p < nx; p++) {
		if (middle.x >= x_coor[p] && middle.x <= x_coor[p + 1])
			Elem.x = p;
	}
	for (int s = 0; s < ny; s++) {
		if (middle.y >= y_coor[s] && middle.y <= y_coor[s + 1])
			Elem.y = s;
	}

	for (int q = 0; q < nz; q++) {
		if (middle.z >= z_coor[q] && middle.z <= z_coor[q + 1])
			Elem.z = q;
	}
	return Elem;
}

vector<double> FEM::getSolution(Point elem) {
	vector<double> solution(12);
	for (int i = 0; i < 12; i++) 
		solution[i] = getGlobalNumber(nx, ny, elem.x, elem.y, elem.z, i);
	return solution;
}

double FEM::getTau(Point value, Point edge0, Point edge) {
	double distance = sqrt((edge.x - edge0.x) * (edge.x - edge0.x) + (edge.y - edge0.y) * (edge.y - edge0.y) + (edge.z - edge0.z) * (edge.z - edge0.z));
	return value.x * (edge.x - edge0.x) / distance + value.y * (edge.y - edge0.y) / distance + value.z * (edge.z - edge0.z) / distance;
}

double FEM::valueAtEdge(Point edge0, Point edge) {
	Point elem = findElem(edge0, edge);
	Point middle = getMiddleEdge(edge0, edge);
	Point ksi_eta_zeta((middle.x - x_coor[elem.x]) / (x_coor[elem.x + 1] - x_coor[elem.x]),
		(middle.y - y_coor[elem.y]) / (y_coor[elem.y + 1] - y_coor[elem.y]),
		(middle.z - z_coor[elem.z]) / (z_coor[elem.z + 1] - z_coor[elem.z]));
	double hx = x_coor[elem.x + 1] - x_coor[elem.x],
	hy = y_coor[elem.y + 1] - y_coor[elem.y],
		hz = z_coor[elem.z + 1] - z_coor[elem.z];
	vector<double> solution = getSolution(elem);
	Point value;
	value.x = q[solution[0]] * ((y_coor[elem.y + 1] - middle.y) * (z_coor[elem.z + 1] - middle.z)) / (hy * hz)
		+ q[solution[1]] * ((middle.y - y_coor[elem.y]) * (z_coor[elem.z + 1] - middle.z)) / (hy * hz)
		+ q[solution[2]] * ((y_coor[elem.y + 1] - middle.y) * (middle.z - z_coor[elem.z])) / (hy * hz)
		+ q[solution[3]] * ((middle.y - y_coor[elem.y]) * (middle.z - z_coor[elem.z])) / (hy * hz);
	value.y = q[solution[4]] * ((x_coor[elem.x + 1] - middle.x) * (z_coor[elem.z + 1] - middle.z)) / (hx * hz)
		+ q[solution[5]] * ((middle.x - x_coor[elem.x]) * (z_coor[elem.z + 1] - middle.z)) / (hx * hz)
		+ q[solution[6]] * ((x_coor[elem.x + 1] - middle.x) * (middle.z - z_coor[elem.z])) / (hx * hz)
		+ q[solution[7]] * ((middle.x - x_coor[elem.x]) * (middle.z - z_coor[elem.z])) / (hx * hz);
	value.z = q[solution[8]] * ((x_coor[elem.x + 1] - middle.x) * (y_coor[elem.y + 1] - middle.y)) / (hx * hy)
		+ q[solution[9]] * ((middle.x - x_coor[elem.x]) * (y_coor[elem.y + 1] - middle.y)) / (hx * hy)
		+ q[solution[10]] * ((x_coor[elem.x + 1] - middle.x) * (middle.y - y_coor[elem.y])) / (hx * hy)
		+ q[solution[11]] * ((middle.x - x_coor[elem.x]) * (middle.y - y_coor[elem.y])) / (hx * hy);

	//for (int i = 0; i < 12; i++)
	//{
	//	Point Psi = getPsi(i, ksi_eta_zeta);
	//	value.x += q[solution[i]] * Psi.x;
	//	value.y += q[solution[i]] * Psi.y;
	//	value.z += q[solution[i]] * Psi.z;
	//}
	return getTau(value, edge0, edge);
}

Point FEM::getPsi(int i, Point point) {
	switch (i)
	{
	case 0: return Point((1 - point.y) * (1 - point.z) / 4.0, 0, 0);
	case 1: return Point((1 + point.y) * (1 - point.z) / 4.0, 0, 0);
	case 2: return Point((1 - point.y) * (1 + point.z) / 4.0, 0, 0);
	case 3: return Point((1 + point.y) * (1 + point.z) / 4.0, 0, 0);
	case 4: return Point(0, (1 - point.x) * (1 - point.z) / 4.0, 0);
	case 5: return Point(0, (1 + point.x) * (1 - point.z) / 4.0, 0);
	case 6: return Point(0, (1 - point.x) * (1 + point.z) / 4.0, 0);
	case 7: return Point(0, (1 + point.x) * (1 + point.z) / 4.0, 0);
	case 8: return Point(0, 0, (1 - point.x) * (1 - point.y) / 4.0);
	case 9: return Point(0, 0, (1 + point.x) * (1 - point.y) / 4.0);
	case 10: return Point(0, 0, (1 - point.x) * (1 + point.y) / 4.0);
	case 11: return Point(0, 0, (1 + point.x) * (1 + point.y) / 4.0);
	}
}



