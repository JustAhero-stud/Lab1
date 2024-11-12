#pragma once
#ifndef COMMON_H
#define COMMON_H

#include <vector>
#include <set>
#include <iostream>
#include <string>  
#include <fstream>
#include <cmath>
#include <vector>


using namespace std;

struct Matrix
{
	int N = 0;
	std::vector<int> IA, JA;
	std::vector<double> DI, AL;

};

struct Point {
	double x, y, z;
	Point(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {};
	Point() : x(0), y(0), z(0) {}
};

extern int getGlobalNumber(int xNum, int yNum, int p, int s, int q, int i);
extern Point getMiddleEdge();
#endif