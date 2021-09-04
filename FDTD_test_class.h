#pragma once

#define _USE_MATH_DEFINES
#include <algorithm>
#include <iostream>
#include <fstream>
#include <clocale>
#include <chrono>
#include <vector>
#include <cmath>
#include <omp.h>

using namespace std;
#define Nx 64 
#define Ny 64
#define Nz 64
template<class ftype>
class MaxwellsSolver {
private:
	int Nt;
	ftype dx, dy, dz, dt, dt_x, dt_y, dt_z;
	ftype T;
	pair<ftype, ftype> ab, cd, fg;
	ftype (*answerEx)(ftype x, ftype y, ftype z, ftype t);
	vector<vector<vector<Component<ftype>>>> cube(Nx + 2, vector<vector<Component<ftype>>>(Ny + 2, vector<Component<ftype>>(Nz + 2)));
	vector<vector<vector<Component<ftype>>>> cube2(Nx + 2, vector<vector<Component<ftype>>>(Ny + 2, vector<Component<ftype>>(Nz + 2)));

public:
	MaxwellsSolver()
	{
		ab.first = 0.0;  ab.second = 2.0 * M_PI;
		cd.first = 0.0;  cd.second = 2.0 * M_PI;
		fg.first = 0.0;  fg.second = 2.0 * M_PI;
		dx = (ab.second - ab.first) / (ftype)Nx; dy = (cd.second - cd.first) / (ftype)Ny;
		dz = (fg.second - fg.first) / (ftype)Nz;
		Nt = T / dt;
	}
};
