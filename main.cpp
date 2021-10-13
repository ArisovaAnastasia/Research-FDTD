#include "half.hpp"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "PML_3dimen_one_array.h"
#include "ionization.h"

using namespace std;
using half_float::half;
using namespace half_float::literal;

int main() {

	setlocale(LC_ALL, "english");

	double sigma_x = 7.;
	double sigma_y = 7.;
	double n = 4, time;

	//time = Ionization<double, double>(n, 48, 48, 48, 1000, (half)0.05, (half)0.1, (half)0.1, (half)0.1, 8, 8, 8, sigma_x, sigma_y, "NoKahan", "Solution Ey.csv", ".csv");
	

	double R1 = 2.67169e-3; // значение R для опт сигмы = 18,86 при T=14pi
	double R2 = 4.52533e-7; // значение R для опт сигмы = 46.5 при T=8pi

	double T = 8.0 * M_PI;
	int delta_x = 8, delta_y = 8, delta_z = 0;
	int Nx = 128, Ny = 64, Nz = 1;
	cout << T << endl;
	FDTD_3D_PML_one_array<double, double>(n, Nx, Ny, Nz, (double)T, (double)0.01,
		delta_x, delta_y, delta_z, 46.5, 46.5, "No Kahan", "Energy.csv", ".csv", ".csv");

	return 0;
}
