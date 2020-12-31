#include "Maxwells_template.h"
#include "Maxwells_two_dimensional_space.h"
#include "FDTD_YeeGrid.h"
#include "FDTD_YeeGrid_2dimen.h"
#include "half.hpp"
#include <iostream>
#include "PML_2dimen.h"
#include "PML_3dimen.h"
using namespace std;
using half_float::half;
using namespace half_float::literal;

int main() {

	int delta = 8;
	double R = 1.0e-2;
	double T = 6. * 4.0 * M_PI / 6. + 4. * M_PI;

	 FDTD_three_dimen_with_PML<double>(128, 16, 16, T, 0.005_h, delta, R, "No Kahan", "Graph for Yee grid_1.csv");

	Check_Curant<float>(6.28 / 64.0, 6.28 / 64.0, 6.28 / 64.0, 0.005_h);

	return 0;
}
