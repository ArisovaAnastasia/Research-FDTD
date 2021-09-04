#include "half.hpp"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "PML_3dimen_one_array.h"
#include "ionization.h"

using namespace std;
using half_float::half;
using namespace half_float::literal;

//int main(int argc, char** argv) {
//
//	if(argc != 12) {
//		printf("Error: found %d arguments. Need exactly 11", argc - 1);
//		exit(1);
//	}
//
//	string ftype(argv[1]), ftypePML(argv[2]);
//	int Nx = atoi(argv[3]), Ny = atoi(argv[4]), Nz = atoi(argv[5]);
//	double n = atof(argv[6]);
//	int delta_x = atoi(argv[7]), delta_y = atoi(argv[8]), delta_z = atoi(argv[9]);
//	double sigma_x = atof(argv[10]);
//	string type_sum(argv[11]);
//
//	double time;
//	if (ftype == "double") {
//		if (ftypePML == "double") {
//			time = FDTD_3D_PML_one_array_scalability<double, double>(Nx, Ny, Nz, n, delta_x, delta_y, delta_z, sigma_x, type_sum);
//		}
//		if (ftypePML == "float") {
//			time = FDTD_3D_PML_one_array_scalability<double, float>(Nx, Ny, Nz, n, delta_x, delta_y, delta_z, sigma_x, type_sum);
//		}
//		if (ftypePML == "half") {
//			time = FDTD_3D_PML_one_array_scalability<double, half>(Nx, Ny, Nz, n, delta_x, delta_y, delta_z, sigma_x, type_sum);
//		}
//	}
//	if (ftype == "float") {
//		if (ftypePML == "double") {
//			time = FDTD_3D_PML_one_array_scalability<float, double>(Nx, Ny, Nz, n, delta_x, delta_y, delta_z, sigma_x, type_sum);
//		}
//		if (ftypePML == "float") {
//			time = FDTD_3D_PML_one_array_scalability<float, float>(Nx, Ny, Nz, n, delta_x, delta_y, delta_z, sigma_x, type_sum);
//		}
//		if (ftypePML == "half") {
//			time = FDTD_3D_PML_one_array_scalability<float, half>(Nx, Ny, Nz, n, delta_x, delta_y, delta_z, sigma_x, type_sum);
//		}
//	}
//	if (ftype == "half") {
//		if (ftypePML == "double") {
//			time = FDTD_3D_PML_one_array_scalability<half, double>(Nx, Ny, Nz, n, delta_x, delta_y, delta_z, sigma_x, type_sum);
//		}
//		if (ftypePML == "float") {
//			time = FDTD_3D_PML_one_array_scalability<half, float>(Nx, Ny, Nz, n, delta_x, delta_y, delta_z, sigma_x, type_sum);
//		}
//		if (ftypePML == "half") {
//			time = FDTD_3D_PML_one_array_scalability<half, half>(Nx, Ny, Nz, n, delta_x, delta_y, delta_z, sigma_x, type_sum);
//		}
//	}
//	printf("time = %f", time);
//	return 0;
//}
//
//template <class charT, charT sep>
//class punct_facet : public std::numpunct<charT> {
//protected:
//	charT do_decimal_point() const { return sep; }
//};

int main() {

	setlocale(LC_ALL, "english");

	double sigma_x = 7.;
	double sigma_y = 7.;
	double n = 4, time;

	//time = Ionization<double, double>(n, 48, 48, 48, 1000, (half)0.05, (half)0.1, (half)0.1, (half)0.1, 8, 8, 8, sigma_x, sigma_y, "NoKahan", "Solution Ey.csv", ".csv");



	double R1 = 2.67169e-3; // значение R для опт сигмы = 18,86 при T=14pi
	double R2 = 4.52533e-7; // значение R для опт сигмы = 46.5 при T=8pi

	double T = 8. * M_PI;
	int delta_x = 8, delta_y = 8, delta_z = 0;
	int Nx = 128, Ny = 16, Nz = 1;

	FDTD_3D_PML_one_array<half, half>(n, Nx, Ny, Nz, (half)T, 0.005_h,
		delta_x, delta_y, delta_z, 46.5, 46.5, "No Kahan", ".csv", ".csv", ".csv");
	FDTD_3D_PML_one_array<half, half>(n, Nx, Ny, Nz, (half)T, 0.005_h,
		delta_x, delta_y, delta_z, 46.5, 46.5, "Kahan", ".csv", ".csv", ".csv");

	//double sigma = 20.0;

	//ofstream numb_energy("График отражения 1D half-float No Kahan.csv");
	//numb_energy << "version code with formula 2.0" << endl;
	//numb_energy << "delta_x= " << delta_x << "  delta_y= " << delta_y << "   delta_z= " << delta_z << "    " << endl;
	//numb_energy << "Nx= " << Nx << "  Ny= " << Ny << "   Nz= " << Nz << "    " << endl;
	//numb_energy << "n = " << n << endl<<endl;
	//while (sigma <= 100.0)
	//{
	//	numb_energy << sigma << ";"
	//		<< FDTD_3D_PML_one_array<half, float>(n, Nx, Ny, Nz, (half)T, 0.005_h,
	//				delta_x, delta_y, delta_z, sigma, sigma_y, "No Kahan", ".csv", ".csv", ".csv") << endl;
	//	sigma += 0.5;
	//}


	//numb_energy << endl << endl;
	//numb_energy.close();

	//Двумерный график отражения в зависимости от sigma и n

	//double sigma = 79.0, res, n = 1.0, step_sigma = 0.5, step_n = 0.5;
	//

	//ofstream numb_energy("Граф отражения от сигма для double-double Kahan 1D_.csv");
	//numb_energy << endl << ";" << "n" << endl;
	//numb_energy << "sigma" << ";" << ";";
	//while (sigma <= 100.0)
	//{
	//	numb_energy << sigma << ";";
	//	sigma += step_sigma;
	//}
	//numb_energy << endl;
	//while (n <= 4.0)
	//{
	//	numb_energy << ";" << n << ";";
	//	sigma = 20.0;
	//	while (sigma <= 100.0)
	//	{
	//		res = FDTD_three_dimen_with_PML_one_array<double, double>(n, 128, 1, 1, 8.0_h * 3.1415_h, 0.005_h, delta_x, delta_y, delta_z, sigma, sigma_y, "Kahan", ".csv", ".csv", ".csv");
	//		numb_energy << res <<";";
	//		sigma += step_sigma;
	//	}
	//	numb_energy<< endl;
	//	n += step_n;
	//}
	//

	//ofstream numb_energy("Граф отражения от сигма для half-half Kahan 2D.csv");
	//numb_energy << endl << ";" << "sigma" << endl;
	//numb_energy << "n" << ";" << ";";
	//while (n <= 4.0)
	//{
	//	numb_energy << n << ";";
	//	n += step_n;
	//}
	//numb_energy << endl;
	//while (sigma <= 100.0)
	//{
	//	numb_energy << ";" << sigma << ";";
	//	n = 1.0;
	//	while (n <= 4.0)
	//	{
	//		res = FDTD_3D_PML_one_array<half, half>(n, 128, 16, 1, 8.0_h * 3.1415_h, 0.005_h, delta_x, delta_y, delta_z, sigma, sigma_y, "Kahan", ".csv", ".csv", ".csv");
	//		numb_energy << res << ";";
	//		n += step_n;
	//	}
	//	numb_energy << endl;
	//	sigma += step_sigma;
	//}
	////
	//numb_energy << endl << endl;
	//numb_energy.close();

	return 0;
}
