#define _USE_MATH_DEFINES
#include <algorithm>
#include <fstream>
#include <clocale>
#include <vector>
#include <cmath>
#include "half.hpp"

using namespace std;

template <class ftypePML>
void Graph_for_Sigma_three_dimen(data3d<SIGMA<double>>& Sigma, int Nx, int Ny, int Nz, int delta, ftypePML dx, string file_sigma)
{
	//����� ������� ���� � ����
	ofstream numb_sigma(file_sigma);

	numb_sigma << ";" << "value sigma_x" << endl;

	for (int s = 1; s < Nx + 2 * delta + 1; s++)
	{
		numb_sigma << dx * (ftypePML)(s - 1) << ";" << Sigma(s, delta / 2, delta / 2).sigmaH_x << endl;
		numb_sigma << dx * ((ftypePML)(s - 1) + 0.5) << ";" << Sigma(s, delta / 2, delta / 2).sigmaE_x << endl;
	}
	numb_sigma << endl << endl;

	numb_sigma.close();
}

template <class ftypePML>
void Graph_for_Coeff_three_dimen(data3d<COEFF<ftypePML>>& Coeff, int Nx, int Ny, int Nz, int delta, ftypePML dx, string file_coeff)
{
	//����� ������� ���� � ����
	ofstream numb_sigma(file_coeff);

	numb_sigma << ";" << "Coeff1" << ";;;" << "Coeff2" << ";;" << "Eyx" << endl;

	for (int s = 1; s < Nx + 2 * delta + 1; s++)
	{
		numb_sigma << dx * (ftypePML)(s - 1) << ";" << std::setprecision(16) << Coeff(s, delta / 2, delta / 2).Bzx1 << ";;";
		numb_sigma << dx * (ftypePML)(s - 1) << ";" << std::setprecision(16) << Coeff(s, delta / 2, delta / 2).Bzx2 << endl;
	}
	numb_sigma << endl << endl;

	numb_sigma.close();
}


template <class ftype>
double CalculateEnergyPML_3dimen(data3d<Component<ftype>>& cube, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	double energy = 0.0;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
			{
				energy += (double)cube(i, j, k).Ex * (double)cube(i, j, k).Ex
					+ (double)cube(i, j, k).Ey * (double)cube(i, j, k).Ey
					+ (double)cube(i, j, k).Ez * (double)cube(i, j, k).Ez;
				energy += (double)cube(i, j, k).Bx * (double)cube(i, j, k).Bx
					+ (double)cube(i, j, k).By * (double)cube(i, j, k).By
					+ (double)cube(i, j, k).Bz * (double)cube(i, j, k).Bz;
			}
	return energy;
}

template <class ftype>
double AverageValueMagnetic(data3d<Component<ftype>>& cube, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	double sum = 0.0;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
			{
				sum += (double)cube(i, j, k).Bx + (double)cube(i, j, k).By + (double)cube(i, j, k).Bz;
			}
	return sum / (double)(3 * (Nx + 2 * delta_x) * (Ny + 2 * delta_y) * (Nz + 2 * delta_z));
}

template <class ftype>
double AverageValueElectric(data3d<Component<ftype>>& cube, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	double sum = 0.0;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
			{
				sum += (double)cube(i, j, k).Ex + (double)cube(i, j, k).Ey + (double)cube(i, j, k).Ez;
			}
	return sum / (double)(3 * (Nx + 2 * delta_x) * (Ny + 2 * delta_y) * (Nz + 2 * delta_z));
}

template <class ftype>
double MaxValueElectric(data3d<Component<ftype>>& cube, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	double sum = DBL_MIN;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
			{
				sum = max(sum, abs((double)cube(i, j, k).Ex));
				sum = max(sum, abs((double)cube(i, j, k).Ey));
				sum = max(sum, abs((double)cube(i, j, k).Ez));
			}
	return sum;
}

template <class ftype>
double MaxValueMagnetic(data3d<Component<ftype>>& cube, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	double sum = DBL_MIN;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
			{
				sum = max(sum, abs((double)cube(i, j, k).Bx));
				sum = max(sum, abs((double)cube(i, j, k).By));
				sum = max(sum, abs((double)cube(i, j, k).Bz));
			}
	return sum;
}

template <class ftype>
double CalculateEnergyPML_3dimen_accurate(data3d<Component<ftype>>& cube, data3d<Component<ftype>>& compensator,
	int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	double energy = 0.0;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
			{
				energy += ((double)cube(i, j, k).Ex - (double)compensator(i, j, k).Ex) * ((double)cube(i, j, k).Ex - (double)compensator(i, j, k).Ex)
					+ ((double)cube(i, j, k).Ey - (double)compensator(i, j, k).Ey) * ((double)cube(i, j, k).Ey - (double)compensator(i, j, k).Ey)
					+ ((double)cube(i, j, k).Ez - (double)compensator(i, j, k).Ez) * ((double)cube(i, j, k).Ez - (double)compensator(i, j, k).Ez);
				energy += ((double)cube(i, j, k).Bx - (double)compensator(i, j, k).Bx) * ((double)cube(i, j, k).Bx - (double)compensator(i, j, k).Bx)
					+ ((double)cube(i, j, k).By - (double)compensator(i, j, k).By) * ((double)cube(i, j, k).By - (double)compensator(i, j, k).By)
					+ ((double)cube(i, j, k).Bz - (double)compensator(i, j, k).Bz) * ((double)cube(i, j, k).Bz - (double)compensator(i, j, k).Bz);
			}
	return energy;
}

template <class ftype>
void Graph_Solution_in_two_planes_3dimen(data3d<Component<ftype>>& cube, int Nx, int Ny, int Nz, int delta, ftype dx, ftype dy, ftype dz, Direction direction, string file_name)
{
	ofstream numbEx("graph_solution_E_x_float-float-sigma 49.6.csv"), numbBx("graph_solution_B_x_float-float-sigma 49.6.csv");
	ofstream numbEy("graph_solution_E_y_float-float-sigma 49.6.csv"), numbBy("graph_solution_B_y_float-float-sigma 49.6.csv");
	ofstream numbEz("graph_solution_E_z_float-float-sigma 49.6.csv"), numbBz("graph_solution_B_z_float-float-sigma 49.6.csv");

	numbEx << direction << ";" << "y" << endl;
	numbEx << "x" << ";" << ";";
	numbBx << direction << ";" << "y" << endl;
	numbBx << "x" << ";" << ";";
	numbEy << direction << ";" << "y" << endl;
	numbEy << "x" << ";" << ";";
	numbBy << direction << ";" << "y" << endl;
	numbBy << "x" << ";" << ";";
	numbEz << direction << ";" << "y" << endl;
	numbEz << "x" << ";" << ";";
	numbBz << direction << ";" << "y" << endl;
	numbBz << "x" << ";" << ";";
	for (int i = 1; i < Nx + 2 * delta + 1; i++)
	{
		numbEx << dx * (ftype)i << ";";
		numbBx << dx * (ftype)i << ";";
		numbEy << dx * (ftype)i << ";";
		numbBy << dx * (ftype)i << ";";
		numbEz << dx * (ftype)i << ";";
		numbBz << dx * (ftype)i << ";";
	}
	numbEx << endl; 		numbBx << endl;	numbEy << endl; 		numbBy << endl;	numbEz << endl; 		numbBz << endl;

	for (int i = 1; i < Ny + 2 * delta + 1; i++)
	{
		ftype y = dy * (ftype)i;
		numbEx << ";" << y << ";"; 			numbBx << ";" << y << ";";			numbEy << ";" << y << ";"; 			numbBy << ";" << y << ";";
		numbEz << ";" << y << ";"; 			numbBz << ";" << y << ";";

		for (int j = 1; j < Nx + 2 * delta + 1; j++)
		{
			numbEx << cube(j, i, Nz / 2 + delta).Ex << ";";
			numbBx << cube(j, i, Nz / 2 + delta).Bx << ";";
			numbEy << cube(j, i, Nz / 2 + delta).Ey << ";";
			numbBy << cube(j, i, Nz / 2 + delta).By << ";";
			numbEz << cube(j, i, Nz / 2 + delta).Ez << ";";
			numbBz << cube(j, i, Nz / 2 + delta).Bz << ";";
		}
		numbEx << endl; 			numbBx << endl;			numbEy << endl; 			numbBy << endl;
		numbEz << endl; 			numbBz << endl;

	}

	numbEx.close();
	numbBx.close();
	numbEy.close();
	numbBy.close();
	numbEz.close();
	numbBz.close();
}

template <class ftype>
void Graph_Solution_in_one_planes_3dimen(data3d<Component<ftype>>& cube, int Nx, int Ny, int Nz, int delta, ftype dx, ftype dy, ftype dz, Direction direction)
{
	ofstream numbEx("graph_solution_E_x_one_dimen_hh-sigma 49.6.csv"), numbBx("graph_solution_B_x_one_dimen_hh-sigma 49.6.csv");
	ofstream numbEy("graph_solution_E_y_one_dimen_hh-sigma 49.6.csv"), numbBy("graph_solution_B_y_one_dimen_hh-sigma 49.6.csv");
	ofstream numbEz("graph_solution_E_z_one_dimen_hh-sigma 49.6.csv"), numbBz("graph_solution_B_z_one_dimen_hh-sigma 49.6.csv");

	for (int i = 1; i < Nx + 2 * delta + 1; i++)
	{
		numbEx << dx * (double)(i) << ";" << cube(i, Ny / 2 + delta, Nz / 2 + delta).Ex << endl;
		numbEy << dx * (double)(i) << ";" << cube(i, Ny / 2 + delta, Nz / 2 + delta).Ey << endl;
		numbEz << dx * (double)(i) << ";" << cube(i, Ny / 2 + delta, Nz / 2 + delta).Ez << endl;

		numbBx << dx * (double)(i) << ";" << cube(i, Ny / 2 + delta, Nz / 2 + delta).Bx << endl;
		numbBy << dx * (double)(i) << ";" << cube(i, Ny / 2 + delta, Nz / 2 + delta).By << endl;
		numbBz << dx * (double)(i) << ";" << cube(i, Ny / 2 + delta, Nz / 2 + delta).Bz << endl;

	}

	numbEx.close();
	numbBx.close();
	numbEy.close();
	numbBy.close();
	numbEz.close();
	numbBz.close();
}

