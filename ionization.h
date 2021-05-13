#define _USE_MATH_DEFINES
#include <iostream>
#include <chrono>
#include <vector>
#include <math.h>
#include <cmath>
#include <omp.h>
#include "half.hpp"
#include "formula PML 3D 1_array.h"
#include "formula PML 3D ionization.h"
#include "optional features 1_array.h"
#include "type data.h"

using namespace std;
using half_float::half;
using namespace half_float::literal;


template <class ftype>
void Add_Currents_3_dimen_ionization(data3d<Component<ftype>>& cube,
	int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z, ftype dx, ftype dy, ftype dz, const ftype dt, int it)
{
	ftype x1, x2, y1, y2, z1, z2, t1, t2;
	ftype f1, g1, f2, g2;

	ftype A = (ftype)0.05; //amplitude
	ftype coeff = (ftype)dt / (ftype)(dz);

	t1 = dt * (ftype)it;
	t2 = t1 - (ftype)0.5 * dt;

	ftype x0 = (ftype)(Nx / 2. + delta_x) * dx;
	ftype y0 = (ftype)(Ny / 2. + delta_y) * dy;

	int offset = 10;
	z1 = (ftype)(dz * (ftype)offset + (ftype)0.5 * dz);
	z2 = (ftype)(dz * (ftype)offset);

	ftype omega = (ftype)7.825, waist = (ftype)1.2;
	ftype tp = (ftype)(10. / 3.3333);
	ftype td_const = (ftype)9.;

	ftype lf = (ftype)3.;

	ftype td;


	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++) {

			if ((i >= delta_x + 1) && (i < Nx + delta_x + 1)
				&& (j >= delta_y + 1) && (j < Ny + delta_y + 1))
			{
				y1 = (ftype)(dy * (ftype)j);
				y2 = (ftype)(dy * (ftype)j - (ftype)0.5 * dy);
				x1 = (ftype)(dx * (ftype)i);
				x2 = (ftype)(dx * (ftype)i - (ftype)0.5 * dx);

				td = td_const - lf * (sqrt(1. + ((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0)) / (lf * lf)) - 1.);

				f1 = exp(-(ftype)2. * (ftype)M_LN2*(t1 - td) * (t1 - td) / (tp * tp)) * sin(omega * (t1 - td));
				f2 = exp(-(ftype)2. * (ftype)M_LN2*(t2 - td) * (t2 - td) / (tp * tp)) * sin(omega * (t2 - td));

				g1 = exp(-((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0)) / ((ftype)2. * waist * waist));
				g2 = exp(-((x2 - x0) * (x2 - x0) + (y2 - y0) * (y2 - y0)) / ((ftype)2. * waist * waist));

				// Ez and Bx has the same y coordinate
				cube(i, j, offset + 1).Ey += -A * coeff * f1 * g1;
				cube(i, j, offset + 1).Bx += A * coeff * f2 * g1;
			}
		}
	//cout << cube(Nx / 2 + delta_x, Ny / 2 + delta_y, 10).Ey << "   " << cube(Nx / 2 + delta_x, Ny / 2 + delta_y, 10).Bx << endl;
}

template <class ftype, class ftypePML>
double Ionization(double n, int Nx, int Ny, int Nz, int Nt, ftype dt, ftype dx, ftype dy, ftype dz, int delta_x, int delta_y, int delta_z, double sigma_x, double sigma_y,
	string type_sum, string file_energy, string file_data)
{
	setlocale(LC_ALL, "Russian");

	int offset_x = 4, offset_y = 4, offset_z = 4;

	data3d<Component<ftype>> cube(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	data3d<ComponentSplit<ftypePML>> cube_split(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	data3d<SIGMA<double>> Sigma(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	data3d<COEFF<ftypePML>> Coeff(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	data3d<ftype> Ne(Nx, Ny, Nz, delta_x, delta_y, delta_z);

	data3d<Component<ftype>> compensator(Nx, Ny, Nz, delta_x, delta_y, delta_z), compensator2(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	data3d<ComponentSplit<ftypePML>> compensatorPML(Nx, Ny, Nz, delta_x, delta_y, delta_z), compensatorPML2(Nx, Ny, Nz, delta_x, delta_y, delta_z);

	pair<ftype, ftype> ab(0.0, 4.0 * M_PI), cd(0.0, 16.0 * M_PI), fg(0.0, 16.0 * M_PI);

	ftype dt_x = dt / (dx), dt_y = dt / (dy), dt_z = dt / (dz);
	ftypePML _1dx = (ftypePML)(1.0 / dx), _1dy = (ftypePML)(1.0 / dy), _1dz = (ftypePML)(1.0 / dz);
	vector<double> vec_energy(Nt + 1), vec_energy_acc(Nt + 1);

	Initializing_cube_split_3_dimen_PML<ftypePML>(cube_split, Nx, Ny, Nz, delta_x, delta_y, delta_z);
	Initializing_Sigma_3_dimen_PML<double>(Sigma, Nx, Ny, Nz, delta_x, delta_y, delta_z, n, sigma_x, sigma_y);
	Initializing_Coeff_3_dimen_PML_correct<ftypePML>(Coeff, Sigma, Nx, Ny, Nz, delta_x, delta_y, delta_z, dt);

	 double t1, t2;
	 t1 = omp_get_wtime();

	if (type_sum == "Kahan")
	{
//		Initializing_Coeff_3_dimen_PML_correct_2_0<ftypePML>(Coeff, Sigma, Nx, Ny, Nz, delta_x, delta_y, delta_z, dt);
//
//		cout << type_sum << endl;
//		for (int it = 0; it < Nt; it++)
//		{
//			//if (it % 1000 == 0) cout << it << endl;
//
//			Add_Currents_3_dimen<ftype>(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z, dx, dy, dz, dt, it, ab, cd, fg);
//
//			vec_energy[it] = CalculateEnergyPML_3dimen(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
//			vec_energy_acc[it] = CalculateEnergyPML_3dimen_accurate(cube, compensator, Nx, Ny, Nz, delta_x, delta_y, delta_z);
//
//			//vec_average_dataE[it] = AverageValueElectric(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
//			//vec_average_dataB[it] = AverageValueMagnetic(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
//			//vec_max_dataE[it] = MaxValueElectric(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
//			//vec_max_dataB[it] = MaxValueMagnetic(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
//
//			//vec_dataEy[it] = (double)cube(Nx/2 + delta_x, Ny/2 + delta_y, Nz/2 + delta_z+1).Ey;
//			//vec_dataBz[it] = (double)cube(Nx / 2 + delta_x, Ny / 2 + delta_y, Nz / 2 + delta_z+1).Bz;
//
//			// cout << cube(Nx, Ny, Nz).Ey << "  " << (double)cube(Nx + 2 * delta - 5, Ny + 2 * delta - 5, Nz + 2 * delta - 5).Bz << endl;
//
//					   
//#pragma omp parallel for collapse(3)
//			for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
//				for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
//					for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
//					{
//						if ((i >= delta_x + 1) && (i < Nx + delta_x + 1) && (j >= delta_y + 1) && (j < Ny + delta_y + 1) && (k >= delta_z + 1) && (k < Nz + delta_z + 1))
//						{
//							Update_electric_field_three_dimen_Kahan<ftype>(cube, compensator, compensator2, dt_x, dt_y, dt_z, i, j, k);
//						}
//						else {
//							//Update_electric_field_three_dimen_PML<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
//							Update_electric_field_three_dimen_PML_Kahan_2_0<ftype, ftypePML>(cube, cube_split, Coeff, compensatorPML, compensatorPML2, compensator, compensator2, _1dx, _1dy, _1dz, i, j, k);
//
//						}
//					}
//			for (int i = 1; i < Nx + 2 * delta_x + 1; i++) 
//				for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
//					for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
//					{
//						compensator(i, j, k).Ex = compensator2(i, j, k).Ex;
//						compensator(i, j, k).Ey = compensator2(i, j, k).Ey;
//						compensator(i, j, k).Ez = compensator2(i, j, k).Ez;
//
//						compensatorPML(i, j, k).Exy = compensatorPML2(i, j, k).Exy;
//						compensatorPML(i, j, k).Exz = compensatorPML2(i, j, k).Exz;
//						compensatorPML(i, j, k).Eyx = compensatorPML2(i, j, k).Eyx;
//						compensatorPML(i, j, k).Eyz = compensatorPML2(i, j, k).Eyz;
//						compensatorPML(i, j, k).Ezx = compensatorPML2(i, j, k).Ezx;
//						compensatorPML(i, j, k).Ezy = compensatorPML2(i, j, k).Ezy;
//					}
//
//#pragma omp parallel for collapse(3)
//			for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
//				for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
//					for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
//					{
//						if ((i >= delta_x + 1) && (i < Nx + delta_x + 1) && (j >= delta_y + 1) && (j < Ny + delta_y + 1) && (k >= delta_z + 1) && (k < Nz + delta_z + 1))
//						{
//							Update_magnetic_field_three_dimen_Kahan<ftype>(cube, compensator, compensator2, dt_x, dt_y, dt_z, i, j, k);
//						}
//						else {
//							//Update_magnetic_field_three_dimen_PML<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
//							Update_magnetic_field_three_dimen_PML_Kahan_2_0<ftype, ftypePML>(cube, cube_split, Coeff, compensatorPML, compensatorPML2, compensator, compensator2, _1dx, _1dy, _1dz, i, j, k);
//
//						}
//					}
//			for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
//				for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
//					for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
//					{
//						compensator(i, j, k).Bx = compensator2(i, j, k).Bx;
//						compensator(i, j, k).By = compensator2(i, j, k).By;
//						compensator(i, j, k).Bz = compensator2(i, j, k).Bz;
//
//						compensatorPML(i, j, k).Bxy = compensatorPML2(i, j, k).Bxy;
//						compensatorPML(i, j, k).Bxz = compensatorPML2(i, j, k).Bxz;
//						compensatorPML(i, j, k).Byx = compensatorPML2(i, j, k).Byx;
//						compensatorPML(i, j, k).Byz = compensatorPML2(i, j, k).Byz;
//						compensatorPML(i, j, k).Bzx = compensatorPML2(i, j, k).Bzx;
//						compensatorPML(i, j, k).Bzy = compensatorPML2(i, j, k).Bzy;
//					}
//		}
	}
	else
	{
		Initializing_Coeff_3_dimen_PML_correct<ftypePML>(Coeff, Sigma, Nx, Ny, Nz, delta_x, delta_y, delta_z, dt);
		cout << type_sum << endl;

		for (int it = 0; it < Nt; it++)
		{
			if ((it + 1) % 100 == 0) {
				FILE* out;
				char tmp[1024];
				sprintf(tmp, "iteration_%d half.txt", it + 1);
				out = fopen(tmp, "w");
				if (NULL == out)
				{
					exit(-1);
				}
				fprintf(out, "%d %d\n", Nx + 2 * delta_x + 2, Ny + 2 * delta_y + 2);
				for (int i = 0; i < Nx + 2 * delta_x + 2; i++)
					for (int j = 0; j < Ny + 2 * delta_y + 2; j++)
					{
						fprintf(out, "%d %d %lf %lf %lf\n", i, j,
							(double)cube(i, 32, j).Ey, (double)cube(i, 32, j).Bx, (double)Ne(i, 32, j));
					}
				fclose(out);
			}
			if (it % 100 == 0) { cout << it << endl; 
			}

			//if (it == 400) {
			//	Graph_E_ionization(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z, dx, dy, dz, "E iteration400_.csv");
			//	Graph_Ne_ionization(Ne, Nx, Ny, Nz, delta_x, delta_y, delta_z, dx, dy, dz, "N iteration400_.csv");
			//}
			//Add_Currents_3_dimen<ftype>(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z, dx, dy, dz, dt, it);

			 Add_Currents_3_dimen_ionization<ftype>(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z, dx, dy, dz, dt, it);
			 vec_energy[it] = cube(Nx/2 + delta_x, Ny/2 + delta_y, Nz/2 + delta_z).Ey;

#pragma omp parallel for collapse(3)
			 for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
				 for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
					 for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
					 {
						 if ((i >= delta_x + 1) && (i < Nx + delta_x + 1) && (j >= delta_y + 1) && (j < Ny + delta_y + 1) && (k >= delta_z + 1) && (k < Nz + delta_z + 1))
						 {
							 Update_magnetic_field_three_dimen<ftype>(cube, dt_x, dt_y, dt_z, i, j, k);
						 }
						 else {
							 Update_magnetic_field_three_dimen_PML<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
						 }
					 }

#pragma omp parallel for collapse(3)
			for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
				for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
					for (int k = 1; k < Nz + 2 * delta_z + 1; k++) {

						if ((i >= delta_x + 1) && (i < Nx + delta_x + 1)
							&& (j >= delta_y + 1) && (j < Ny + delta_y + 1)
							&& (k >= delta_z + 1) && (k < Nz + delta_z + 1))
						{
							Update_electric_field_three_dimen_ionization<ftype>(cube, Ne, dt, dt_x, dt_y, dt_z, i, j, k);

							if ((i >= delta_x  + offset_x + 1) && (i < Nx + delta_x - offset_x + 1)
								&& (j >= delta_y + offset_y + 1) && (j < Ny + delta_y - offset_y + 1)
								&& (k >= delta_z + offset_z + 1) && (k < Nz + delta_z - offset_z + 1))
							{
								Update_electric_currents_three_dimen_ionization(cube, Ne, dt, dt_x, dt_y, dt_z, i, j, k);
							}

						}
						else {
							Update_electric_field_three_dimen_PML<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
						}
					}
			cout << "electric = " << MaxValueElectricModule(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z) << endl;
			//cout << "  max value = " << MaxValueCurrents(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
			
			//cout << "  min value = " << MinValueCurrents(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
			//cout << "  max value Ne = " << MaxValueNe(Ne, Nx, Ny, Nz, delta_x, delta_y, delta_z) << endl;

			}
	}

	t2 = omp_get_wtime();
	vec_energy[Nt] = cube(Nx / 2 + delta_x, Ny / 2 + delta_y, Nz / 2 + delta_z).Ey;

	ofstream numb_energy(file_energy);

	numb_energy << "Nx= " << Nx << ";  Ny= " << Ny << ";   Nz= " << Nz << "    " << endl;
	numb_energy << "delta_x= " << delta_x << ";  delta_y= " << delta_y << ";   delta_z= " << delta_z << "    " << endl;
	numb_energy << sigma_x << ";  " << sigma_y << ";   " << sigma_y << ";    " << endl << endl;


	for (int s = 0; s < Nt + 1; s++)
	{
		numb_energy << dt * (double)(s) << ";" << vec_energy[s] <<endl;
	}
	numb_energy << endl << endl;
	numb_energy.close();


	vec_energy[Nt] = moduleE(cube, Nx / 2 + delta_x, Ny / 2 + delta_y, Nz / 2 + delta_z);

	vec_energy_acc[Nt] = CalculateEnergyPML_3dimen_accurate(cube, compensator, Nx, Ny, Nz, delta_x, delta_y, delta_z);

	cout << "Reflection coefficient = "<< endl << vec_energy[Nt] / vec_energy[2512] << endl << "Accurate reflection coefficient = " << endl << vec_energy_acc[Nt] / vec_energy_acc[2512] << endl << endl;
	cout << "time = " << t2 - t1 << endl << endl;

	return t2-t1;
}

template <class ftype>
void Graph_E_ionization(data3d<Component<ftype>>& cube, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z, ftype dx, ftype dy, ftype dz, string file_name)
{
	ofstream numbE(file_name);

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
	{
		for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
		{
			numbE << moduleE<ftype>(cube, i, Ny / 2 + delta_y, k)<< ";";
		}
		numbE << endl;
	}

	numbE.close();
}

template <class ftype>
void Graph_Ne_ionization(data3d<ftype>& Ne, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z, ftype dx, ftype dy, ftype dz, string file_name)
{
	ofstream numbN(file_name);

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
	{
		for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
		{
			numbN << Ne(i, Ny / 2 + delta_y, k) << ";";
		}
		numbN << endl;
	}

	numbN.close();
}
