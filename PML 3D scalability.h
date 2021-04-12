#define _USE_MATH_DEFINES

#include <iostream>
#include <chrono>
#include <vector>
#include <math.h>
#include <cmath>
#include <omp.h>
#include "half.hpp"
#include "formula PML 3D 1_array.h"
#include "type data.h"

using namespace std;
using half_float::half;
using namespace half_float::literal;



template <class ftypePML>
void Initializing_cube_split_3_dimen_PML(data3d<ComponentSplit<ftypePML>>& cube_split, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	for (int i = 0; i < Nx + 2 * delta_x + 2; i++)
		for (int j = 0; j < Ny + 2 * delta_y + 2; j++)
			for (int k = 0; k < Nz + 2 * delta_z + 2; k++) {
				cube_split(i, j, k).Exy = (ftypePML)0.0;
				cube_split(i, j, k).Exz = (ftypePML)0.0;
				cube_split(i, j, k).Eyx = (ftypePML)0.0;
				cube_split(i, j, k).Eyz = (ftypePML)0.0;
				cube_split(i, j, k).Ezx = (ftypePML)0.0;
				cube_split(i, j, k).Ezy = (ftypePML)0.0;

				cube_split(i, j, k).Bxy = (ftypePML)0.0;
				cube_split(i, j, k).Bxz = (ftypePML)0.0;
				cube_split(i, j, k).Byx = (ftypePML)0.0;
				cube_split(i, j, k).Byz = (ftypePML)0.0;
				cube_split(i, j, k).Bzx = (ftypePML)0.0;
				cube_split(i, j, k).Bzy = (ftypePML)0.0;
			}
}

template <class ftypePML>
void Initializing_Sigma_3_dimen_PML(data3d<SIGMA<ftypePML>>& Sigma, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z, double n,
	ftypePML sigma_x, ftypePML sigma_y)
{
	ftypePML var_Sigma_max_x = sigma_x;
	ftypePML var_Sigma_max_y = sigma_y;
	ftypePML var_Sigma_max_z = sigma_y;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++) {
				Sigma(i, j, k).sigmaE_x = var_Sigma_max_x * powf(distanceE<ftypePML>(Nx, delta_x, i), n);
				Sigma(i, j, k).sigmaH_x = var_Sigma_max_x * powf(distanceH<ftypePML>(Nx, delta_x, i), n);

				Sigma(i, j, k).sigmaE_y = var_Sigma_max_y * powf(distanceE<ftypePML>(Ny, delta_y, j), n);
				Sigma(i, j, k).sigmaH_y = var_Sigma_max_y * powf(distanceH<ftypePML>(Ny, delta_y, j), n);

				Sigma(i, j, k).sigmaE_z = var_Sigma_max_z * powf(distanceE<ftypePML>(Nz, delta_z, k), n);
				Sigma(i, j, k).sigmaH_z = var_Sigma_max_z * powf(distanceH<ftypePML>(Nz, delta_z, k), n);
			}
}

template <class ftypePML>
void Initializing_Coeff_3_dimen_PML_correct_2_0(data3d<COEFF<ftypePML>>& Coeff, data3d<SIGMA<double>>& Sigma, int Nx, int Ny, int Nz,
	int delta_x, int delta_y, int delta_z, ftypePML dt)
{
	double Exy1, Exz1, Ezx1, Ezy1, Eyx1, Eyz1;
	double Bxy1, Bxz1, Bzx1, Bzy1, Byx1, Byz1;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++) {
				if ((i >= delta_x + 1) && (i < Nx + delta_x + 1) && (j >= delta_y + 1) && (j < Ny + delta_y + 1) && (k >= delta_z + 1) && (k < Nz + delta_z + 1))
				{
				}
				else {

					Exy1 = exp(-dt * Sigma(i, j, k).sigmaE_y);
					Exz1 = exp(-dt * Sigma(i, j, k).sigmaE_z);
					Eyx1 = exp(-dt * Sigma(i, j, k).sigmaE_x);
					Eyz1 = exp(-dt * Sigma(i, j, k).sigmaE_z);
					Ezx1 = exp(-dt * Sigma(i, j, k).sigmaE_x);
					Ezy1 = exp(-dt * Sigma(i, j, k).sigmaE_y);

					Bxy1 = exp(-dt * Sigma(i, j, k).sigmaH_y);
					Bxz1 = exp(-dt * Sigma(i, j, k).sigmaH_z);
					Byx1 = exp(-dt * Sigma(i, j, k).sigmaH_x);
					Byz1 = exp(-dt * Sigma(i, j, k).sigmaH_z);
					Bzx1 = exp(-dt * Sigma(i, j, k).sigmaH_x);
					Bzy1 = exp(-dt * Sigma(i, j, k).sigmaH_y);

					Coeff(i, j, k).Exy1 = (ftypePML)Exy1;
					Coeff(i, j, k).Exz1 = (ftypePML)Exz1;
					Coeff(i, j, k).Eyx1 = (ftypePML)Eyx1;
					Coeff(i, j, k).Eyz1 = (ftypePML)Eyz1;
					Coeff(i, j, k).Ezx1 = (ftypePML)Ezx1;
					Coeff(i, j, k).Ezy1 = (ftypePML)Ezy1;

					Coeff(i, j, k).Bxy1 = (ftypePML)Bxy1;
					Coeff(i, j, k).Bxz1 = (ftypePML)Bxz1;
					Coeff(i, j, k).Byx1 = (ftypePML)Byx1;
					Coeff(i, j, k).Byz1 = (ftypePML)Byz1;
					Coeff(i, j, k).Bzx1 = (ftypePML)Bzx1;
					Coeff(i, j, k).Bzy1 = (ftypePML)Bzy1;

					if (Sigma(i, j, k).sigmaE_x != (ftypePML)0.0) {
						Coeff(i, j, k).Eyx2 = (1.0 / Sigma(i, j, k).sigmaE_x - Eyx1 / Sigma(i, j, k).sigmaE_x) / Eyx1;
						Coeff(i, j, k).Ezx2 = (1.0 / Sigma(i, j, k).sigmaE_x - Ezx1 / Sigma(i, j, k).sigmaE_x) / Ezx1;

					}
					else {
						Coeff(i, j, k).Eyx2 = dt;
						Coeff(i, j, k).Ezx2 = dt;
					}
					if (Sigma(i, j, k).sigmaE_y != (ftypePML)0.0) {
						Coeff(i, j, k).Exy2 = (1.0 / Sigma(i, j, k).sigmaE_y - Exy1 / Sigma(i, j, k).sigmaE_y) / Exy1;
						Coeff(i, j, k).Ezy2 = (1.0 / Sigma(i, j, k).sigmaE_y - Ezy1 / Sigma(i, j, k).sigmaE_y) / Ezy1;
					}
					else {
						Coeff(i, j, k).Exy2 = dt;
						Coeff(i, j, k).Ezy2 = dt;
					}
					if (Sigma(i, j, k).sigmaE_z != (ftypePML)0.0)
					{
						Coeff(i, j, k).Exz2 = (1.0 / Sigma(i, j, k).sigmaE_z - Exz1 / Sigma(i, j, k).sigmaE_z) / Exz1;
						Coeff(i, j, k).Eyz2 = (1.0 / Sigma(i, j, k).sigmaE_z - Eyz1 / Sigma(i, j, k).sigmaE_z) / Eyz1;
					}
					else {
						Coeff(i, j, k).Exz2 = dt;
						Coeff(i, j, k).Eyz2 = dt;
					}
					if (Sigma(i, j, k).sigmaH_x != (ftypePML)0.0)
					{
						Coeff(i, j, k).Byx2 = (1.0 / Sigma(i, j, k).sigmaH_x - Byx1 / Sigma(i, j, k).sigmaH_x) / Byx1;
						Coeff(i, j, k).Bzx2 = (1.0 / Sigma(i, j, k).sigmaH_x - Bzx1 / Sigma(i, j, k).sigmaH_x) / Bzx1;
					}
					else {
						Coeff(i, j, k).Byx2 = dt;
						Coeff(i, j, k).Bzx2 = dt;
					}
					if (Sigma(i, j, k).sigmaH_y != (ftypePML)0.0)
					{
						Coeff(i, j, k).Bxy2 = (1.0 / Sigma(i, j, k).sigmaH_y - Bxy1 / Sigma(i, j, k).sigmaH_y) / Bxy1;
						Coeff(i, j, k).Bzy2 = (1.0 / Sigma(i, j, k).sigmaH_y - Bzy1 / Sigma(i, j, k).sigmaH_y) / Bzy1;
					}
					else {
						Coeff(i, j, k).Bxy2 = dt;
						Coeff(i, j, k).Bzy2 = dt;
					}
					if (Sigma(i, j, k).sigmaH_z != (ftypePML)0.0)
					{
						Coeff(i, j, k).Bxz2 = (1.0 / Sigma(i, j, k).sigmaH_z - Bxz1 / Sigma(i, j, k).sigmaH_z) / Bxz1;
						Coeff(i, j, k).Byz2 = (1.0 / Sigma(i, j, k).sigmaH_z - Byz1 / Sigma(i, j, k).sigmaH_z) / Byz1;
					}
					else {
						Coeff(i, j, k).Bxz2 = dt;
						Coeff(i, j, k).Byz2 = dt;
					}
				}
			}
}

template <class ftype>
void Add_Currents_3_dimen(data3d<Component<ftype>>& cube,
	int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z, ftype dx, ftype dy, ftype dz, const ftype dt, int it,
	pair <ftype, ftype> ab, pair <ftype, ftype> cd, pair <ftype, ftype> fg)
{
	ftype x1, x2, y1, y2, z1, z2, t1, t2;

	ftype ax_transverse = ab.second * (ftype)3. / (ftype)4.; // поперечное
	ftype ay_transverse = cd.second * (ftype)3. / (ftype)4.;
	ftype az_transverse = fg.second * (ftype)3. / (ftype)4.;

	ftype tp_x = ab.second / (ftype)6.;
	ftype tp_y = cd.second / (ftype)6.;
	ftype tp_z = fg.second / (ftype)6.;

	ftype lambda_ = (ftype)2. * (ftype)M_PI / (ftype)4.; // wvelength               //?
	ftype k_ = (ftype)2. * (ftype)M_PI / lambda_; // wavenumber = 2 pi/ wavelength
	ftype omega_ = (ftype)2. * (ftype)M_PI / lambda_; // 2 pi c /wavelength 
	ftype A = (ftype)1.; //amplitude

	t1 = dt * (ftype)it;
	t2 = t1 + (ftype)0.5 * dt;

	ftype x0 = (ab.second - ab.first) / (ftype)2.;
	ftype y0 = (cd.second - cd.first) / (ftype)2.;
	ftype z0 = (fg.second - fg.first) / (ftype)2.;
	ftype t0_x = (ftype)3. * tp_x;
	ftype t0_y = (ftype)3. * tp_y;
	ftype t0_z = (ftype)3. * tp_z;

	int index_start_x = delta_x + 1, index_start_y = delta_y + 1, index_start_z = delta_z + 1;

	int offset = 5;
	x1 = (ftype)(ab.first + dx * (ftype)offset + (ftype)0.5 * dx);
	x2 = (ftype)(ab.first + dx * (ftype)offset);

	for (int j = 0; j < Ny; j++)
		for (int k = 0; k < Nz; k++) {
			y1 = (ftype)(cd.first + dy * (ftype)j + (ftype)0.5 * dy);
			y2 = (ftype)(cd.first + dy * (ftype)j);

			z1 = (ftype)(fg.first + dz * (ftype)k + (ftype)0.5 * dz);
			z2 = (ftype)(fg.first + dz * (ftype)k);

			// Ez and By has the same y coordinate
			cube(offset + index_start_x, j + index_start_y, k + index_start_z).Ey += A * exp(-(t1 - t0_x) * (t1 - t0_x) / (tp_x * tp_x))
				* exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse)) *
				sin(omega_ * t1 - k_ * x1); // sin(phase), phase = omega * t - k * x
			cube(offset + index_start_x, j + index_start_y, k + index_start_z).Bz += A * exp(-(t2 - t0_x) * (t2 - t0_x) / (tp_x * tp_x))
				* exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse)) *
				sin(omega_ * t2 - k_ * x2); // sin(phase), phase = omega * t - k * x
		}
}

template <class ftype, class ftypePML>
double FDTD_3D_PML_one_array_scalability(int Nx, int Ny, int Nz, double n, int delta_x, int delta_y, int delta_z, double sigma_x, string type_sum)
{
	setlocale(LC_ALL, "Russian");

	double sigma_y = 1.45312;
	ftype T = 8.0_h * (ftype)M_PI, dt = 0.005_h;
	int Nt = T / dt;

	data3d<Component<ftype>> cube(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	data3d<ComponentSplit<ftypePML>> cube_split(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	data3d<SIGMA<double>> Sigma(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	data3d<COEFF<ftypePML>> Coeff(Nx, Ny, Nz, delta_x, delta_y, delta_z);

	data3d<Component<ftype>> compensator(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	data3d<Component<ftype>> compensator2(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	data3d<ComponentSplit<ftypePML>> compensatorPML(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	data3d<ComponentSplit<ftypePML>> compensatorPML2(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	
	pair<ftype, ftype> ab(0.0, 4.0 * M_PI), cd(0.0, 16.0 * M_PI), fg(0.0, 16.0 * M_PI);
	ftype dx = (ab.second - ab.first) / (ftype)Nx, dy = (cd.second - cd.first) / (ftype)Ny, dz = (fg.second - fg.first) / (ftype)Nz;

	ftype dt_x = dt / (dx), dt_y = dt / (dy), dt_z = dt / (dz);
	ftypePML _1dx = (ftypePML)1.0 / (ftypePML)dx, _1dy = (ftypePML)1.0 / (ftypePML)dy, _1dz = (ftypePML)1.0 / (ftypePML)dz;

	Initializing_cube_split_3_dimen_PML<ftypePML>(cube_split, Nx, Ny, Nz, delta_x, delta_y, delta_z);
	Initializing_Sigma_3_dimen_PML<double>(Sigma, Nx, Ny, Nz, delta_x, delta_y, delta_z, n, sigma_x, sigma_y);
	Initializing_Coeff_3_dimen_PML_correct_2_0<ftypePML>(Coeff, Sigma, Nx, Ny, Nz, delta_x, delta_y, delta_z, (ftypePML)dt);

	double t1, t2;
	t1 = omp_get_wtime();

	if (type_sum == "Kahan")
	{
		for (int it = 0; it < Nt; it++)
		{
			Add_Currents_3_dimen<ftype>(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z, dx, dy, dz, dt, it, ab, cd, fg);

#pragma omp parallel for collapse(3)
			for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
				for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
					for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
					{
						if ((i >= delta_x + 1) && (i < Nx + delta_x + 1) && (j >= delta_y + 1) && (j < Ny + delta_y + 1) && (k >= delta_z + 1) && (k < Nz + delta_z + 1))
						{
							Update_electric_field_three_dimen_Kahan<ftype>(cube, compensator, compensator2, dt_x, dt_y, dt_z, i, j, k);
						}
						else {
							Update_electric_field_three_dimen_PML_Kahan_2_0<ftype, ftypePML>(cube, cube_split, Coeff, compensatorPML, compensatorPML2, compensator, compensator2, _1dx, _1dy, _1dz, i, j, k);
						}
					}
			for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
				for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
					for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
					{
						compensator(i, j, k).Ex = compensator2(i, j, k).Ex;
						compensator(i, j, k).Ey = compensator2(i, j, k).Ey;
						compensator(i, j, k).Ez = compensator2(i, j, k).Ez;

						compensatorPML(i, j, k).Exy = compensatorPML2(i, j, k).Exy;
						compensatorPML(i, j, k).Exz = compensatorPML2(i, j, k).Exz;
						compensatorPML(i, j, k).Eyx = compensatorPML2(i, j, k).Eyx;
						compensatorPML(i, j, k).Eyz = compensatorPML2(i, j, k).Eyz;
						compensatorPML(i, j, k).Ezx = compensatorPML2(i, j, k).Ezx;
						compensatorPML(i, j, k).Ezy = compensatorPML2(i, j, k).Ezy;
					}

#pragma omp parallel for collapse(3)
			for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
				for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
					for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
					{
						if ((i >= delta_x + 1) && (i < Nx + delta_x + 1) && (j >= delta_y + 1) && (j < Ny + delta_y + 1) && (k >= delta_z + 1) && (k < Nz + delta_z + 1))
						{
							Update_magnetic_field_three_dimen_Kahan<ftype>(cube, compensator, compensator2, dt_x, dt_y, dt_z, i, j, k);
						}
						else {
							Update_magnetic_field_three_dimen_PML_Kahan_2_0<ftype, ftypePML>(cube, cube_split, Coeff, compensatorPML, compensatorPML2, compensator, compensator2, _1dx, _1dy, _1dz, i, j, k);
						}
					}
			for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
				for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
					for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
					{
						compensator(i, j, k).Bx = compensator2(i, j, k).Bx;
						compensator(i, j, k).By = compensator2(i, j, k).By;
						compensator(i, j, k).Bz = compensator2(i, j, k).Bz;

						compensatorPML(i, j, k).Bxy = compensatorPML2(i, j, k).Bxy;
						compensatorPML(i, j, k).Bxz = compensatorPML2(i, j, k).Bxz;
						compensatorPML(i, j, k).Byx = compensatorPML2(i, j, k).Byx;
						compensatorPML(i, j, k).Byz = compensatorPML2(i, j, k).Byz;
						compensatorPML(i, j, k).Bzx = compensatorPML2(i, j, k).Bzx;
						compensatorPML(i, j, k).Bzy = compensatorPML2(i, j, k).Bzy;
					}
		}
	}
	else
	{
		for (int it = 0; it < Nt; it++)
		{
			Add_Currents_3_dimen<ftype>(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z, dx, dy, dz, dt, it, ab, cd, fg);
			
#pragma omp parallel for collapse(3)
			for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
				for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
					for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
					{
						if ((i >= delta_x + 1) && (i < Nx + delta_x + 1) && (j >= delta_y + 1) && (j < Ny + delta_y + 1) && (k >= delta_z + 1) && (k < Nz + delta_z + 1))
						{
							Update_electric_field_three_dimen<ftype>(cube, dt_x, dt_y, dt_z, i, j, k);
						}
						else {
							Update_electric_field_three_dimen_PML_2_0<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
						}
					}

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
							Update_magnetic_field_three_dimen_PML_2_0<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
						}
					}
		}
	}

	t2 = omp_get_wtime();
	double time = t2 - t1;

	return time;
}
