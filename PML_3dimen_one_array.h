#define _USE_MATH_DEFINES
#include <iostream>
#include <chrono>
#include <vector>
#include <math.h>
#include <cmath>
#include <omp.h>
#include "half.hpp"
#include "formula PML 3D 1_array.h"
#include "optional features 1_array.h"
#include "type data.h"

using namespace std;
using half_float::half;
using namespace half_float::literal;




template <class ftype>
void Initializing_FDTD_3_dimen_Gauss_PML(data3d<Component<ftype>>& cube,
	int Nx, int Ny, int Nz, int delta, ftype dx, ftype dy, ftype dz, const ftype dt, int index_start, Direction direction,
	pair <ftype, ftype> ab, pair <ftype, ftype> cd, pair <ftype, ftype> fg)
{
	ftype x1, x2, y1, y2, z1, z2, t;
	//ftype ax = 3. * M_PI / 3.0; // size of beam in x direction // продольное направление
	//ftype ay = 8. * M_PI / 3.0; // size of beam in y direction // поперечное направление
	//ftype az = 8. * M_PI / 3.0; // size of beam in z direction  // поперечное направление

	ftype ax_longitudinal = ab.second / 6.; // продольное
	ftype ay_longitudinal = cd.second / 6.;
	ftype az_longitudinal = fg.second / 6.;
	ftype ax_transverse = ab.second * 3. / 4.; // поперечное
	ftype ay_transverse = cd.second * 3. / 4.;
	ftype az_transverse = fg.second * 3. / 4.;

	double lambda_ = 2. * M_PI / 4.; // wvelength               //?
	double k_ = 2. * M_PI / lambda_; // wavenumber = 2 pi/ wavelength
	double omega_ = 2. * M_PI / lambda_; // 2 pi c /wavelength 
	double t0_ = 0; // initial moment of time
	t = (ftype)0.5 * dt; // timeshift between components

	float x0 = (ab.second - ab.first) / 2.;
	float y0 = (cd.second - cd.first) / 2.;
	float z0 = (fg.second - fg.first) / 2.;

	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			for (int k = 0; k < Nz; k++) {
				x1 = (ftype)(ab.first + dx * (ftype)i + (ftype)0.5 * dx);
				x2 = (ftype)(ab.first + dx * (ftype)i);

				y1 = (ftype)(cd.first + dy * (ftype)j + (ftype)0.5 * dy);
				y2 = (ftype)(cd.first + dy * (ftype)j);

				z1 = (ftype)(fg.first + dz * (ftype)k + (ftype)0.5 * dz);
				z2 = (ftype)(fg.first + dz * (ftype)k);

				switch (direction)
				{
				case x_comp_zy: {
					// Ez and By has the same y coordinate
					cube(i + index_start, j + index_start, k + index_start).Ez = exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse))
						* exp(-(x1 - x0) * (x1 - x0) / (ax_longitudinal * ax_longitudinal)) * exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse)) *
						sin(omega_ * t0_ - k_ * x1); // sin(phase), phase = omega * t - k * x
					cube(i + index_start, j + index_start, k + index_start).By = exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_longitudinal * ax_longitudinal)) * exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse)) *
						sin(omega_ * (t0_ + t) - k_ * x2); // sin(phase), phase = omega * t - k * x
					break;
				}
				case x_comp_yz: {
					// Ey and Bz has the same x coordinate
					cube(i + index_start, j + index_start, k + index_start).Ey = exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse))
						* exp(-(x1 - x0) * (x1 - x0) / (ax_longitudinal * ax_longitudinal)) * exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse)) *
						sin(omega_ * t0_ - k_ * x1); // sin(phase), phase = omega * t - k * x
					cube(i + index_start, j + index_start, k + index_start).Bz = exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_longitudinal * ax_longitudinal)) * exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse)) *
						sin(omega_ * (t0_ + t) - k_ * x2); // sin(phase), phase = omega * t - k * x
					break;
				}
				case y_comp_xz: {
					// Ex and Bz has the same y coordinate
					cube(i + index_start, j + index_start, k + index_start).Ex = exp(-(y1 - y0) * (y1 - y0) / (ay_longitudinal * ay_longitudinal))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse)) * exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse)) *
						sin(omega_ * t0_ - k_ * y1); // sin(phase), phase = omega * t - k * x
					cube(i + index_start, j + index_start, k + index_start).Bz = exp(-(y2 - y0) * (y2 - y0) / (ay_longitudinal * ay_longitudinal))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse)) * exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse)) *
						sin(omega_ * (t0_ + t) - k_ * y2); // sin(phase), phase = omega * t - k * x
					break;
				}
				case y_comp_zx: {
					// Ez and Bx has the same y coordinate
					cube(i + index_start, j + index_start, k + index_start).Ez = exp(-(y1 - y0) * (y1 - y0) / (ay_longitudinal * ay_longitudinal))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse)) * exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse)) *
						sin(omega_ * t0_ - k_ * y1); // sin(phase), phase = omega * t - k * x
					cube(i + index_start, j + index_start, k + index_start).Bx = exp(-(y2 - y0) * (y2 - y0) / (ay_longitudinal * ay_longitudinal))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse)) * exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse)) *
						sin(omega_ * (t0_ + t) - k_ * y2); // sin(phase), phase = omega * t - k * x
					break;
				}
				case z_comp_xy: {
					// Ex and By has the same z coordinate
					cube(i + index_start, j + index_start, k + index_start).Ex = -exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse)) * exp(-(z1 - z0) * (z1 - z0) / (az_longitudinal * az_longitudinal)) *
						sin(omega_ * t0_ + k_ * z1); // sin(phase), phase = omega * t - k * x
					cube(i + index_start, j + index_start, k + index_start).By = exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse)) * exp(-(z2 - z0) * (z2 - z0) / (az_longitudinal * az_longitudinal)) *
						sin(omega_ * (t0_ + t) + k_ * z2); // sin(phase), phase = omega * t - k * x
					break;
				}
				case z_comp_yx: {
					// Ey and Bx has the same z coordinate
					cube(i + index_start, j + index_start, k + index_start).Ey = exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse)) * exp(-(z1 - z0) * (z1 - z0) / (az_longitudinal * az_longitudinal)) *
						sin(omega_ * t0_ - k_ * z1); // sin(phase), phase = omega * t - k * x
					cube(i + index_start, j + index_start, k + index_start).Bx = exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse)) * exp(-(z2 - z0) * (z2 - z0) / (az_longitudinal * az_longitudinal)) *
						sin(omega_ * (t0_ + t) - k_ * z2); // sin(phase), phase = omega * t - k * x
					break;
				}
				default:
					break;
				}
	}
	// Border conditions
	for (int i = 0; i < Ny + 2 * delta + 2; i++)
		for (int j = 0; j < Nz + 2 * delta + 2; j++) {
		cube(0, i, j).Ex = cube(Nx + 2 * delta, i, j).Ex;
		cube(0, i, j).Ey = cube(Nx + 2 * delta, i, j).Ey;
		cube(0, i, j).Ez = cube(Nx + 2 * delta, i, j).Ez;
		cube(0, i, j).Bx = cube(Nx + 2 * delta, i, j).Bx;
		cube(0, i, j).By = cube(Nx + 2 * delta, i, j).By;
		cube(0, i, j).Bz = cube(Nx + 2 * delta, i, j).Bz;
		cube(Nx + 2 * delta + 1, i, j).Ex = cube(1, i, j).Ex;
		cube(Nx + 2 * delta + 1, i, j).Ey = cube(1, i, j).Ey;
		cube(Nx + 2 * delta + 1, i, j).Ez = cube(1, i, j).Ez;
		cube(Nx + 2 * delta + 1, i, j).Bx = cube(1, i, j).Bx;
		cube(Nx + 2 * delta + 1, i, j).By = cube(1, i, j).By;
		cube(Nx + 2 * delta + 1, i, j).Bz = cube(1, i, j).Bz;
	}
	for (int i = 0; i < Nx + 2 * delta + 2; i++)
		for (int j = 0; j < Nz + 2 * delta + 2; j++) {
		cube(i, 0, j).Ex = cube(i, Ny + 2 * delta, j).Ex;
		cube(i, 0, j).Ey = cube(i, Ny + 2 * delta, j).Ey;
		cube(i, 0, j).Ez = cube(i, Ny + 2 * delta, j).Ez;
		cube(i, 0, j).Bx = cube(i, Ny + 2 * delta, j).Bx;
		cube(i, 0, j).By = cube(i, Ny + 2 * delta, j).By;
		cube(i, 0, j).Bz = cube(i, Ny + 2 * delta, j).Bz;
		cube(i, Ny + 2 * delta + 1, j).Ex = cube(i, 1, j).Ex;
		cube(i, Ny + 2 * delta + 1, j).Ey = cube(i, 1, j).Ey;
		cube(i, Ny + 2 * delta + 1, j).Ez = cube(i, 1, j).Ez;
		cube(i, Ny + 2 * delta + 1, j).Bx = cube(i, 1, j).Bx;
		cube(i, Ny + 2 * delta + 1, j).By = cube(i, 1, j).By;
		cube(i, Ny + 2 * delta + 1, j).Bz = cube(i, 1, j).Bz;
	}
	for (int i = 0; i < Nx + 2 * delta + 2; i++)
		for (int j = 0; j < Ny + 2 * delta + 2; j++) {
			cube(i, j, 0).Ex = cube(i, j, Nz + 2 * delta).Ex;
			cube(i, j, 0).Ey = cube(i, j, Nz + 2 * delta).Ey;
			cube(i, j, 0).Ez = cube(i, j, Nz + 2 * delta).Ez;
			cube(i, j, 0).Bx = cube(i, j, Nz + 2 * delta).Bx;
			cube(i, j, 0).By = cube(i, j, Nz + 2 * delta).By;
			cube(i, j, 0).Bz = cube(i, j, Nz + 2 * delta).Bz;
			cube(i, j, Nz + 2 * delta + 1).Ex = cube(i, j, 1).Ex;
			cube(i, j, Nz + 2 * delta + 1).Ey = cube(i, j, 1).Ey;
			cube(i, j, Nz + 2 * delta + 1).Ez = cube(i, j, 1).Ez;
			cube(i, j, Nz + 2 * delta + 1).Bx = cube(i, j, 1).Bx;
			cube(i, j, Nz + 2 * delta + 1).By = cube(i, j, 1).By;
			cube(i, j, Nz + 2 * delta + 1).Bz = cube(i, j, 1).Bz;
		}
}

template <class ftypePML>
void Initializing_cube_split_3_dimen_PML(data3d<ComponentSplit<ftypePML>>& cube_split, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	for (int i = 0; i < Nx + 2 * delta_x +2; i++)
		for (int j = 0; j < Ny + 2 * delta_y +2; j++)
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

	//cout << "delta_x= "<< delta_x << "  delta_y= " << delta_y << "   delta_z= " << delta_z << "    " << endl;
	//cout << "Nx= " << Nx << "  Ny= " << Ny << "   Nz= " << Nz << "    " << endl;
	//cout << "n = " << n << endl;
	cout << var_Sigma_max_x << "  " << var_Sigma_max_y << "   " << var_Sigma_max_z << "    " << endl;
	
	//cout << func_R(n, var_Sigma_max_x, delta, dx) << endl;

	for(int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++) {
				Sigma(i, j, k).sigmaE_x = var_Sigma_max_x * powf(distanceE<ftypePML>(Nx, delta_x, i), n);
				Sigma(i, j, k).sigmaH_x = var_Sigma_max_x * powf(distanceH<ftypePML>(Nx, delta_x, i), n);
	
				Sigma(i, j, k).sigmaE_y = var_Sigma_max_y * powf(distanceE<ftypePML>(Ny, delta_y, j), n);
				Sigma(i, j, k).sigmaH_y = var_Sigma_max_y * powf(distanceH<ftypePML>(Ny, delta_y, j), n);

				Sigma(i, j, k).sigmaE_z = var_Sigma_max_z * powf(distanceE<ftypePML>(Nz, delta_z, k), n);
				Sigma(i, j, k).sigmaH_z = var_Sigma_max_z * powf(distanceH<ftypePML>(Nz, delta_z, k), n);
			}
	//cout << sizeof(Sigma(10, 18, 11).sigmaH_x) << sizeof(double) << sizeof(float) << endl;
}

template <class ftypePML>
void Initializing_Coeff_3_dimen_PML_old(data3d<COEFF<ftypePML>>& Coeff, data3d<SIGMA<double>>& Sigma, int Nx, int Ny, int Nz,
	int delta_x, int delta_y, int delta_z, ftypePML dt)
{
	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++) {
				if ((i >= delta_x + 1) && (i < Nx + delta_x + 1) && (j >= delta_y + 1) && (j < Ny + delta_y + 1) && (k >= delta_z + 1) && (k < Nz + delta_z + 1))
				{
				}
				else {
					Coeff(i, j, k).Exy1 = (ftypePML)exp(-dt * Sigma(i, j, k).sigmaE_y);
					Coeff(i, j, k).Exz1 = (ftypePML)exp(-dt * Sigma(i, j, k).sigmaE_z);
					Coeff(i, j, k).Eyx1 = (ftypePML)exp(-dt * Sigma(i, j, k).sigmaE_x);
					Coeff(i, j, k).Eyz1 = (ftypePML)exp(-dt * Sigma(i, j, k).sigmaE_z);
					Coeff(i, j, k).Ezx1 = (ftypePML)exp(-dt * Sigma(i, j, k).sigmaE_x);
					Coeff(i, j, k).Ezy1 = (ftypePML)exp(-dt * Sigma(i, j, k).sigmaE_y);

					Coeff(i, j, k).Bxy1 = (ftypePML)exp(-dt * Sigma(i, j, k).sigmaH_y);
					Coeff(i, j, k).Bxz1 = (ftypePML)exp(-dt * Sigma(i, j, k).sigmaH_z);
					Coeff(i, j, k).Byx1 = (ftypePML)exp(-dt * Sigma(i, j, k).sigmaH_x);
					Coeff(i, j, k).Byz1 = (ftypePML)exp(-dt * Sigma(i, j, k).sigmaH_z);
					Coeff(i, j, k).Bzx1 = (ftypePML)exp(-dt * Sigma(i, j, k).sigmaH_x);
					Coeff(i, j, k).Bzy1 = (ftypePML)exp(-dt * Sigma(i, j, k).sigmaH_y);

					if (Sigma(i, j, k).sigmaE_x != (ftypePML)0.0) {
						Coeff(i, j, k).Eyx2 = (ftypePML)1.0 / Sigma(i, j, k).sigmaE_x - Coeff(i, j, k).Eyx1 / Sigma(i, j, k).sigmaE_x;
						Coeff(i, j, k).Ezx2 = (ftypePML)1.0 / Sigma(i, j, k).sigmaE_x - Coeff(i, j, k).Ezx1 / Sigma(i, j, k).sigmaE_x;
					}
					else {
						Coeff(i, j, k).Eyx2 = dt;
						Coeff(i, j, k).Ezx2 = dt;
					}
					if (Sigma(i, j, k).sigmaE_y != (ftypePML)0.0) {
						Coeff(i, j, k).Exy2 = (ftypePML)1.0 / Sigma(i, j, k).sigmaE_y - Coeff(i, j, k).Exy1 / Sigma(i, j, k).sigmaE_y;
						Coeff(i, j, k).Ezy2 = (ftypePML)1.0 / Sigma(i, j, k).sigmaE_y - Coeff(i, j, k).Ezy1 / Sigma(i, j, k).sigmaE_y;
					}
					else {
						Coeff(i, j, k).Exy2 = dt;
						Coeff(i, j, k).Ezy2 = dt;
					}
					if (Sigma(i, j, k).sigmaE_z != (ftypePML)0.0)
					{
						Coeff(i, j, k).Exz2 = (ftypePML)1.0 / Sigma(i, j, k).sigmaE_z - Coeff(i, j, k).Exz1 / Sigma(i, j, k).sigmaE_z;
						Coeff(i, j, k).Eyz2 = (ftypePML)1.0 / Sigma(i, j, k).sigmaE_z - Coeff(i, j, k).Eyz1 / Sigma(i, j, k).sigmaE_z;
					}
					else {
						Coeff(i, j, k).Exz2 = dt;
						Coeff(i, j, k).Eyz2 = dt;
					}
					if (Sigma(i, j, k).sigmaH_x != (ftypePML)0.0)
					{
						Coeff(i, j, k).Byx2 = (ftypePML)1.0 / Sigma(i, j, k).sigmaH_x - Coeff(i, j, k).Byx1 / Sigma(i, j, k).sigmaH_x;
						Coeff(i, j, k).Bzx2 = (ftypePML)1.0 / Sigma(i, j, k).sigmaH_x - Coeff(i, j, k).Bzx1 / Sigma(i, j, k).sigmaH_x;
					}
					else {
						Coeff(i, j, k).Byx2 = dt;
						Coeff(i, j, k).Bzx2 = dt;
					}
					if (Sigma(i, j, k).sigmaH_y != (ftypePML)0.0)
					{
						Coeff(i, j, k).Bxy2 = (ftypePML)1.0 / Sigma(i, j, k).sigmaH_y - Coeff(i, j, k).Bxy1 / Sigma(i, j, k).sigmaH_y;
						Coeff(i, j, k).Bzy2 = (ftypePML)1.0 / Sigma(i, j, k).sigmaH_y - Coeff(i, j, k).Bzy1 / Sigma(i, j, k).sigmaH_y;
					}
					else {
						Coeff(i, j, k).Bxy2 = dt;
						Coeff(i, j, k).Bzy2 = dt;
					}
					if (Sigma(i, j, k).sigmaH_z != (ftypePML)0.0)
					{
						Coeff(i, j, k).Bxz2 = (ftypePML)1.0 / Sigma(i, j, k).sigmaH_z - Coeff(i, j, k).Bxz1 / Sigma(i, j, k).sigmaH_z;
						Coeff(i, j, k).Byz2 = (ftypePML)1.0 / Sigma(i, j, k).sigmaH_z - Coeff(i, j, k).Byz1 / Sigma(i, j, k).sigmaH_z;
					}
					else {
						Coeff(i, j, k).Bxz2 = dt;
						Coeff(i, j, k).Byz2 = dt;
					}
				}
		}
	// cout << sizeof(Coeff(1, 5, 4).Bxy2) << sizeof(Coeff(2, 4, 3).Exz1) << endl;

}

template <class ftypePML, class ftype = double>
void Initializing_Coeff_3_dimen_PML_correct(data3d<COEFF<ftypePML>>& Coeff, data3d<SIGMA<double>>& Sigma, int Nx, int Ny, int Nz,
	int delta_x, int delta_y, int delta_z, ftype dt) {

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
						Coeff(i, j, k).Eyx2 = 1.0 / Sigma(i, j, k).sigmaE_x - Eyx1 / Sigma(i, j, k).sigmaE_x;
						Coeff(i, j, k).Ezx2 = 1.0 / Sigma(i, j, k).sigmaE_x - Ezx1 / Sigma(i, j, k).sigmaE_x;
					}
					else {
						Coeff(i, j, k).Eyx2 = dt;
						Coeff(i, j, k).Ezx2 = dt;
					}
					if (Sigma(i, j, k).sigmaE_y != (ftypePML)0.0) {
						Coeff(i, j, k).Exy2 = 1.0 / Sigma(i, j, k).sigmaE_y - Exy1 / Sigma(i, j, k).sigmaE_y;
						Coeff(i, j, k).Ezy2 = 1.0 / Sigma(i, j, k).sigmaE_y - Ezy1 / Sigma(i, j, k).sigmaE_y;
					}
					else {
						Coeff(i, j, k).Exy2 = dt;
						Coeff(i, j, k).Ezy2 = dt;
					}
					if (Sigma(i, j, k).sigmaE_z != (ftypePML)0.0)
					{
						Coeff(i, j, k).Exz2 = 1.0 / Sigma(i, j, k).sigmaE_z - Exz1 / Sigma(i, j, k).sigmaE_z;
						Coeff(i, j, k).Eyz2 = 1.0 / Sigma(i, j, k).sigmaE_z - Eyz1 / Sigma(i, j, k).sigmaE_z;
					}
					else {
						Coeff(i, j, k).Exz2 = dt;
						Coeff(i, j, k).Eyz2 = dt;
					}
					if (Sigma(i, j, k).sigmaH_x != (ftypePML)0.0)
					{
						Coeff(i, j, k).Byx2 = 1.0 / Sigma(i, j, k).sigmaH_x - Byx1 / Sigma(i, j, k).sigmaH_x;
						Coeff(i, j, k).Bzx2 = 1.0 / Sigma(i, j, k).sigmaH_x - Bzx1 / Sigma(i, j, k).sigmaH_x;
					}
					else {
						Coeff(i, j, k).Byx2 = dt;
						Coeff(i, j, k).Bzx2 = dt;
					}
					if (Sigma(i, j, k).sigmaH_y != (ftypePML)0.0)
					{
						Coeff(i, j, k).Bxy2 = 1.0 / Sigma(i, j, k).sigmaH_y - Bxy1 / Sigma(i, j, k).sigmaH_y;
						Coeff(i, j, k).Bzy2 = 1.0 / Sigma(i, j, k).sigmaH_y - Bzy1 / Sigma(i, j, k).sigmaH_y;
					}
					else {
						Coeff(i, j, k).Bxy2 = dt;
						Coeff(i, j, k).Bzy2 = dt;
					}
					if (Sigma(i, j, k).sigmaH_z != (ftypePML)0.0)
					{
						Coeff(i, j, k).Bxz2 = 1.0 / Sigma(i, j, k).sigmaH_z - Bxz1 / Sigma(i, j, k).sigmaH_z;
						Coeff(i, j, k).Byz2 = 1.0 / Sigma(i, j, k).sigmaH_z - Byz1 / Sigma(i, j, k).sigmaH_z;
					}
					else {
						Coeff(i, j, k).Bxz2 = dt;
						Coeff(i, j, k).Byz2 = dt;
					}
				}
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
						Coeff(i, j, k).Exy2 = (1.0 / Sigma(i, j, k).sigmaE_y - Exy1 / Sigma(i, j, k).sigmaE_y)/Exy1;
						Coeff(i, j, k).Ezy2 = (1.0 / Sigma(i, j, k).sigmaE_y - Ezy1 / Sigma(i, j, k).sigmaE_y)/Ezy1;
					}
					else {
						Coeff(i, j, k).Exy2 = dt;
						Coeff(i, j, k).Ezy2 = dt;
					}
					if (Sigma(i, j, k).sigmaE_z != (ftypePML)0.0)
					{
						Coeff(i, j, k).Exz2 = (1.0 / Sigma(i, j, k).sigmaE_z - Exz1 / Sigma(i, j, k).sigmaE_z)/Exz1;
						Coeff(i, j, k).Eyz2 = (1.0 / Sigma(i, j, k).sigmaE_z - Eyz1 / Sigma(i, j, k).sigmaE_z)/Eyz1;
					}
					else {
						Coeff(i, j, k).Exz2 = dt;
						Coeff(i, j, k).Eyz2 = dt;
					}
					if (Sigma(i, j, k).sigmaH_x != (ftypePML)0.0)
					{
						Coeff(i, j, k).Byx2 = (1.0 / Sigma(i, j, k).sigmaH_x - Byx1 / Sigma(i, j, k).sigmaH_x)/Byx1;
						Coeff(i, j, k).Bzx2 = (1.0 / Sigma(i, j, k).sigmaH_x - Bzx1 / Sigma(i, j, k).sigmaH_x)/Bzx1;
					}
					else {
						Coeff(i, j, k).Byx2 = dt;
						Coeff(i, j, k).Bzx2 = dt;
					}
					if (Sigma(i, j, k).sigmaH_y != (ftypePML)0.0)
					{
						Coeff(i, j, k).Bxy2 = (1.0 / Sigma(i, j, k).sigmaH_y - Bxy1 / Sigma(i, j, k).sigmaH_y)/Bxy1;
						Coeff(i, j, k).Bzy2 = (1.0 / Sigma(i, j, k).sigmaH_y - Bzy1 / Sigma(i, j, k).sigmaH_y)/Bzy1;
					}
					else {
						Coeff(i, j, k).Bxy2 = dt;
						Coeff(i, j, k).Bzy2 = dt;
					}
					if (Sigma(i, j, k).sigmaH_z != (ftypePML)0.0)
					{
						Coeff(i, j, k).Bxz2 = (1.0 / Sigma(i, j, k).sigmaH_z - Bxz1 / Sigma(i, j, k).sigmaH_z)/Bxz1;
						Coeff(i, j, k).Byz2 = (1.0 / Sigma(i, j, k).sigmaH_z - Byz1 / Sigma(i, j, k).sigmaH_z)/Byz1;
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

	//ftype ax_transverse = ab.second * (ftype)1. / (ftype)2.; // поперечное
	//ftype ay_transverse = cd.second * (ftype)1. / (ftype)2.;
	//ftype az_transverse = fg.second * (ftype)1. / (ftype)2.;

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
double FDTD_3D_PML_one_array(double n, int Nx, int Ny, int Nz, ftype T, ftype dt, int delta_x, int delta_y, int delta_z, double sigma_x, double sigma_y,
	string type_sum, string file_energy, string file_data, string file_sigma)
{
	setlocale(LC_ALL, "Russian");

	int Nt = T / dt;
	// delta  =  width boundary layer

	data3d<Component<ftype>> cube(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	data3d<ComponentSplit<ftypePML>> cube_split(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	data3d<SIGMA<double>> Sigma(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	data3d<COEFF<ftypePML>> Coeff(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	data3d<COEFF<half>> CoeffHalf(Nx, Ny, Nz, delta_x, delta_y, delta_z);

	data3d<Component<ftype>> compensator(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	data3d<Component<ftype>> compensator2(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	data3d<ComponentSplit<ftypePML>> compensatorPML(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	data3d<ComponentSplit<ftypePML>> compensatorPML2(Nx, Ny, Nz, delta_x, delta_y, delta_z);

	pair<ftype, ftype> ab(0.0, 4.0 * M_PI), cd(0.0, 16.0 * M_PI), fg(0.0, 16.0 * M_PI);
	ftype dx = (ab.second - ab.first) / (ftype)Nx, dy = (cd.second - cd.first) / (ftype)Ny, dz = (fg.second - fg.first) / (ftype)Nz;

	ftype dt_x = dt / (dx), dt_y = dt / (dy), dt_z = dt / (dz);
	ftypePML _1dx = (ftypePML)1.0 / dx, _1dy = (ftypePML)1.0 / dy, _1dz = (ftypePML)1.0 / dz;
	vector<double> vec_energy(Nt + 1), vec_energy_acc(Nt + 1);

	//vector<double> vec_max_dataE(Nt + 1), vec_max_dataB(Nt + 1);
	//vector<double> vec_average_dataE(Nt + 1), vec_average_dataB(Nt + 1);

	//vector<double> vec_dataEy(Nt + 1), vec_dataBz(Nt + 1);


	// Initializing_FDTD_3_dimen_Gauss_PML<ftype>(cube, Nx, Ny, Nz, delta, dx, dy, dz, dt, delta + 1, direction, ab, cd, fg);
	Initializing_cube_split_3_dimen_PML<ftypePML>(cube_split, Nx, Ny, Nz, delta_x, delta_y, delta_z);
	Initializing_Sigma_3_dimen_PML<double>(Sigma, Nx, Ny, Nz, delta_x, delta_y, delta_z, n, sigma_x, sigma_y);
	
	Initializing_Coeff_3_dimen_PML_correct<ftypePML>(Coeff, Sigma, Nx, Ny, Nz, delta_x, delta_y, delta_z, dt);

	//Graph_for_Sigma_three_dimen<ftypePML>(Sigma, Nx, Ny, Nz, delta_x, delta_y, delta_z, dx, file_sigma);

//	Initializing_Coeff_3_dimen_PML_correct<half>(CoeffHalf, Sigma, Nx, Ny, Nz, delta, (half)dt);
//	for (int i = 0; i < Nx + 2 * delta + 2; i++)
//		for (int j = 0; j < Ny + 2 * delta + 2; j++)
//			for (int k = 0; k < Nz + 2 * delta + 2; k++)
//			{
//				Coeff(i, j, k).Exy1 = ((CoeffHalf(i, j, k).Exy1));
//				Coeff(i, j, k).Exz1 = ((CoeffHalf(i, j, k).Exz1));
//				Coeff(i, j, k).Eyx1 = ((CoeffHalf(i, j, k).Eyx1));
//				Coeff(i, j, k).Eyz1 = ((CoeffHalf(i, j, k).Eyz1));
//				Coeff(i, j, k).Ezx1 = ((CoeffHalf(i, j, k).Ezx1));
//				Coeff(i, j, k).Ezy1 = ((CoeffHalf(i, j, k).Ezy1));
//
//				Coeff(i, j, k).Bxy1 = ((CoeffHalf(i, j, k).Bxy1));
//				Coeff(i, j, k).Bxz1 = ((CoeffHalf(i, j, k).Bxz1));
//				Coeff(i, j, k).Byx1 = ((CoeffHalf(i, j, k).Byx1));
//				Coeff(i, j, k).Byz1 = ((CoeffHalf(i, j, k).Byz1));
//				Coeff(i, j, k).Bzx1 = ((CoeffHalf(i, j, k).Bzx1));
//				Coeff(i, j, k).Bzy1 = ((CoeffHalf(i, j, k).Bzy1));
//
//				Coeff(i, j, k).Exy2 = ((CoeffHalf(i, j, k).Exy2));
//				Coeff(i, j, k).Exz2 = ((CoeffHalf(i, j, k).Exz2));
//				Coeff(i, j, k).Eyx2 = ((CoeffHalf(i, j, k).Eyx2));
//				Coeff(i, j, k).Eyz2 = ((CoeffHalf(i, j, k).Eyz2));
//				Coeff(i, j, k).Ezx2 = ((CoeffHalf(i, j, k).Ezx2));
//				Coeff(i, j, k).Ezy2 = ((CoeffHalf(i, j, k).Ezy2));
//
//				Coeff(i, j, k).Bxy2 = ((CoeffHalf(i, j, k).Bxy2));
//				Coeff(i, j, k).Bxz2 = ((CoeffHalf(i, j, k).Bxz2));
//				Coeff(i, j, k).Byx2 = ((CoeffHalf(i, j, k).Byx2));
//				Coeff(i, j, k).Byz2 = ((CoeffHalf(i, j, k).Byz2));
//				Coeff(i, j, k).Bzx2 = ((CoeffHalf(i, j, k).Bzx2));
//				Coeff(i, j, k).Bzy2 = ((CoeffHalf(i, j, k).Bzy2));
////			}
////		}
//// cout << sizeof(Coeff(1, 5, 4).Bxy2) << sizeof(Coeff(2, 4, 3).Exz1) << endl;
//			}

	 double t1, t2;
	 t1 = omp_get_wtime();

	if (type_sum == "Kahan")
	{
		Initializing_Coeff_3_dimen_PML_correct_2_0<ftypePML>(Coeff, Sigma, Nx, Ny, Nz, delta_x, delta_y, delta_z, dt);

		cout << type_sum << endl;
		for (int it = 0; it < Nt; it++)
		{
			//if (it % 1000 == 0) cout << it << endl;

			Add_Currents_3_dimen<ftype>(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z, dx, dy, dz, dt, it, ab, cd, fg);

			vec_energy[it] = CalculateEnergyPML_3dimen(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
			vec_energy_acc[it] = CalculateEnergyPML_3dimen_accurate(cube, compensator, Nx, Ny, Nz, delta_x, delta_y, delta_z);

			//vec_average_dataE[it] = AverageValueElectric(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
			//vec_average_dataB[it] = AverageValueMagnetic(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
			//vec_max_dataE[it] = MaxValueElectric(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
			//vec_max_dataB[it] = MaxValueMagnetic(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);

			//vec_dataEy[it] = (double)cube(Nx/2 + delta_x, Ny/2 + delta_y, Nz/2 + delta_z+1).Ey;
			//vec_dataBz[it] = (double)cube(Nx / 2 + delta_x, Ny / 2 + delta_y, Nz / 2 + delta_z+1).Bz;

			// cout << cube(Nx, Ny, Nz).Ey << "  " << (double)cube(Nx + 2 * delta - 5, Ny + 2 * delta - 5, Nz + 2 * delta - 5).Bz << endl;

					   
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
							//Update_electric_field_three_dimen_PML<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
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
							//Update_magnetic_field_three_dimen_PML<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
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
		Initializing_Coeff_3_dimen_PML_correct<ftypePML>(Coeff, Sigma, Nx, Ny, Nz, delta_x, delta_y, delta_z, dt);

		cout << type_sum << endl;
		for (int it = 0; it < Nt; it++)
		{
			if (it % 1000 == 0)
				cout << it << endl;
			//if(it== (int)(Nt*0.25))
			//	Graph_Solution_in_two_planes_3dimen(cube, Nx, Ny, Nz, delta, dx, dy, dz, direction, "Bz PML T=8pi 1.csv");
			//if (it == (int)(Nt * 0.5))
			//	Graph_Solution_in_two_planes_3dimen(cube, Nx, Ny, Nz, delta, dx, dy, dz, direction, "Bz PML T=8pi 2.csv");
			//if (it == (int)(Nt * 0.75))
			//	Graph_Solution_in_two_planes_3dimen(cube, Nx, Ny, Nz, delta, dx, dy, dz, direction, "Bz PML T=8pi 3.csv");
			//if (it == 2188)
			//	Graph_Solution_in_two_planes_3dimen(cube, Nx, Ny, Nz, delta, dx, dy, dz, direction, "Bz PML T=8pi 3-5.csv");
			
			 Add_Currents_3_dimen<ftype>(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z, dx, dy, dz, dt, it, ab, cd, fg);
			 vec_energy[it] = CalculateEnergyPML_3dimen(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);

			 //vec_average_dataE[0] = AverageValueElectric(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
			 //vec_average_dataB[0] = AverageValueMagnetic(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
			 //vec_max_dataE[0] = MaxValueElectric(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
			 //vec_max_dataB[0] = MaxValueMagnetic(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);

#pragma omp parallel for collapse(3)
			for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
				for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
					for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
					{
						if ((i >= delta_x + 1) && (i < Nx + delta_x + 1)&& (j >= delta_y + 1) && (j < Ny + delta_y + 1) && (k >= delta_z + 1) && (k < Nz + delta_z + 1))
						{
							Update_electric_field_three_dimen<ftype>(cube, dt_x, dt_y, dt_z, i, j, k);
						}
						else {
							Update_electric_field_three_dimen_PML<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
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
							Update_magnetic_field_three_dimen_PML<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
						}
					}
			

			//vec_average_dataE[it + 1] = AverageValueElectric(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
			//vec_average_dataB[it + 1] = AverageValueMagnetic(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
			//vec_max_dataE[it + 1] = MaxValueElectric(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
			//vec_max_dataB[it + 1] = MaxValueMagnetic(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
		}
	}

	t2 = omp_get_wtime();


	vec_energy[Nt] = CalculateEnergyPML_3dimen(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
	vec_energy_acc[Nt] = CalculateEnergyPML_3dimen_accurate(cube, compensator, Nx, Ny, Nz, delta_x, delta_y, delta_z);

	double result = vec_energy[Nt] / vec_energy[2512];

	//vec_average_dataE[Nt] = AverageValueElectric(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
	//vec_average_dataB[Nt] = AverageValueMagnetic(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
	//vec_max_dataE[Nt] = MaxValueElectric(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
	//vec_max_dataB[Nt] = MaxValueMagnetic(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);


	// Graph_Solution_in_one_planes_3dimen(cube, Nx, Ny, Nz, delta, dx, dy, dz, direction);

	cout << "Reflection coefficient = "<< endl << vec_energy[Nt] / vec_energy[2512] << endl << "Accurate reflection coefficient = " << endl << vec_energy_acc[Nt] / vec_energy_acc[2512] << endl << endl;
	cout << "time = " << t2 - t1 << endl << endl;

	// Graph_Solution_in_two_planes_3dimen(cube, Nx, Ny, Nz, delta, dx, dy, dz, direction, "Bz PML T=8pi 4.csv");
	// Graph_Solution_in_one_planes_3dimen(cube, Nx, Ny, Nz, delta, dx, dy, dz, direction);

	//����� ������� ������� � ����

	ofstream numb_energy(file_energy);

	numb_energy << "Nx= " << Nx << ";  Ny= " << Ny << ";   Nz= " << Nz << "    " << endl;
	numb_energy << "delta_x= " << delta_x << ";  delta_y= " << delta_y << ";   delta_z= " << delta_z << "    " << endl;
	numb_energy << sigma_x << ";  " << sigma_y << ";   " << sigma_y << ";    " << endl<<endl;

	numb_energy << ";" << "R" << ";" << ";" << ";" << "R_accurate"  << endl;

	for (int s = 0; s < Nt + 1; s++)
	{
		numb_energy << dt * (double)(s) << ";" << vec_energy[s] <<";"<< ";" << dt * (double)(s) << ";" << vec_energy_acc[s] << endl;
	}
	numb_energy << endl << endl;
	numb_energy.close();

	//ofstream numb_data(file_data);
	//numb_data << "Nx= " << Nx << ";  Ny= " << Ny << ";   Nz= " << Nz << "    " << endl;
	//numb_data << "delta_x= " << delta_x << ";  delta_y= " << delta_y << ";   delta_z= " << delta_z << "    " << endl;
	//numb_data << sigma_x << ";  " << sigma_y << ";   " << sigma_y << ";    " << endl << endl;

	//numb_data << ";" << file_data << ";;;" << file_data << endl;

	//for (int s = 0; s < Nt + 1; s++)
	//{
	//	numb_data << dt * (double)(s) << ";" << vec_dataEy[s] << ";" << ";" << dt * (double)(s) << ";" << vec_dataBz[s] << endl;
	//}
	//numb_data << endl << endl;
	//numb_data.close();

	//ofstream numb_charact(file_charact);
	//numb_charact <<endl;
	//numb_charact << ";" << "maxData Electric" << ";;;;" << "maxData Magnetic" << ";;;" << "averageData Electric" << ";;;" << "averageDataMagnetic" << endl;
	//for (int s = 0; s < Nt + 1; s++)
	//{
	//	numb_charact << dt * (double)(s) << ";" << vec_max_dataE[s] << ";" << ";" << dt * (double)(s) << ";" << vec_max_dataB[s]<<";;;";
	//	numb_charact << dt * (double)(s) << ";" << vec_average_dataE[s] << ";" << ";" << dt * (double)(s) << ";" << vec_average_dataB[s];
	//	numb_charact << endl;

	//}
	//numb_charact.close();

	return result;
}
