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

template <class ftypePML>
ftypePML func_R(int n, ftypePML sigma_max, int delta, ftypePML step)
{
	return (ftypePML)exp(-(2.0 * (ftypePML)sigma_max * (ftypePML)delta * step)/ (ftypePML)(n + 1));
}

template <class ftype>
void Initializing_FDTD_3_dimen_Gauss_PML(vector< vector<vector<Component<ftype>>>>& cube,
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
					cube[i + index_start][j + index_start][k + index_start].Ez = exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse))
						* exp(-(x1 - x0) * (x1 - x0) / (ax_longitudinal * ax_longitudinal)) * exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse)) *
						sin(omega_ * t0_ - k_ * x1); // sin(phase), phase = omega * t - k * x
					cube[i + index_start][j + index_start][k + index_start].By = exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_longitudinal * ax_longitudinal)) * exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse)) *
						sin(omega_ * (t0_ + t) - k_ * x2); // sin(phase), phase = omega * t - k * x
					break;
				}
				case x_comp_yz: {
					// Ey and Bz has the same x coordinate
					cube[i + index_start][j + index_start][k + index_start].Ey = exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse))
						* exp(-(x1 - x0) * (x1 - x0) / (ax_longitudinal * ax_longitudinal)) * exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse)) *
						sin(omega_ * t0_ - k_ * x1); // sin(phase), phase = omega * t - k * x
					cube[i + index_start][j + index_start][k + index_start].Bz = exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_longitudinal * ax_longitudinal)) * exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse)) *
						sin(omega_ * (t0_ + t) - k_ * x2); // sin(phase), phase = omega * t - k * x
					break;
				}
				case y_comp_xz: {
					// Ex and Bz has the same y coordinate
					cube[i + index_start][j + index_start][k + index_start].Ex = exp(-(y1 - y0) * (y1 - y0) / (ay_longitudinal * ay_longitudinal))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse)) * exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse)) *
						sin(omega_ * t0_ - k_ * y1); // sin(phase), phase = omega * t - k * x
					cube[i + index_start][j + index_start][k + index_start].Bz = exp(-(y2 - y0) * (y2 - y0) / (ay_longitudinal * ay_longitudinal))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse)) * exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse)) *
						sin(omega_ * (t0_ + t) - k_ * y2); // sin(phase), phase = omega * t - k * x
					break;
				}
				case y_comp_zx: {
					// Ez and Bx has the same y coordinate
					cube[i + index_start][j + index_start][k + index_start].Ez = exp(-(y1 - y0) * (y1 - y0) / (ay_longitudinal * ay_longitudinal))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse)) * exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse)) *
						sin(omega_ * t0_ - k_ * y1); // sin(phase), phase = omega * t - k * x
					cube[i + index_start][j + index_start][k + index_start].Bx = exp(-(y2 - y0) * (y2 - y0) / (ay_longitudinal * ay_longitudinal))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse)) * exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse)) *
						sin(omega_ * (t0_ + t) - k_ * y2); // sin(phase), phase = omega * t - k * x
					break;
				}
				case z_comp_xy: {
					// Ex and By has the same z coordinate
					cube[i + index_start][j + index_start][k + index_start].Ex = -exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse)) * exp(-(z1 - z0) * (z1 - z0) / (az_longitudinal * az_longitudinal)) *
						sin(omega_ * t0_ + k_ * z1); // sin(phase), phase = omega * t - k * x
					cube[i + index_start][j + index_start][k + index_start].By = exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse)) * exp(-(z2 - z0) * (z2 - z0) / (az_longitudinal * az_longitudinal)) *
						sin(omega_ * (t0_ + t) + k_ * z2); // sin(phase), phase = omega * t - k * x
					break;
				}
				case z_comp_yx: {
					// Ey and Bx has the same z coordinate
					cube[i + index_start][j + index_start][k + index_start].Ey = exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse)) * exp(-(z1 - z0) * (z1 - z0) / (az_longitudinal * az_longitudinal)) *
						sin(omega_ * t0_ - k_ * z1); // sin(phase), phase = omega * t - k * x
					cube[i + index_start][j + index_start][k + index_start].Bx = exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse))
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
		cube[0][i][j].Ex = cube[Nx + 2 * delta][i][j].Ex;
		cube[0][i][j].Ey = cube[Nx + 2 * delta][i][j].Ey;
		cube[0][i][j].Ez = cube[Nx + 2 * delta][i][j].Ez;
		cube[0][i][j].Bx = cube[Nx + 2 * delta][i][j].Bx;
		cube[0][i][j].By = cube[Nx + 2 * delta][i][j].By;
		cube[0][i][j].Bz = cube[Nx + 2 * delta][i][j].Bz;
		cube[Nx + 2 * delta + 1][i][j].Ex = cube[1][i][j].Ex;
		cube[Nx + 2 * delta + 1][i][j].Ey = cube[1][i][j].Ey;
		cube[Nx + 2 * delta + 1][i][j].Ez = cube[1][i][j].Ez;
		cube[Nx + 2 * delta + 1][i][j].Bx = cube[1][i][j].Bx;
		cube[Nx + 2 * delta + 1][i][j].By = cube[1][i][j].By;
		cube[Nx + 2 * delta + 1][i][j].Bz = cube[1][i][j].Bz;
	}
	for (int i = 0; i < Nx + 2 * delta + 2; i++)
		for (int j = 0; j < Nz + 2 * delta + 2; j++) {
		cube[i][0][j].Ex = cube[i][Ny + 2 * delta][j].Ex;
		cube[i][0][j].Ey = cube[i][Ny + 2 * delta][j].Ey;
		cube[i][0][j].Ez = cube[i][Ny + 2 * delta][j].Ez;
		cube[i][0][j].Bx = cube[i][Ny + 2 * delta][j].Bx;
		cube[i][0][j].By = cube[i][Ny + 2 * delta][j].By;
		cube[i][0][j].Bz = cube[i][Ny + 2 * delta][j].Bz;
		cube[i][Ny + 2 * delta + 1][j].Ex = cube[i][1][j].Ex;
		cube[i][Ny + 2 * delta + 1][j].Ey = cube[i][1][j].Ey;
		cube[i][Ny + 2 * delta + 1][j].Ez = cube[i][1][j].Ez;
		cube[i][Ny + 2 * delta + 1][j].Bx = cube[i][1][j].Bx;
		cube[i][Ny + 2 * delta + 1][j].By = cube[i][1][j].By;
		cube[i][Ny + 2 * delta + 1][j].Bz = cube[i][1][j].Bz;
	}
	for (int i = 0; i < Nx + 2 * delta + 2; i++)
		for (int j = 0; j < Ny + 2 * delta + 2; j++) {
			cube[i][j][0].Ex = cube[i][j][Nz + 2 * delta].Ex;
			cube[i][j][0].Ey = cube[i][j][Nz + 2 * delta].Ey;
			cube[i][j][0].Ez = cube[i][j][Nz + 2 * delta].Ez;
			cube[i][j][0].Bx = cube[i][j][Nz + 2 * delta].Bx;
			cube[i][j][0].By = cube[i][j][Nz + 2 * delta].By;
			cube[i][j][0].Bz = cube[i][j][Nz + 2 * delta].Bz;
			cube[i][j][Nz + 2 * delta + 1].Ex = cube[i][j][1].Ex;
			cube[i][j][Nz + 2 * delta + 1].Ey = cube[i][j][1].Ey;
			cube[i][j][Nz + 2 * delta + 1].Ez = cube[i][j][1].Ez;
			cube[i][j][Nz + 2 * delta + 1].Bx = cube[i][j][1].Bx;
			cube[i][j][Nz + 2 * delta + 1].By = cube[i][j][1].By;
			cube[i][j][Nz + 2 * delta + 1].Bz = cube[i][j][1].Bz;
		}
}

template <class ftypePML>
void Initializing_cube_split_3_dimen_PML(vector< vector<vector<ComponentSplit<ftypePML>>>>& cube_split, int Nx, int Ny, int Nz, int delta)
{
	for (int i = 0; i < Nx + 2 * delta +2; i++)
		for (int j = 0; j < Ny + 2 * delta +2; j++)
			for (int k = 0; k < Nz + 2 * delta + 2; k++) {			
				cube_split[i][j][k].Exy = (ftypePML)0.0;
				cube_split[i][j][k].Exz = (ftypePML)0.0;
				cube_split[i][j][k].Eyx = (ftypePML)0.0;
				cube_split[i][j][k].Eyz = (ftypePML)0.0;
				cube_split[i][j][k].Ezx = (ftypePML)0.0;
				cube_split[i][j][k].Ezy = (ftypePML)0.0;

				cube_split[i][j][k].Bxy = (ftypePML)0.0;
				cube_split[i][j][k].Bxz = (ftypePML)0.0;
				cube_split[i][j][k].Byx = (ftypePML)0.0;
				cube_split[i][j][k].Byz = (ftypePML)0.0;
				cube_split[i][j][k].Bzx = (ftypePML)0.0;
				cube_split[i][j][k].Bzy = (ftypePML)0.0;
			}
}

template <class ftypePML>
void Initializing_Sigma_3_dimen_PML(vector< vector<vector<SIGMA<ftypePML>>>>& Sigma, int Nx, int Ny, int Nz, int delta, int n,
	ftypePML sigma_x, ftypePML sigma_y)
{
	ftypePML var_Sigma_max_x = sigma_x;
	ftypePML var_Sigma_max_y = sigma_y;
	ftypePML var_Sigma_max_z = sigma_y;

	cout << var_Sigma_max_x << "  " << var_Sigma_max_y << "   " << var_Sigma_max_z << "    " << endl;

	//cout << func_R(n, var_Sigma_max_x, delta, dx) << endl;

	for(int i = 1; i < Nx + 2 * delta + 1; i++)
		for (int j = 1; j < Ny + 2 * delta + 1; j++)
			for (int k = 1; k < Nz + 2 * delta + 1; k++) {
				Sigma[i][j][k].sigmaE_x = var_Sigma_max_x * pow(distanceE<ftypePML>(Nx, delta, i), n); 
				Sigma[i][j][k].sigmaH_x = var_Sigma_max_x * pow(distanceH<ftypePML>(Nx, delta, i), n);
	
				Sigma[i][j][k].sigmaE_y = var_Sigma_max_y * pow(distanceE<ftypePML>(Ny, delta, j), n);
				Sigma[i][j][k].sigmaH_y = var_Sigma_max_y * pow(distanceH<ftypePML>(Ny, delta, j), n);

				Sigma[i][j][k].sigmaE_z = var_Sigma_max_z * pow(distanceE<ftypePML>(Nz, delta, k), n);
				Sigma[i][j][k].sigmaH_z = var_Sigma_max_z * pow(distanceH<ftypePML>(Nz, delta, k), n);
			}
	//cout << sizeof(Sigma[10][18][11].sigmaH_x) << sizeof(double) << sizeof(float) << endl;
}

template <class ftypePML>
void Initializing_Coeff_3_dimen_PML(vector< vector<vector<COEFF<ftypePML>>>>& Coeff, vector< vector<vector<SIGMA<ftypePML>>>>& Sigma, int Nx, int Ny, int Nz, int delta, ftypePML dt)
{
	for (int i = 1; i < Nx + 2 * delta + 1; i++)
		for (int j = 1; j < Ny + 2 * delta + 1; j++)
			for (int k = 1; k < Nz + 2 * delta + 1; k++) {
				if ((i >= delta + 1) && (i < Nx + delta + 1) && (j >= delta + 1) && (j < Ny + delta + 1) && (k >= delta + 1) && (k < Nz + delta + 1))
				{
				}
				else {
					Coeff[i][j][k].Exy1 = (ftypePML)exp(-dt * Sigma[i][j][k].sigmaE_y);
					Coeff[i][j][k].Exz1 = (ftypePML)exp(-dt * Sigma[i][j][k].sigmaE_z);
					Coeff[i][j][k].Eyx1 = (ftypePML)exp(-dt * Sigma[i][j][k].sigmaE_x);
					Coeff[i][j][k].Eyz1 = (ftypePML)exp(-dt * Sigma[i][j][k].sigmaE_z);
					Coeff[i][j][k].Ezx1 = (ftypePML)exp(-dt * Sigma[i][j][k].sigmaE_x);
					Coeff[i][j][k].Ezy1 = (ftypePML)exp(-dt * Sigma[i][j][k].sigmaE_y);

					Coeff[i][j][k].Bxy1 = (ftypePML)exp(-dt * Sigma[i][j][k].sigmaH_y);
					Coeff[i][j][k].Bxz1 = (ftypePML)exp(-dt * Sigma[i][j][k].sigmaH_z);
					Coeff[i][j][k].Byx1 = (ftypePML)exp(-dt * Sigma[i][j][k].sigmaH_x);
					Coeff[i][j][k].Byz1 = (ftypePML)exp(-dt * Sigma[i][j][k].sigmaH_z);
					Coeff[i][j][k].Bzx1 = (ftypePML)exp(-dt * Sigma[i][j][k].sigmaH_x);
					Coeff[i][j][k].Bzy1 = (ftypePML)exp(-dt * Sigma[i][j][k].sigmaH_y);

					if (Sigma[i][j][k].sigmaE_x != (ftypePML)0.0) {
						Coeff[i][j][k].Eyx2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaE_x * ((ftypePML)1.0 - Coeff[i][j][k].Eyx1);
						Coeff[i][j][k].Ezx2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaE_x * ((ftypePML)1.0 - Coeff[i][j][k].Ezx1);
					}
					else {
						Coeff[i][j][k].Eyx2 = dt;
						Coeff[i][j][k].Ezx2 = dt;
					}
					if (Sigma[i][j][k].sigmaE_y != (ftypePML)0.0) {
						Coeff[i][j][k].Exy2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaE_y * ((ftypePML)1.0 - Coeff[i][j][k].Exy1);
						Coeff[i][j][k].Ezy2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaE_y * ((ftypePML)1.0 - Coeff[i][j][k].Ezy1);
					}
					else {
						Coeff[i][j][k].Exy2 = dt;
						Coeff[i][j][k].Ezy2 = dt;
					}
					if (Sigma[i][j][k].sigmaE_z != (ftypePML)0.0)
					{
						Coeff[i][j][k].Exz2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaE_z * ((ftypePML)1.0 - Coeff[i][j][k].Exz1);
						Coeff[i][j][k].Eyz2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaE_z * ((ftypePML)1.0 - Coeff[i][j][k].Eyz1);
					}
					else {
						Coeff[i][j][k].Exz2 = dt;
						Coeff[i][j][k].Eyz2 = dt;
					}
					if (Sigma[i][j][k].sigmaH_x != (ftypePML)0.0)
					{
						Coeff[i][j][k].Byx2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaH_x * ((ftypePML)1.0 - Coeff[i][j][k].Byx1);
						Coeff[i][j][k].Bzx2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaH_x * ((ftypePML)1.0 - Coeff[i][j][k].Bzx1);
					}
					else {
						Coeff[i][j][k].Byx2 = dt;
						Coeff[i][j][k].Bzx2 = dt;
					}
					if (Sigma[i][j][k].sigmaH_y != (ftypePML)0.0)
					{
						Coeff[i][j][k].Bxy2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaH_y * ((ftypePML)1.0 - Coeff[i][j][k].Bxy1);
						Coeff[i][j][k].Bzy2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaH_y * ((ftypePML)1.0 - Coeff[i][j][k].Bzy1);
					}
					else {
						Coeff[i][j][k].Bxy2 = dt;
						Coeff[i][j][k].Bzy2 = dt;
					}
					if (Sigma[i][j][k].sigmaH_z != (ftypePML)0.0)
					{
						Coeff[i][j][k].Bxz2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaH_z * ((ftypePML)1.0 - Coeff[i][j][k].Bxz1);
						Coeff[i][j][k].Byz2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaH_z * ((ftypePML)1.0 - Coeff[i][j][k].Byz1);
					}
					else {
						Coeff[i][j][k].Bxz2 = dt;
						Coeff[i][j][k].Byz2 = dt;
					}
				}
		}
	// cout << sizeof(Coeff[1][5][4].Bxy2) << sizeof(Coeff[2][4][3].Exz1) << endl;

}

template <class ftype, class ftypePML>
void Update_electric_field_three_dimen_PML(vector< vector<vector<Component<ftype>>>>& cube,
	vector< vector<vector<ComponentSplit<ftypePML>>>>& cube_split, vector< vector<vector<COEFF<ftypePML>>>>& Coeff,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k)
{
	ftypePML tExy, tExz, tEyx, tEyz, tEzx, tEzy;

	tExy = cube_split[i][j][k].Exy * Coeff[i][j][k].Exy1 + ((ftypePML)cube[i][j + 1][k].Bz - (ftypePML)cube[i][j][k].Bz) * Coeff[i][j][k].Exy2 * (_1dy);

	tExz = cube_split[i][j][k].Exz * Coeff[i][j][k].Exz1 - ((ftypePML)cube[i][j][k + 1].By - (ftypePML)cube[i][j][k].By) * Coeff[i][j][k].Exz2 * (_1dz);

	tEyx = cube_split[i][j][k].Eyx * Coeff[i][j][k].Eyx1 - ((ftypePML)cube[i + 1][j][k].Bz - (ftypePML)cube[i][j][k].Bz) * Coeff[i][j][k].Eyx2 * (_1dx);

	tEyz = cube_split[i][j][k].Eyz * Coeff[i][j][k].Eyz1 + ((ftypePML)cube[i][j][k + 1].Bx - (ftypePML)cube[i][j][k].Bx) * Coeff[i][j][k].Eyz2 * (_1dz);

	tEzx = cube_split[i][j][k].Ezx * Coeff[i][j][k].Ezx1 + ((ftypePML)cube[i + 1][j][k].By - (ftypePML)cube[i][j][k].By) * Coeff[i][j][k].Ezx2 * (_1dx);

	tEzy = cube_split[i][j][k].Ezy * Coeff[i][j][k].Ezy1 - ((ftypePML)cube[i][j + 1][k].Bx - (ftypePML)cube[i][j][k].Bx) * Coeff[i][j][k].Ezy2 * (_1dy);

	cube_split[i][j][k].Exy = tExy;
	cube_split[i][j][k].Exz = tExz;
	cube_split[i][j][k].Eyx = tEyx;
	cube_split[i][j][k].Eyz = tEyz;
	cube_split[i][j][k].Ezx = tEzx;
	cube_split[i][j][k].Ezy = tEzy;

	cube[i][j][k].Ex = (ftype)(tExy + tExz);
	cube[i][j][k].Ey = (ftype)(tEyx + tEyz);
	cube[i][j][k].Ez = (ftype)(tEzx + tEzy);
}

template <class ftype, class ftypePML>
void Update_magnetic_field_three_dimen_PML(vector< vector<vector<Component<ftype>>>>& cube,
	vector< vector<vector<ComponentSplit<ftypePML>>>>& cube_split, vector< vector<vector<COEFF<ftypePML>>>>& Coeff,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k)
{
	ftypePML tBxy, tBxz, tByx, tByz, tBzx, tBzy;

	tBxy = cube_split[i][j][k].Bxy * Coeff[i][j][k].Bxy1 - ((ftypePML)cube[i][j][k].Ez - (ftypePML)cube[i][j - 1][k].Ez) * Coeff[i][j][k].Bxy2 * (_1dy);

	tBxz = cube_split[i][j][k].Bxz * Coeff[i][j][k].Bxz1 + ((ftypePML)cube[i][j][k].Ey - (ftypePML)cube[i][j][k - 1].Ey) * Coeff[i][j][k].Bxz2 * (_1dz);

	tByx = cube_split[i][j][k].Byx * Coeff[i][j][k].Byx1 + ((ftypePML)cube[i][j][k].Ez - (ftypePML)cube[i - 1][j][k].Ez) * Coeff[i][j][k].Byx2 * (_1dx);

	tByz = cube_split[i][j][k].Byz * Coeff[i][j][k].Byz1 - ((ftypePML)cube[i][j][k].Ex - (ftypePML)cube[i][j][k - 1].Ex) * Coeff[i][j][k].Byz2 * (_1dz);

	tBzx = cube_split[i][j][k].Bzx * Coeff[i][j][k].Bzx1 - ((ftypePML)cube[i][j][k].Ey - (ftypePML)cube[i - 1][j][k].Ey) * Coeff[i][j][k].Bzx2 * (_1dx);

	tBzy = cube_split[i][j][k].Bzy * Coeff[i][j][k].Bzy1 + ((ftypePML)cube[i][j][k].Ex - (ftypePML)cube[i][j - 1][k].Ex) * Coeff[i][j][k].Bzy2 * (_1dy);
	
	cube_split[i][j][k].Bxy = tBxy;
	cube_split[i][j][k].Bxz = tBxz;
	cube_split[i][j][k].Byx = tByx;
	cube_split[i][j][k].Byz = tByz;
	cube_split[i][j][k].Bzx = tBzx;
	cube_split[i][j][k].Bzy = tBzy;

	cube[i][j][k].Bx = (ftype)(tBxy + tBxz);
	cube[i][j][k].By = (ftype)(tByx + tByz);
	cube[i][j][k].Bz = (ftype)(tBzx + tBzy);
}

template <class ftype>
void Add_Currents_3_dimen(vector< vector<vector<Component<ftype>>>>& cube,
	int Nx, int Ny, int Nz, int delta, ftype dx, ftype dy, ftype dz, const ftype dt, int it, int index_start, Direction direction,
	pair <ftype, ftype> ab, pair <ftype, ftype> cd, pair <ftype, ftype> fg)
{
	ftype x1, x2, y1, y2, z1, z2, t1, t2;

	ftype ax_transverse = ab.second * 3. / 4.; // поперечное
	ftype ay_transverse = cd.second * 3. / 4.;
	ftype az_transverse = fg.second * 3. / 4.;

	ftype tp_x = ab.second / 6.;
	ftype tp_y = cd.second / 6.;
	ftype tp_z = fg.second / 6.;

	ftype lambda_ = 2. * M_PI / 4.; // wvelength               //?
	ftype k_ = 2. * M_PI / lambda_; // wavenumber = 2 pi/ wavelength
	ftype omega_ = 2. * M_PI / lambda_; // 2 pi c /wavelength 
	ftype A = 1.; //amplitude

	t1 = dt * (ftype)it;
	t2 = t1 + (ftype)0.5 * dt;

	ftype x0 = (ab.second - ab.first) / 2.;
	ftype y0 = (cd.second - cd.first) / 2.;
	ftype z0 = (fg.second - fg.first) / 2.;
	ftype t0_x = 3. * tp_x;
	ftype t0_y = 3. * tp_y;
	ftype t0_z = 3. * tp_z;

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
			cube[offset + index_start][j + index_start][k + index_start].Ey += A * exp(-(t1 - t0_x) * (t1 - t0_x) / (tp_x * tp_x))
				* exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse)) *
				sin(omega_ * t1 - k_ * x1); // sin(phase), phase = omega * t - k * x
			cube[offset + index_start][j + index_start][k + index_start].Bz += A * exp(-(t2 - t0_x) * (t2 - t0_x) / (tp_x * tp_x))
				* exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse)) *
				sin(omega_ * t2 - k_ * x2); // sin(phase), phase = omega * t - k * x
		}

}

template <class ftype, class ftypePML>
double FDTD_three_dimen_with_PML(int Nx, int Ny, int Nz, ftype T, ftype dt, int delta, ftypePML sigma_x, ftypePML sigma_y, string type_sum, string name_file)
{
	setlocale(LC_ALL, "Russian");

	int Nt = T / dt; int n = 4;
	// delta  =  width boundary layer

	vector<vector<vector<Component<ftype>>>> cube(Nx + 2 * delta + 2, vector<vector<Component<ftype>>>(Ny + 2 * delta + 2, vector<Component<ftype>>(Nz + 2 * delta + 2)));
	vector< vector<vector<ComponentSplit<ftypePML>>>> cube_split(Nx + 2 * delta + 2, vector< vector<ComponentSplit<ftypePML>>>(Ny + 2 * delta + 2, vector<ComponentSplit<ftypePML>>(Nz + 2 * delta + 2)));
	vector< vector<vector<SIGMA<ftypePML>>>> Sigma(Nx + 2 * delta + 2, vector< vector<SIGMA<ftypePML>>>(Ny + 2 * delta + 2, vector<SIGMA<ftypePML>>(Nz + 2 * delta + 2)));
	vector< vector<vector<COEFF<ftypePML>>>> Coeff(Nx + 2 * delta + 2, vector< vector<COEFF<ftypePML>>>(Ny + 2 * delta + 2, vector<COEFF<ftypePML>>(Nz + 2 * delta + 2)));

	pair<ftype, ftype> ab(0.0, 4.0 * M_PI), cd(0.0, 16.0 * M_PI), fg(0.0, 16.0 * M_PI);
	ftype dx = (ab.second - ab.first) / (ftype)Nx, dy = (cd.second - cd.first) / (ftype)Ny, dz = (fg.second - fg.first) / (ftype)Nz;

	ftype dt_x = dt / (dx), dt_y = dt / (dy), dt_z = dt / (dz);
	ftypePML _1dx = (ftypePML)1.0 / dx, _1dy = (ftypePML)1.0 / dy, _1dz = (ftypePML)1.0 / dz;
	vector<double> vec_energy(Nt + 1);

	Direction direction = x_comp_yz;

	// Initializing_FDTD_3_dimen_Gauss_PML<ftype>(cube, Nx, Ny, Nz, delta, dx, dy, dz, dt, delta + 1, direction, ab, cd, fg);
	Initializing_cube_split_3_dimen_PML<ftypePML>(cube_split, Nx, Ny, Nz, delta);
	Initializing_Sigma_3_dimen_PML<ftypePML>(Sigma, Nx, Ny, Nz, delta, n, sigma_x, sigma_y);
	Initializing_Coeff_3_dimen_PML<ftypePML>(Coeff, Sigma, Nx, Ny, Nz, delta, dt);

	// Graph_Solution_in_two_planes_3dimen(cube, Nx, Ny, Nz, delta, dx, dy, dz, direction, "Bz PML T=8pi 0.csv");

	if (type_sum == "Kahan")
	{
	}
	else
	{

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
			
			 Add_Currents_3_dimen<ftype>(cube, Nx, Ny, Nz, delta, dx, dy, dz, dt, it, delta + 1, direction, ab, cd, fg);
			 vec_energy[0] = CalculateEnergyPML_3dimen(cube, Nx, Ny, Nz, delta);

#pragma omp parallel for collapse(3)
			for (int i = 1; i < Nx + 2 * delta + 1; i++)
				for (int j = 1; j < Ny + 2 * delta + 1; j++)
					for (int k = 1; k < Nz + 2 * delta + 1; k++)
					{
						if ((i >= delta + 1) && (i < Nx + delta + 1)&& (j >= delta + 1) && (j < Ny + delta + 1) && (k >= delta + 1) && (k < Nz + delta + 1))
						{
							Update_electric_field_three_dimen<ftype>(cube, dt_x, dt_y, dt_z, i, j, k);
						}
						else {
							Update_electric_field_three_dimen_PML<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
						}
					}

#pragma omp parallel for collapse(3)
			for (int i = 1; i < Nx + 2 * delta + 1; i++)
				for (int j = 1; j < Ny + 2 * delta + 1; j++)
					for (int k = 1; k < Nz + 2 * delta + 1; k++)
					{
						if ((i >= delta + 1) && (i < Nx + delta + 1) && (j >= delta + 1) && (j < Ny + delta + 1) && (k >= delta + 1) && (k < Nz + delta + 1))
						{
							Update_magnetic_field_three_dimen<ftype>(cube, dt_x, dt_y, dt_z, i, j, k);
						}
						else {
							Update_magnetic_field_three_dimen_PML<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
						}
					}
			vec_energy[it + 1] = CalculateEnergyPML_3dimen(cube, Nx, Ny, Nz, delta);
		}
	}

	double result = vec_energy[Nt] / vec_energy[2512];

	cout << endl << result << endl;
	// Graph_Solution_in_two_planes_3dimen(cube, Nx, Ny, Nz, delta, dx, dy, dz, direction, "Bz PML T=8pi 4.csv");

	//����� ������� ������� � ����
	//ofstream numb_energy(name_file);
	//for (int s = 0; s < Nt; s++)
	//{
	//	numb_energy << dt * (double)(s + 1) << ";" << vec_energy[s] << endl;
	//}
	//numb_energy << endl << endl;
	//numb_energy.close();

	return result;
}

template <class ftype>
double CalculateEnergyPML_3dimen(vector<vector<vector<Component<ftype>>>>& cube, int Nx, int Ny, int Nz, int delta)
{
	double energy = 0.0;

	for (int i =  1; i < Nx +  2 * delta + 1; i++)
		for (int j =  1; j < Ny + 2 * delta + 1; j++)
			for (int k = 1; k < Nz + 2 * delta + 1; k++)
			{
				energy += (double)cube[i][j][k].Ex * (double)cube[i][j][k].Ex
					+ (double)cube[i][j][k].Ey * (double)cube[i][j][k].Ey
					+ (double)cube[i][j][k].Ez * (double)cube[i][j][k].Ez;
				energy += (double)cube[i][j][k].Bx * (double)cube[i][j][k].Bx
					+ (double)cube[i][j][k].By * (double)cube[i][j][k].By
					+ (double)cube[i][j][k].Bz * (double)cube[i][j][k].Bz;
			}
	return energy;
}

template <class ftype>
void Graph_Solution_in_two_planes_3dimen(vector< vector<vector<Component<ftype>>>>& cube, int Nx, int Ny, int Nz, int delta, ftype dx, ftype dy, ftype dz, Direction direction, string file_name)
{
	ofstream numbEx("graph_solution_E_x.csv"), numbBx("graph_solution_B_x.csv");
	ofstream numbEy("graph_solution_E_y.csv"), numbBy("graph_solution_B_y.csv");
	ofstream numbEz("graph_solution_E_z.csv"), numbBz(file_name);

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
	ftype z0 = 8.0 * M_PI;

	for (int i = 1; i < Ny + 2 * delta + 1; i++)
	{
		ftype y = dy * (ftype)i;
		numbEx << ";" << y << ";"; 			numbBx << ";" << y << ";";			numbEy << ";" << y << ";"; 			numbBy << ";" << y << ";";
		numbEz << ";" << y << ";"; 			numbBz << ";" << y << ";";

		for (int j = 1; j < Nx + 2 * delta + 1; j++)
		{
			numbEx << cube[j][i][Nz / 2 + delta].Ex << ";";
			numbBx << cube[j][i][Nz / 2 + delta].Bx << ";";
			numbEy << cube[j][i][Nz / 2 + delta].Ey << ";";
			numbBy << cube[j][i][Nz / 2 + delta].By << ";";
			numbEz << cube[j][i][Nz / 2 + delta].Ez << ";";
			numbBz << cube[j][i][Nz / 2 + delta].Bz << ";";
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