#define _USE_MATH_DEFINES
#include <algorithm>
#include <iostream>
#include <iomanip>
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

template <class ftypePML>
void Graph_for_Sigma_three_dimen(vector<vector<vector<SIGMA<double>>>>& Sigma, int Nx, int Ny, int Nz, int delta, ftypePML dx, string file_sigma)
{
	//����� ������� ���� � ����
	ofstream numb_sigma(file_sigma);

	numb_sigma << ";"<< "value sigma_x" << endl;

	for (int s = 1; s < Nx + 2 * delta + 1; s++)
	{
		numb_sigma << dx * (ftypePML)(s - 1) << ";" << Sigma[s][delta/2][delta/2].sigmaH_x << endl;
		numb_sigma << dx * ((ftypePML)(s - 1) + 0.5) << ";" << Sigma[s][delta / 2][delta / 2].sigmaE_x << endl;
	}
	numb_sigma << endl << endl;

	numb_sigma.close();
}

template <class ftypePML>
void Graph_for_Coeff_three_dimen(vector<vector<vector<COEFF<ftypePML>>>>& Coeff, int Nx, int Ny, int Nz, int delta, ftypePML dx, string file_coeff)
{
	//����� ������� ���� � ����
	ofstream numb_sigma(file_coeff);

	numb_sigma << ";" << "Coeff1" <<";;;"<<"Coeff2"<< endl;

	for (int s = 1; s < Nx + 2 * delta + 1; s++)
	{
		numb_sigma << dx * (ftypePML)(s - 1) << ";" << std::setprecision(16)<< Coeff[s][delta / 2][delta / 2].Ezx1 << ";;";
		numb_sigma << dx * (ftypePML)(s - 1) << ";" << std::setprecision(16) << Coeff[s][delta / 2][delta / 2].Ezx2 << endl;

		//numb_sigma << dx * ((ftypePML)(s - 1) + 0.5) << ";" << Sigma[s][delta / 2][delta / 2].sigmaE_x << endl;
	}
	numb_sigma << endl << endl;

	numb_sigma.close();
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
				Sigma[i][j][k].sigmaE_x = var_Sigma_max_x * powf(distanceE<ftypePML>(Nx, delta, i), n);
				Sigma[i][j][k].sigmaH_x = var_Sigma_max_x * powf(distanceH<ftypePML>(Nx, delta, i), n);
	
				Sigma[i][j][k].sigmaE_y = var_Sigma_max_y * powf(distanceE<ftypePML>(Ny, delta, j), n);
				Sigma[i][j][k].sigmaH_y = var_Sigma_max_y * powf(distanceH<ftypePML>(Ny, delta, j), n);

				Sigma[i][j][k].sigmaE_z = var_Sigma_max_z * powf(distanceE<ftypePML>(Nz, delta, k), n);
				Sigma[i][j][k].sigmaH_z = var_Sigma_max_z * powf(distanceH<ftypePML>(Nz, delta, k), n);
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
						Coeff[i][j][k].Eyx2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaE_x - Coeff[i][j][k].Eyx1 / Sigma[i][j][k].sigmaE_x;
						Coeff[i][j][k].Ezx2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaE_x - Coeff[i][j][k].Ezx1 / Sigma[i][j][k].sigmaE_x;
					}
					else {
						Coeff[i][j][k].Eyx2 = dt;
						Coeff[i][j][k].Ezx2 = dt;
					}
					if (Sigma[i][j][k].sigmaE_y != (ftypePML)0.0) {
						Coeff[i][j][k].Exy2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaE_y - Coeff[i][j][k].Exy1 / Sigma[i][j][k].sigmaE_y;
						Coeff[i][j][k].Ezy2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaE_y - Coeff[i][j][k].Ezy1 / Sigma[i][j][k].sigmaE_y;
					}
					else {
						Coeff[i][j][k].Exy2 = dt;
						Coeff[i][j][k].Ezy2 = dt;
					}
					if (Sigma[i][j][k].sigmaE_z != (ftypePML)0.0)
					{
						Coeff[i][j][k].Exz2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaE_z - Coeff[i][j][k].Exz1 / Sigma[i][j][k].sigmaE_z;
						Coeff[i][j][k].Eyz2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaE_z - Coeff[i][j][k].Eyz1 / Sigma[i][j][k].sigmaE_z;
					}
					else {
						Coeff[i][j][k].Exz2 = dt;
						Coeff[i][j][k].Eyz2 = dt;
					}
					if (Sigma[i][j][k].sigmaH_x != (ftypePML)0.0)
					{
						Coeff[i][j][k].Byx2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaH_x - Coeff[i][j][k].Byx1 / Sigma[i][j][k].sigmaH_x;
						Coeff[i][j][k].Bzx2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaH_x - Coeff[i][j][k].Bzx1 / Sigma[i][j][k].sigmaH_x;
					}
					else {
						Coeff[i][j][k].Byx2 = dt;
						Coeff[i][j][k].Bzx2 = dt;
					}
					if (Sigma[i][j][k].sigmaH_y != (ftypePML)0.0)
					{
						Coeff[i][j][k].Bxy2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaH_y - Coeff[i][j][k].Bxy1 / Sigma[i][j][k].sigmaH_y;
						Coeff[i][j][k].Bzy2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaH_y - Coeff[i][j][k].Bzy1 / Sigma[i][j][k].sigmaH_y;
					}
					else {
						Coeff[i][j][k].Bxy2 = dt;
						Coeff[i][j][k].Bzy2 = dt;
					}
					if (Sigma[i][j][k].sigmaH_z != (ftypePML)0.0)
					{
						Coeff[i][j][k].Bxz2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaH_z - Coeff[i][j][k].Bxz1 / Sigma[i][j][k].sigmaH_z;
						Coeff[i][j][k].Byz2 = (ftypePML)1.0 / Sigma[i][j][k].sigmaH_z - Coeff[i][j][k].Byz1 / Sigma[i][j][k].sigmaH_z;
					}
					else {
						Coeff[i][j][k].Bxz2 = dt;
						Coeff[i][j][k].Byz2 = dt;
					}
				}
		}
	// cout << sizeof(Coeff[1][5][4].Bxy2) << sizeof(Coeff[2][4][3].Exz1) << endl;

}

template <class ftypePML>
void Initializing_Coeff_3_dimen_PML_correct(vector< vector<vector<COEFF<ftypePML>>>>& Coeff, vector< vector<vector<SIGMA<double>>>>& Sigma, int Nx, int Ny, int Nz, int delta, ftypePML dt)
{
	double Exy1, Exz1, Ezx1, Ezy1, Eyx1, Eyz1;
	double Bxy1, Bxz1, Bzx1, Bzy1, Byx1, Byz1;

	for (int i = 1; i < Nx + 2 * delta + 1; i++)
		for (int j = 1; j < Ny + 2 * delta + 1; j++)
			for (int k = 1; k < Nz + 2 * delta + 1; k++) {
				if ((i >= delta + 1) && (i < Nx + delta + 1) && (j >= delta + 1) && (j < Ny + delta + 1) && (k >= delta + 1) && (k < Nz + delta + 1))
				{
				}
				else {

					Exy1 = exp(-dt * Sigma[i][j][k].sigmaE_y);
					Exz1 = exp(-dt * Sigma[i][j][k].sigmaE_z);
					Eyx1 = exp(-dt * Sigma[i][j][k].sigmaE_x);
					Eyz1 = exp(-dt * Sigma[i][j][k].sigmaE_z);
					Ezx1 = exp(-dt * Sigma[i][j][k].sigmaE_x);
					Ezy1 = exp(-dt * Sigma[i][j][k].sigmaE_y);

					Bxy1 = exp(-dt * Sigma[i][j][k].sigmaH_y);
					Bxz1 = exp(-dt * Sigma[i][j][k].sigmaH_z);
					Byx1 = exp(-dt * Sigma[i][j][k].sigmaH_x);
					Byz1 = exp(-dt * Sigma[i][j][k].sigmaH_z);
					Bzx1 = exp(-dt * Sigma[i][j][k].sigmaH_x);
					Bzy1 = exp(-dt * Sigma[i][j][k].sigmaH_y);

					Coeff[i][j][k].Exy1 = (ftypePML)Exy1;
					Coeff[i][j][k].Exz1 = (ftypePML)Exz1;
					Coeff[i][j][k].Eyx1 = (ftypePML)Eyx1;
					Coeff[i][j][k].Eyz1 = (ftypePML)Eyz1;
					Coeff[i][j][k].Ezx1 = (ftypePML)Ezx1;
					Coeff[i][j][k].Ezy1 = (ftypePML)Ezy1;

					Coeff[i][j][k].Bxy1 = (ftypePML)Bxy1;
					Coeff[i][j][k].Bxz1 = (ftypePML)Bxz1;
					Coeff[i][j][k].Byx1 = (ftypePML)Byx1;
					Coeff[i][j][k].Byz1 = (ftypePML)Byz1;
					Coeff[i][j][k].Bzx1 = (ftypePML)Bzx1;
					Coeff[i][j][k].Bzy1 = (ftypePML)Bzy1;

					if (Sigma[i][j][k].sigmaE_x != (ftypePML)0.0) {
						Coeff[i][j][k].Eyx2 = 1.0 / Sigma[i][j][k].sigmaE_x - Eyx1 / Sigma[i][j][k].sigmaE_x;
						Coeff[i][j][k].Ezx2 = 1.0 / Sigma[i][j][k].sigmaE_x - Ezx1 / Sigma[i][j][k].sigmaE_x;
					}
					else {
						Coeff[i][j][k].Eyx2 = dt;
						Coeff[i][j][k].Ezx2 = dt;
					}
					if (Sigma[i][j][k].sigmaE_y != (ftypePML)0.0) {
						Coeff[i][j][k].Exy2 = 1.0 / Sigma[i][j][k].sigmaE_y - Exy1 / Sigma[i][j][k].sigmaE_y;
						Coeff[i][j][k].Ezy2 = 1.0 / Sigma[i][j][k].sigmaE_y - Ezy1 / Sigma[i][j][k].sigmaE_y;
					}
					else {
						Coeff[i][j][k].Exy2 = dt;
						Coeff[i][j][k].Ezy2 = dt;
					}
					if (Sigma[i][j][k].sigmaE_z != (ftypePML)0.0)
					{
						Coeff[i][j][k].Exz2 = 1.0 / Sigma[i][j][k].sigmaE_z - Exz1 / Sigma[i][j][k].sigmaE_z;
						Coeff[i][j][k].Eyz2 = 1.0 / Sigma[i][j][k].sigmaE_z - Eyz1 / Sigma[i][j][k].sigmaE_z;
					}
					else {
						Coeff[i][j][k].Exz2 = dt;
						Coeff[i][j][k].Eyz2 = dt;
					}
					if (Sigma[i][j][k].sigmaH_x != (ftypePML)0.0)
					{
						Coeff[i][j][k].Byx2 = 1.0 / Sigma[i][j][k].sigmaH_x - Byx1 / Sigma[i][j][k].sigmaH_x;
						Coeff[i][j][k].Bzx2 = 1.0 / Sigma[i][j][k].sigmaH_x - Bzx1 / Sigma[i][j][k].sigmaH_x;
					}
					else {
						Coeff[i][j][k].Byx2 = dt;
						Coeff[i][j][k].Bzx2 = dt;
					}
					if (Sigma[i][j][k].sigmaH_y != (ftypePML)0.0)
					{
						Coeff[i][j][k].Bxy2 = 1.0 / Sigma[i][j][k].sigmaH_y - Bxy1 / Sigma[i][j][k].sigmaH_y;
						Coeff[i][j][k].Bzy2 = 1.0 / Sigma[i][j][k].sigmaH_y - Bzy1 / Sigma[i][j][k].sigmaH_y;
					}
					else {
						Coeff[i][j][k].Bxy2 = dt;
						Coeff[i][j][k].Bzy2 = dt;
					}
					if (Sigma[i][j][k].sigmaH_z != (ftypePML)0.0)
					{
						Coeff[i][j][k].Bxz2 = 1.0 / Sigma[i][j][k].sigmaH_z - Bxz1 / Sigma[i][j][k].sigmaH_z;
						Coeff[i][j][k].Byz2 = 1.0 / Sigma[i][j][k].sigmaH_z - Byz1 / Sigma[i][j][k].sigmaH_z;
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

template <class ftype, class ftypePML>
void Update_electric_field_three_dimen_PML_Kahan(vector< vector<vector<Component<ftype>>>>& cube,
	vector< vector<vector<ComponentSplit<ftypePML>>>>& cube_split, vector< vector<vector<COEFF<ftypePML>>>>& Coeff,
	vector<vector<vector<ComponentSplit<ftypePML>>>>& compensatorPML, vector<vector<vector<ComponentSplit<ftypePML>>>>& compensatorPML2,
	vector<vector<vector<Component<ftype>>>>& compensator, vector<vector<vector<Component<ftype>>>>& compensator2,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k)
{
	ftypePML tExy, tExz, tEyx, tEyz, tEzx, tEzy;
	ftype sum, y;
	//
	sum = cube_split[i][j][k].Exy * Coeff[i][j][k].Exy1;
	y = ((ftypePML)cube[i][j + 1][k].Bz - (ftypePML)cube[i][j][k].Bz
			- ((ftypePML)compensator[i][j + 1][k].Bz - (ftypePML)compensator[i][j][k].Bz)) * Coeff[i][j][k].Exy2 * (_1dy)
				- compensatorPML[i][j][k].Exy * Coeff[i][j][k].Exy1;
	tExy = sum + y;
	compensatorPML2[i][j][k].Exy = (tExy - sum) - y;
	//
	sum = cube_split[i][j][k].Exz * Coeff[i][j][k].Exz1;
	y = - ((ftypePML)cube[i][j][k + 1].By - (ftypePML)cube[i][j][k].By
			- ((ftypePML)compensator[i][j][k + 1].By - (ftypePML)compensator[i][j][k].By)) * Coeff[i][j][k].Exz2 * (_1dz)
				- compensatorPML[i][j][k].Exz * Coeff[i][j][k].Exz1;
	tExz = sum + y;
	compensatorPML2[i][j][k].Exz = (tExz - sum) - y;
	//
	sum = cube_split[i][j][k].Eyx * Coeff[i][j][k].Eyx1;
	y = - ((ftypePML)cube[i + 1][j][k].Bz - (ftypePML)cube[i][j][k].Bz
			- ((ftypePML)compensator[i + 1][j][k].Bz - (ftypePML)compensator[i][j][k].Bz)) * Coeff[i][j][k].Eyx2 * (_1dx)
				- compensatorPML[i][j][k].Eyx * Coeff[i][j][k].Eyx1;
	tEyx = sum + y;
	compensatorPML2[i][j][k].Eyx = (tEyx - sum) - y;
	//
	sum = cube_split[i][j][k].Eyz * Coeff[i][j][k].Eyz1;
	y = ((ftypePML)cube[i][j][k + 1].Bx - (ftypePML)cube[i][j][k].Bx
			- ((ftypePML)compensator[i][j][k + 1].Bx - (ftypePML)compensator[i][j][k].Bx)) * Coeff[i][j][k].Eyz2 * (_1dz)
				- compensatorPML[i][j][k].Eyz * Coeff[i][j][k].Eyz1;
	tEyz = sum + y;
	compensatorPML2[i][j][k].Eyz = (tEyz - sum) - y;
	//
	sum = cube_split[i][j][k].Ezx * Coeff[i][j][k].Ezx1;
	y = ((ftypePML)cube[i + 1][j][k].By - (ftypePML)cube[i][j][k].By
			- ((ftypePML)compensator[i + 1][j][k].By - (ftypePML)compensator[i][j][k].By)) * Coeff[i][j][k].Ezx2 * (_1dx)
				- compensatorPML[i][j][k].Ezx * Coeff[i][j][k].Ezx1;
	tEzx = sum + y;
	compensatorPML2[i][j][k].Ezx = (tEzx - sum) - y;
	//
	sum = cube_split[i][j][k].Ezy * Coeff[i][j][k].Ezy1;
	y = - ((ftypePML)cube[i][j + 1][k].Bx - (ftypePML)cube[i][j][k].Bx
			- ((ftypePML)compensator[i][j + 1][k].Bx - (ftypePML)compensator[i][j][k].Bx)) * Coeff[i][j][k].Ezy2 * (_1dy)
				- compensatorPML[i][j][k].Ezy * Coeff[i][j][k].Ezy1;
	tEzy = sum + y;
	compensatorPML2[i][j][k].Ezy = (tEzy - sum) - y;
	//
	cube_split[i][j][k].Exy = tExy;
	cube_split[i][j][k].Exz = tExz;
	cube_split[i][j][k].Eyx = tEyx;
	cube_split[i][j][k].Eyz = tEyz;
	cube_split[i][j][k].Ezx = tEzx;
	cube_split[i][j][k].Ezy = tEzy;

	cube[i][j][k].Ex = (ftype)(tExy + tExz);
	cube[i][j][k].Ey = (ftype)(tEyx + tEyz);
	cube[i][j][k].Ez = (ftype)(tEzx + tEzy);

	compensator2[i][j][k].Ex = (ftype)(((ftypePML)cube[i][j][k].Ex - tExy) - tExz);
	compensator2[i][j][k].Ey = (ftype)(((ftypePML)cube[i][j][k].Ey - tEyx) - tEyz);
	compensator2[i][j][k].Ez = (ftype)(((ftypePML)cube[i][j][k].Ez - tEzx) - tEzy);
}

template <class ftype, class ftypePML>
void Update_magnetic_field_three_dimen_PML_Kahan(vector< vector<vector<Component<ftype>>>>& cube,
	vector< vector<vector<ComponentSplit<ftypePML>>>>& cube_split, vector< vector<vector<COEFF<ftypePML>>>>& Coeff,
	vector<vector<vector<ComponentSplit<ftypePML>>>>& compensatorPML, vector<vector<vector<ComponentSplit<ftypePML>>>>& compensatorPML2,
	vector<vector<vector<Component<ftype>>>>& compensator, vector<vector<vector<Component<ftype>>>>& compensator2,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k)
{
	ftypePML tBxy, tBxz, tByx, tByz, tBzx, tBzy;
	ftype sum, y;
	//
	sum = cube_split[i][j][k].Bxy * Coeff[i][j][k].Bxy1;
	y = - ((ftypePML)cube[i][j][k].Ez - (ftypePML)cube[i][j - 1][k].Ez
			- ((ftypePML)compensator[i][j][k].Ez - (ftypePML)compensator[i][j - 1][k].Ez))* Coeff[i][j][k].Bxy2 * (_1dy)
				- compensatorPML[i][j][k].Bxy * Coeff[i][j][k].Bxy1;
	tBxy = sum + y;
	compensatorPML2[i][j][k].Bxy = (tBxy - sum) - y;
	//
	sum = cube_split[i][j][k].Bxz * Coeff[i][j][k].Bxz1;
	y = ((ftypePML)cube[i][j][k].Ey - (ftypePML)cube[i][j][k - 1].Ey
		- ((ftypePML)compensator[i][j][k].Ey - (ftypePML)compensator[i][j][k - 1].Ey)) * Coeff[i][j][k].Bxz2 * (_1dz)
				- compensatorPML[i][j][k].Bxz * Coeff[i][j][k].Bxz1;
	tBxz = sum + y;
	compensatorPML2[i][j][k].Bxz = (tBxz - sum) - y;
	//
	sum = cube_split[i][j][k].Byx * Coeff[i][j][k].Byx1;
	y = ((ftypePML)cube[i][j][k].Ez - (ftypePML)cube[i - 1][j][k].Ez
		- ((ftypePML)compensator[i][j][k].Ez - (ftypePML)compensator[i - 1][j][k].Ez)) * Coeff[i][j][k].Byx2 * (_1dx)
				- compensatorPML[i][j][k].Byx * Coeff[i][j][k].Byx1;
	tByx = sum + y;
	compensatorPML2[i][j][k].Byx = (tByx - sum) - y;
	//
	sum = cube_split[i][j][k].Byz * Coeff[i][j][k].Byz1;
	y = - ((ftypePML)cube[i][j][k].Ex - (ftypePML)cube[i][j][k - 1].Ex
			- ((ftypePML)compensator[i][j][k].Ex - (ftypePML)compensator[i][j][k - 1].Ex)) * Coeff[i][j][k].Byz2 * (_1dz)
				- compensatorPML[i][j][k].Byz * Coeff[i][j][k].Byz1;
	tByz = sum + y;
	compensatorPML2[i][j][k].Byz = (tByz - sum) - y;
	//
	sum = cube_split[i][j][k].Bzx * Coeff[i][j][k].Bzx1;
	y = - ((ftypePML)cube[i][j][k].Ey - (ftypePML)cube[i - 1][j][k].Ey
			- ((ftypePML)compensator[i][j][k].Ey - (ftypePML)compensator[i - 1][j][k].Ey)) * Coeff[i][j][k].Bzx2 * (_1dx)
				- compensatorPML[i][j][k].Bzx * Coeff[i][j][k].Bzx1;
	tBzx = sum + y;
	compensatorPML2[i][j][k].Bzx = (tBzx - sum) - y;
	//
	sum = cube_split[i][j][k].Bzy * Coeff[i][j][k].Bzy1;
	y = ((ftypePML)cube[i][j][k].Ex - (ftypePML)cube[i][j - 1][k].Ex
		- ((ftypePML)compensator[i][j][k].Ex - (ftypePML)compensator[i][j - 1][k].Ex)) * Coeff[i][j][k].Bzy2 * (_1dy)
				- compensatorPML[i][j][k].Bzy * Coeff[i][j][k].Bzy1;
	tBzy = sum + y;
	compensatorPML2[i][j][k].Bzy = (tBzy - sum) - y;
	//
	cube_split[i][j][k].Bxy = tBxy;
	cube_split[i][j][k].Bxz = tBxz;
	cube_split[i][j][k].Byx = tByx;
	cube_split[i][j][k].Byz = tByz;
	cube_split[i][j][k].Bzx = tBzx;
	cube_split[i][j][k].Bzy = tBzy;

	cube[i][j][k].Bx = (ftype)(tBxy + tBxz);
	cube[i][j][k].By = (ftype)(tByx + tByz);
	cube[i][j][k].Bz = (ftype)(tBzx + tBzy);

	compensator2[i][j][k].Bx = (ftype)(((ftypePML)cube[i][j][k].Bx - tBxy) - tBxz);
	compensator2[i][j][k].By = (ftype)(((ftypePML)cube[i][j][k].By - tByx) - tByz);
	compensator2[i][j][k].Bz = (ftype)(((ftypePML)cube[i][j][k].Bz - tBzx) - tBzy);
}

template <class ftype>
void Update_electric_field_three_dimen_Kahan(vector<vector<vector<Component<ftype>>>>& cube,
	vector<vector<vector<Component<ftype>>>>& compensator, vector<vector<vector<Component<ftype>>>>& compensator2,
	ftype dt_x, ftype dt_y, ftype dt_z, int i, int j, int k)
{
	ftype tEx, tEy, tEz;
	ftype sum, y;
	//
	sum = cube[i][j][k].Ex;
	y = (cube[i][j + 1][k].Bz - cube[i][j][k].Bz - (compensator[i][j + 1][k].Bz - compensator[i][j][k].Bz)) * (dt_y)
		-(cube[i][j][k + 1].By - cube[i][j][k].By - (compensator[i][j][k + 1].By - compensator[i][j][k].By)) * (dt_z) - compensator[i][j][k].Ex;
	tEx = sum + y;
	compensator2[i][j][k].Ex = (tEx - sum) - y;
	//
	sum = cube[i][j][k].Ey;
	y = (cube[i][j][k + 1].Bx - cube[i][j][k].Bx - (compensator[i][j][k + 1].Bx - compensator[i][j][k].Bx)) * (dt_z)
		-(cube[i + 1][j][k].Bz - cube[i][j][k].Bz - (compensator[i + 1][j][k].Bz - compensator[i][j][k].Bz)) * (dt_x) - compensator[i][j][k].Ey;
	tEy = sum + y;
	compensator2[i][j][k].Ey = (tEy - sum) - y;
	//
	sum = cube[i][j][k].Ez;
	y = (cube[i + 1][j][k].By - cube[i][j][k].By - (compensator[i + 1][j][k].By - compensator[i][j][k].By)) * (dt_x)
		-(cube[i][j + 1][k].Bx - cube[i][j][k].Bx - (compensator[i][j + 1][k].Bx - compensator[i][j][k].Bx)) * (dt_y) - compensator[i][j][k].Ez;
	tEz = sum + y;
	compensator2[i][j][k].Ez = (tEz - sum) - y;

	cube[i][j][k].Ex = tEx;
	cube[i][j][k].Ey = tEy;
	cube[i][j][k].Ez = tEz;
}

template <class ftype>
void Update_magnetic_field_three_dimen_Kahan(vector<vector<vector<Component<ftype>>>>& cube, vector<vector<vector<Component<ftype>>>>& compensator, vector<vector<vector<Component<ftype>>>>& compensator2,
	ftype dt_x, ftype dt_y, ftype dt_z, int i, int j, int k)
{
	ftype tBx, tBy, tBz;
	ftype sum, y;

	sum = cube[i][j][k].Bx;
	y = -(cube[i][j][k].Ez - cube[i][j - 1][k].Ez - (compensator[i][j][k].Ez - compensator[i][j - 1][k].Ez)) * (dt_y)
		+(cube[i][j][k].Ey - cube[i][j][k - 1].Ey - (compensator[i][j][k].Ey - compensator[i][j][k - 1].Ey)) * (dt_z) - compensator[i][j][k].Bx;
	tBx = sum + y;
	compensator2[i][j][k].Bx = (tBx - sum) - y;
	//
	sum = cube[i][j][k].By;
	y = -(cube[i][j][k].Ex - cube[i][j][k - 1].Ex - (compensator[i][j][k].Ex - compensator[i][j][k - 1].Ex)) * (dt_z)
		+(cube[i][j][k].Ez - cube[i - 1][j][k].Ez - (compensator[i][j][k].Ez - compensator[i - 1][j][k].Ez)) * (dt_x) - compensator[i][j][k].By;
	tBy = sum + y;
	compensator2[i][j][k].By = (tBy - sum) - y;
	//
	sum = cube[i][j][k].Bz;
	y = -(cube[i][j][k].Ey - cube[i - 1][j][k].Ey - (compensator[i][j][k].Ey - compensator[i - 1][j][k].Ey)) * (dt_x)
		+(cube[i][j][k].Ex - cube[i][j - 1][k].Ex - (compensator[i][j][k].Ex - compensator[i][j - 1][k].Ex)) * (dt_y) - compensator[i][j][k].Bz;
	tBz = sum + y;
	compensator2[i][j][k].Bz = (tBz - sum) - y;

	cube[i][j][k].Bx = tBx;
	cube[i][j][k].By = tBy;
	cube[i][j][k].Bz = tBz;
}

template <class ftype>
void Add_Currents_3_dimen(vector< vector<vector<Component<ftype>>>>& cube,
	int Nx, int Ny, int Nz, int delta, ftype dx, ftype dy, ftype dz, const ftype dt, int it, int index_start, Direction direction,
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
double FDTD_three_dimen_with_PML(int Nx, int Ny, int Nz, ftype T, ftype dt, int delta, double sigma_x, double sigma_y,
	string type_sum, string file_energy, string file_coeff, string file_sigma)
{
	setlocale(LC_ALL, "Russian");

	int Nt = T / dt; int n = 4;
	// delta  =  width boundary layer

	vector<vector<vector<Component<ftype>>>> cube(Nx + 2 * delta + 2, vector<vector<Component<ftype>>>(Ny + 2 * delta + 2, vector<Component<ftype>>(Nz + 2 * delta + 2)));
	vector< vector<vector<ComponentSplit<ftypePML>>>> cube_split(Nx + 2 * delta + 2, vector< vector<ComponentSplit<ftypePML>>>(Ny + 2 * delta + 2, vector<ComponentSplit<ftypePML>>(Nz + 2 * delta + 2)));
	vector< vector<vector<SIGMA<double>>>> Sigma(Nx + 2 * delta + 2, vector< vector<SIGMA<double>>>(Ny + 2 * delta + 2, vector<SIGMA<double>>(Nz + 2 * delta + 2)));
	vector< vector<vector<COEFF<ftypePML>>>> Coeff(Nx + 2 * delta + 2, vector< vector<COEFF<ftypePML>>>(Ny + 2 * delta + 2, vector<COEFF<ftypePML>>(Nz + 2 * delta + 2)));

	vector<vector<vector<Component<ftype>>>> compensator(Nx + 2 + 2 * delta, vector<vector<Component<ftype>>>(Ny + 2 + 2 * delta, vector<Component<ftype>>(Nz + 2 + 2 * delta)));
	vector<vector<vector<Component<ftype>>>> compensator2(Nx + 2 + 2 * delta, vector<vector<Component<ftype>>>(Ny + 2 + 2 * delta, vector<Component<ftype>>(Nz + 2 + 2 * delta)));
	vector<vector<vector<ComponentSplit<ftypePML>>>> compensatorPML(Nx + 2 + 2 * delta, vector<vector<ComponentSplit<ftypePML>>>(Ny + 2 + 2 * delta, vector<ComponentSplit<ftypePML>>(Nz + 2 + 2 * delta)));
	vector<vector<vector<ComponentSplit<ftypePML>>>> compensatorPML2(Nx + 2 + 2 * delta, vector<vector<ComponentSplit<ftypePML>>>(Ny + 2 + 2 * delta, vector<ComponentSplit<ftypePML>>(Nz + 2 + 2 * delta)));

	pair<ftype, ftype> ab(0.0, 4.0 * M_PI), cd(0.0, 16.0 * M_PI), fg(0.0, 16.0 * M_PI);
	ftype dx = (ab.second - ab.first) / (ftype)Nx, dy = (cd.second - cd.first) / (ftype)Ny, dz = (fg.second - fg.first) / (ftype)Nz;

	ftype dt_x = dt / (dx), dt_y = dt / (dy), dt_z = dt / (dz);
	ftypePML _1dx = (ftypePML)1.0 / dx, _1dy = (ftypePML)1.0 / dy, _1dz = (ftypePML)1.0 / dz;
	vector<double> vec_energy(Nt + 1), vec_energy_acc(Nt + 1);

	Direction direction = x_comp_yz;


	// Initializing_FDTD_3_dimen_Gauss_PML<ftype>(cube, Nx, Ny, Nz, delta, dx, dy, dz, dt, delta + 1, direction, ab, cd, fg);
	Initializing_cube_split_3_dimen_PML<ftypePML>(cube_split, Nx, Ny, Nz, delta);
	Initializing_Sigma_3_dimen_PML<double>(Sigma, Nx, Ny, Nz, delta, n, sigma_x, sigma_y);
	Initializing_Coeff_3_dimen_PML_correct<ftypePML>(Coeff, Sigma, Nx, Ny, Nz, delta, dt);

	// Graph_for_Coeff_three_dimen<ftypePML>(Coeff, Nx, Ny, Nz, delta, dx, file_coeff);
    // Graph_for_Sigma_three_dimen<ftypePML>(Sigma, Nx, Ny, Nz, delta, dx, file_sigma);

	// Graph_Solution_in_two_planes_3dimen(cube, Nx, Ny, Nz, delta, dx, dy, dz, direction, "Bz PML T=8pi 0.csv");

	 double t1, t2;
	 t1 = omp_get_wtime();

	if (type_sum == "Kahan")
	{
		cout << type_sum << endl;
		for (int it = 0; it < Nt; it++)
		{
			if (it % 1000 == 0) cout << it << endl;

			Add_Currents_3_dimen<ftype>(cube, Nx, Ny, Nz, delta, dx, dy, dz, dt, it, delta + 1, direction, ab, cd, fg);
			vec_energy[0] = CalculateEnergyPML_3dimen(cube, Nx, Ny, Nz, delta);
			vec_energy_acc[0] = CalculateEnergyPML_3dimen_accurate(cube, compensator, Nx, Ny, Nz, delta);

#pragma omp parallel for collapse(3)
			for (int i = 1; i < Nx + 2 * delta + 1; i++)
				for (int j = 1; j < Ny + 2 * delta + 1; j++)
					for (int k = 1; k < Nz + 2 * delta + 1; k++)
					{
						if ((i >= delta + 1) && (i < Nx + delta + 1) && (j >= delta + 1) && (j < Ny + delta + 1) && (k >= delta + 1) && (k < Nz + delta + 1))
						{
							Update_electric_field_three_dimen_Kahan<ftype>(cube, compensator, compensator2, dt_x, dt_y, dt_z, i, j, k);
						}
						else {
							 // Update_electric_field_three_dimen_PML<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
							Update_electric_field_three_dimen_PML_Kahan<ftype, ftypePML>(cube, cube_split, Coeff, compensatorPML, compensatorPML2, compensator, compensator2, _1dx, _1dy, _1dz, i, j, k);

						}
					}
			for (int i = 1; i < Nx + 2 * delta + 1; i++) 
				for (int j = 1; j < Ny + 2 * delta + 1; j++)
					for (int k = 1; k < Nz + 2 * delta + 1; k++)
					{
						compensator[i][j][k].Ex = compensator2[i][j][k].Ex;
						compensator[i][j][k].Ey = compensator2[i][j][k].Ey;
						compensator[i][j][k].Ez = compensator2[i][j][k].Ez;

						compensatorPML[i][j][k].Exy = compensatorPML2[i][j][k].Exy;
						compensatorPML[i][j][k].Exz = compensatorPML2[i][j][k].Exz;
						compensatorPML[i][j][k].Eyx = compensatorPML2[i][j][k].Eyx;
						compensatorPML[i][j][k].Eyz = compensatorPML2[i][j][k].Eyz;
						compensatorPML[i][j][k].Ezx = compensatorPML2[i][j][k].Ezx;
						compensatorPML[i][j][k].Ezy = compensatorPML2[i][j][k].Ezy;
					}

#pragma omp parallel for collapse(3)
			for (int i = 1; i < Nx + 2 * delta + 1; i++)
				for (int j = 1; j < Ny + 2 * delta + 1; j++)
					for (int k = 1; k < Nz + 2 * delta + 1; k++)
					{
						if ((i >= delta + 1) && (i < Nx + delta + 1) && (j >= delta + 1) && (j < Ny + delta + 1) && (k >= delta + 1) && (k < Nz + delta + 1))
						{
							Update_magnetic_field_three_dimen_Kahan<ftype>(cube, compensator, compensator2, dt_x, dt_y, dt_z, i, j, k);
						}
						else {
							// Update_magnetic_field_three_dimen_PML<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
							Update_magnetic_field_three_dimen_PML_Kahan<ftype, ftypePML>(cube, cube_split, Coeff, compensatorPML, compensatorPML2, compensator, compensator2, _1dx, _1dy, _1dz, i, j, k);

						}
					}
			for (int i = 1; i < Nx + 2 * delta + 1; i++)
				for (int j = 1; j < Ny + 2 * delta + 1; j++)
					for (int k = 1; k < Nz + 2 * delta + 1; k++)
					{
						compensator[i][j][k].Bx = compensator2[i][j][k].Bx;
						compensator[i][j][k].By = compensator2[i][j][k].By;
						compensator[i][j][k].Bz = compensator2[i][j][k].Bz;

						compensatorPML[i][j][k].Bxy = compensatorPML2[i][j][k].Bxy;
						compensatorPML[i][j][k].Bxz = compensatorPML2[i][j][k].Bxz;
						compensatorPML[i][j][k].Byx = compensatorPML2[i][j][k].Byx;
						compensatorPML[i][j][k].Byz = compensatorPML2[i][j][k].Byz;
						compensatorPML[i][j][k].Bzx = compensatorPML2[i][j][k].Bzx;
						compensatorPML[i][j][k].Bzy = compensatorPML2[i][j][k].Bzy;
					}

			vec_energy[it + 1] = CalculateEnergyPML_3dimen(cube, Nx, Ny, Nz, delta);
			vec_energy_acc[it + 1] = CalculateEnergyPML_3dimen_accurate(cube, compensator, Nx, Ny, Nz, delta);
		}
	}
	else
	{
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

	t2 = omp_get_wtime();
	double result = vec_energy[Nt] / vec_energy[2512];

	cout << "Reflection coefficient = "<< endl << vec_energy[Nt] / vec_energy[2512] << endl << "Accurate reflection coefficient = " << endl << vec_energy_acc[Nt] / vec_energy_acc[2512] << endl << endl;
	cout << "time = " << t2 - t1 << endl << endl;

	// Graph_Solution_in_two_planes_3dimen(cube, Nx, Ny, Nz, delta, dx, dy, dz, direction, "Bz PML T=8pi 4.csv");
	// Graph_Solution_in_one_planes_3dimen(cube, Nx, Ny, Nz, delta, dx, dy, dz, direction);

	//����� ������� ������� � ����
	ofstream numb_energy(file_energy);
	numb_energy << ";" << "R" << ";" << ";" << ";" << "R_accurate"  << endl;

	for (int s = 0; s < Nt; s++)
	{
		numb_energy << dt * (double)(s + 1) << ";" << vec_energy[s] <<";"<< ";" << dt * (double)(s + 1) << ";" << vec_energy_acc[s] << endl;
	}
	numb_energy << endl << endl;
	numb_energy.close();

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
double CalculateEnergyPML_3dimen_accurate(vector<vector<vector<Component<ftype>>>>& cube, vector<vector<vector<Component<ftype>>>>& compensator, int Nx, int Ny, int Nz, int delta)
{
	double energy = 0.0;

	for (int i = 1; i < Nx + 2 * delta + 1; i++)
		for (int j = 1; j < Ny + 2 * delta + 1; j++)
			for (int k = 1; k < Nz + 2 * delta + 1; k++)
			{
				energy += ((double)cube[i][j][k].Ex - (double)compensator[i][j][k].Ex) * ((double)cube[i][j][k].Ex - (double)compensator[i][j][k].Ex)
					+ ((double)cube[i][j][k].Ey - (double)compensator[i][j][k].Ey) * ((double)cube[i][j][k].Ey - (double)compensator[i][j][k].Ey)
					+ ((double)cube[i][j][k].Ez - (double)compensator[i][j][k].Ez) * ((double)cube[i][j][k].Ez - (double)compensator[i][j][k].Ez);
				energy += ((double)cube[i][j][k].Bx - (double)compensator[i][j][k].Bx) * ((double)cube[i][j][k].Bx - (double)compensator[i][j][k].Bx)
					+ ((double)cube[i][j][k].By - (double)compensator[i][j][k].By) * ((double)cube[i][j][k].By - (double)compensator[i][j][k].By)
					+ ((double)cube[i][j][k].Bz - (double)compensator[i][j][k].Bz) * ((double)cube[i][j][k].Bz - (double)compensator[i][j][k].Bz);
			}
	return energy;
}

template <class ftype>
void Graph_Solution_in_two_planes_3dimen(vector< vector<vector<Component<ftype>>>>& cube, int Nx, int Ny, int Nz, int delta, ftype dx, ftype dy, ftype dz, Direction direction, string file_name)
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

template <class ftype>
void Graph_Solution_in_one_planes_3dimen(vector< vector<vector<Component<ftype>>>>& cube, int Nx, int Ny, int Nz, int delta, ftype dx, ftype dy, ftype dz, Direction direction)
{
	ofstream numbEx("graph_solution_E_x_one_dimen_float-float-sigma 49.6.csv"), numbBx("graph_solution_B_x_one_dimen_float-float-sigma 49.6.csv");
	ofstream numbEy("graph_solution_E_y_one_dimen_float-float-sigma 49.6.csv"), numbBy("graph_solution_B_y_one_dimen_float-float-sigma 49.6.csv");
	ofstream numbEz("graph_solution_E_z_one_dimen_float-float-sigma 49.6.csv"), numbBz("graph_solution_B_z_one_dimen_float-float-sigma 49.6.csv");

	for (int i = 1; i < Nx + 2 * delta + 1; i++)
	{
		numbEx << dx * (double)(i) << ";" << cube[i][Ny / 2 + delta][Nz / 2 + delta].Ex << endl;
		numbEy << dx * (double)(i) << ";" << cube[i][Ny / 2 + delta][Nz / 2 + delta].Ey << endl;
		numbEz << dx * (double)(i) << ";" << cube[i][Ny / 2 + delta][Nz / 2 + delta].Ez << endl;

		numbBx << dx * (double)(i) << ";" << cube[i][Ny / 2 + delta][Nz / 2 + delta].Bx << endl;
		numbBy << dx * (double)(i) << ";" << cube[i][Ny / 2 + delta][Nz / 2 + delta].By << endl;
		numbBz << dx * (double)(i) << ";" << cube[i][Ny / 2 + delta][Nz / 2 + delta].Bz << endl;

	}

	numbEx.close();
	numbBx.close();
	numbEy.close();
	numbBy.close();
	numbEz.close();
	numbBz.close();
}

