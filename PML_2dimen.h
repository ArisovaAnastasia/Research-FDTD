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
enum Direction {
	x_comp_yz, x_comp_zy, y_comp_xz, y_comp_zx, z_comp_xy, z_comp_yx,
};
template <class ftype>
class SIGMA
{
public:
	ftype sigmaE_x, sigmaE_y, sigmaE_z, sigmaH_x, sigmaH_y, sigmaH_z;
};

template <class ftype>
class COEFF
{
public:
	ftype Exy1, Exz1, Eyx1, Eyz1, Ezx1, Ezy1, Bxy1, Bxz1, Byx1, Byz1, Bzx1, Bzy1,
		Exy2, Exz2, Eyx2, Eyz2, Ezx2, Ezy2, Bxy2, Bxz2, Byx2, Byz2, Bzx2, Bzy2;
};

template <class ftype>
class ComponentSplit
{
public:
	ftype Exy, Exz, Eyx, Eyz, Ezx, Ezy,
		Bxy, Bxz, Byx, Byz, Bzx, Bzy;
};

// #define R 1.0e-8

template <class ftype>
ftype Sigma_max(int n, double R, int delta, ftype step)
{
	return -(ftype)(n + 1)*(ftype)log ((ftype)R)/(2.0 * (ftype)delta * step);
}

template <class ftype>
ftype distanceH(int N, int delta, int i)
{
	ftype dist = 0.0;
	if (i < delta + 1)
		dist = (ftype)(delta + 1 - i);

	if(i > N +  delta)
		dist = (ftype)(i - delta - N) - (ftype)0.5;

	return dist / (ftype)delta;
}

template <class ftype>
ftype distanceE(int N, int delta, int i)
{ //���� ���������� ������������ PML-����
	ftype dist = 0.0;
	if (i < delta + 1)
		dist = (ftype)(delta + 1 - i) - (ftype)0.5;

	if (i > N + delta)
		dist = (ftype)(i - delta - N);

	return dist / (ftype)delta;
}

template <class ftype>
void Initializing_FDTD_two_dimen_Gauss_PML(vector<vector<Component<ftype>>>& cube, vector<vector<ComponentSplit<ftype>>>& cube_split, int Nx, int Ny, int delta, ftype dx, ftype dy, const ftype dt,
	int index_start, Direction direction,	pair <ftype, ftype> ab, pair <ftype, ftype> cd, double& abs_err_t, double& rel_err_t, double& smape_err_t)
{
	ftype x1, x2, y1, y2, t;
	//ftype ay = 8. * M_PI / 3.0; // size of beam in y direction
	//ftype ax = 3. * M_PI / 3.0; // size of beam in x direction

	ftype ax_longitudinal = ab.second / 6.; // продольное
	ftype ay_longitudinal = cd.second / 6.;
	ftype ax_transverse = ab.second * 3. / 4.; // поперечное
	ftype ay_transverse = cd.second * 3. / 4.;
	double lambda_ = 2. * M_PI / 4.; // wvelength
	double k_ = 2. * M_PI / lambda_; // wavenumber = 2 pi/ wavelength
	double omega_ = 2. * M_PI / lambda_; // 2 pi c /wavelength 
	double t0_ = 0; // initial moment of time
	t = (ftype)0.5 * dt; // timeshift between components

	float x0 = (ab.second - ab.first) / 2.;
	float y0 = (cd.second - cd.first) / 2.;

	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
		{
			x1 = (ftype)(ab.first + dx * (ftype)i + (ftype)0.5 * dx);
			x2 = (ftype)(ab.first + dx * (ftype)i);

			y1 = (ftype)(cd.first + dy * (ftype)j + (ftype)0.5 * dy);
			y2 = (ftype)(cd.first + dy * (ftype)j);

			switch (direction)
			{
			case x_comp_zy: {
				// Ez and By has the same y coordinate
				cube[i + index_start][j + index_start].Ez = exp(-(y2 - y0) * (y2 - y0) / (ax_transverse * ax_transverse)) * exp(-(x1 - x0) * (x1 - x0) / (ax_longitudinal * ax_longitudinal)) *
					sin(omega_ * t0_ + k_ * x1); // sin(phase), phase = omega * t - k * x
				cube[i + index_start][j + index_start].By = exp(-(y2 - y0) * (y2 - y0) / (ax_transverse * ax_transverse)) * exp(-(x2 - x0) * (x2 - x0) / (ax_longitudinal * ax_longitudinal)) *
					sin(omega_ * (t0_ + t) + k_ * x2); // sin(phase), phase = omega * t - k * x
				break;
			}
			case x_comp_yz: {
				// Ez and By has the same y coordinate
				cube[i + index_start][j + index_start].Ey = exp(-(y2 - y0) * (y2 - y0) / (ax_transverse * ax_transverse)) * exp(-(x1 - x0) * (x1 - x0) / (ax_longitudinal * ax_longitudinal)) *
					sin(omega_ * t0_ + k_ * x1); // sin(phase), phase = omega * t - k * x
				cube[i + index_start][j + index_start].Bz = exp(-(y2 - y0) * (y2 - y0) / (ax_transverse * ax_transverse)) * exp(-(x2 - x0) * (x2 - x0) / (ax_longitudinal * ax_longitudinal)) *
					sin(omega_ * (t0_ + t) + k_ * x2); // sin(phase), phase = omega * t - k * x
				break;
			}
			case y_comp_zx: {
				// Ex and Bz has the same x coordinate
				cube[i + index_start][j + index_start].Ez = exp(-(y1 - y0) * (y1 - y0) / (ay_longitudinal * ay_longitudinal)) * exp(-(x2 - x0) * (x2 - x0) / (ay_transverse * ay_transverse)) *
					sin(omega_ * t0_ - k_ * y1); // sin(phase), phase = omega * t - k * x
				cube[i + index_start][j + index_start].Bx = exp(-(y2 - y0) * (y2 - y0) / (ay_longitudinal * ay_longitudinal)) * exp(-(x2 - x0) * (x2 - x0) / (ay_transverse * ay_transverse)) *
					sin(omega_ * (t0_ + t) - k_ * y2); // sin(phase), phase = omega * t - k * x
				break;
			}
			case y_comp_xz: {
				// Ex and Bz has the same x coordinate
				cube[i + index_start][j + index_start].Ex = exp(-(y1 - y0) * (y1 - y0) / (ay_longitudinal * ay_longitudinal)) * exp(-(x2 - x0) * (x2 - x0) / (ay_transverse * ay_transverse)) *
					sin(omega_ * t0_ - k_ * y1); // sin(phase), phase = omega * t - k * x
				cube[i + index_start][j + index_start].Bz = exp(-(y2 - y0) * (y2 - y0) / (ay_longitudinal * ay_longitudinal)) * exp(-(x2 - x0) * (x2 - x0) / (ay_transverse * ay_transverse)) *
					sin(omega_ * (t0_ + t) - k_ * y2); // sin(phase), phase = omega * t - k * x
				break;
			}
			default:
				break;
			}
		}
	// Border conditions
	for (int i = 0; i < Ny + 2 * delta + 2; i++) {
		cube[0][i].Ex = cube[Nx + 2 * delta][i].Ex;
		cube[0][i].Ey = cube[Nx + 2 * delta][i].Ey;
		cube[0][i].Ez = cube[Nx + 2 * delta][i].Ez;
		cube[0][i].Bx = cube[Nx + 2 * delta][i].Bx;
		cube[0][i].By = cube[Nx + 2 * delta][i].By;
		cube[0][i].Bz = cube[Nx + 2 * delta][i].Bz;
		cube[Nx + 2 * delta + 1][i].Ex = cube[1][i].Ex;
		cube[Nx + 2 * delta + 1][i].Ey = cube[1][i].Ey;
		cube[Nx + 2 * delta + 1][i].Ez = cube[1][i].Ez;
		cube[Nx + 2 * delta + 1][i].Bx = cube[1][i].Bx;
		cube[Nx + 2 * delta + 1][i].By = cube[1][i].By;
		cube[Nx + 2 * delta + 1][i].Bz = cube[1][i].Bz;
	}
	for (int i = 0; i < Nx + 2 * delta + 2; i++) {
		cube[i][0].Ex = cube[i][Ny + 2 * delta].Ex;
		cube[i][0].Ey = cube[i][Ny + 2 * delta].Ey;
		cube[i][0].Ez = cube[i][Ny + 2 * delta].Ez;
		cube[i][0].Bx = cube[i][Ny + 2 * delta].Bx;
		cube[i][0].By = cube[i][Ny + 2 * delta].By;
		cube[i][0].Bz = cube[i][Ny + 2 * delta].Bz;
		cube[i][Ny + 2 * delta + 1].Ex = cube[i][1].Ex;
		cube[i][Ny + 2 * delta + 1].Ey = cube[i][1].Ey;
		cube[i][Ny + 2 * delta + 1].Ez = cube[i][1].Ez;
		cube[i][Ny + 2 * delta + 1].Bx = cube[i][1].Bx;
		cube[i][Ny + 2 * delta + 1].By = cube[i][1].By;
		cube[i][Ny + 2 * delta + 1].Bz = cube[i][1].Bz;
	}
}

template <class ftype>
void Initializing_cube_split_two_dimen_PML(vector<vector<ComponentSplit<ftype>>>& cube_split, int Nx, int Ny, int delta)
{
	for (int i = 0; i < Nx + 2 * delta +2; i++)
		for (int j = 0; j < Ny + 2 * delta +2; j++)
		{
			cube_split[i][j].Ezx = (ftype)0.0;
			cube_split[i][j].Ezy = (ftype)0.0;
			cube_split[i][j].Bzx = (ftype)0.0;
			cube_split[i][j].Bzy = (ftype)0.0;
		}
}

template <class ftype>
void Initializing_Sigma_two_dimen_PML(vector<vector<SIGMA<ftype>>>& Sigma, int Nx, int Ny, ftype dx, ftype dy,
	int delta, int n, double R)
{
	ftype dz = 0.0;
	ftype var_Sigma_max_x = Sigma_max<ftype>(n, R, delta, dx);
	ftype var_Sigma_max_y = Sigma_max<ftype>(n, R, delta, dy);
	ftype var_Sigma_max_z = (ftype)0.0;

	cout << var_Sigma_max_x << "  " << var_Sigma_max_y << "   " << var_Sigma_max_z << "    " << endl;

	for(int i=1; i< Nx + 2 * delta + 1; i++)
		for (int j = 1; j < Ny + 2 * delta + 1; j++)
		{
			Sigma[i][j].sigmaE_x = var_Sigma_max_x * pow(distanceE<ftype>(Nx, delta, i), n);
			Sigma[i][j].sigmaH_x = var_Sigma_max_x * pow(distanceH<ftype>(Nx, delta, i), n);

			Sigma[i][j].sigmaE_y = var_Sigma_max_y * pow(distanceE<ftype>(Ny, delta, j), n);
			Sigma[i][j].sigmaH_y = var_Sigma_max_y * pow(distanceH<ftype>(Ny, delta, j), n);

			Sigma[i][j].sigmaE_z = var_Sigma_max_z;
			Sigma[i][j].sigmaH_z = var_Sigma_max_z;
		}

	//for (int i = 1; i < delta + 1; i++)//����� ������������� ������� PML
	//{
	//	for (int j = 1; j < delta + 1; j++)
	//	{
	//		Sigma[i][j].sigmaE_x = var_Sigma_max_x * pow(distanceE_x<ftype>(Nx, delta, i), n);
	//		Sigma[i][j].sigmaH_x = var_Sigma_max_x * pow(distanceH_x<ftype>(Nx, delta, i), n);

	//		Sigma[i][j].sigmaE_y = var_Sigma_max_y * pow(distanceE_y<ftype>(Ny, delta, j), n);
	//		Sigma[i][j].sigmaH_y = var_Sigma_max_y * pow(distanceH_y<ftype>(Ny, delta, j), n);
	//	}
	//	for (int j = delta + 1; j < Ny + delta + 1; j++)
	//	{
	//		Sigma[i][j].sigmaE_x = var_Sigma_max_x * pow(distanceE_x<ftype>(Nx, delta, i), n);
	//		Sigma[i][j].sigmaH_x = var_Sigma_max_x * pow(distanceH_x<ftype>(Nx, delta, i), n);
	//	}
	//	for (int j = Ny + delta + 1; j < Ny + 2 * delta + 1; j++)
	//	{
	//		Sigma[i][j].sigmaE_x = var_Sigma_max_x * pow(distanceE_x<ftype>(Nx, delta, i), n);
	//		Sigma[i][j].sigmaH_x = var_Sigma_max_x * pow(distanceH_x<ftype>(Nx, delta, i), n);

	//		Sigma[i][j].sigmaE_y = var_Sigma_max_y * pow(distanceE_y<ftype>(Ny, delta, j), n);
	//		Sigma[i][j].sigmaH_y = var_Sigma_max_y * pow(distanceH_y<ftype>(Ny, delta, j), n);
	//	}
	//}
	//for (int i = delta + 1; i < Nx + delta + 1; i++)//������� � ������ ������������� ������� PML ��� �����
	//{
	//	for (int j = 1; j < delta + 1; j++)
	//	{
	//		Sigma[i][j].sigmaE_y = var_Sigma_max_y * pow(distanceE_y<ftype>(Ny, delta, j), n);
	//		Sigma[i][j].sigmaH_y = var_Sigma_max_y * pow(distanceH_y<ftype>(Ny, delta, j), n);
	//	}
	//	for (int j = Ny + delta + 1; j < Ny + 2 * delta + 1; j++)
	//	{
	//		Sigma[i][j].sigmaE_y = var_Sigma_max_y * pow(distanceE_y<ftype>(Ny, delta, j), n);
	//		Sigma[i][j].sigmaH_y = var_Sigma_max_y * pow(distanceH_y<ftype>(Ny, delta, j), n);
	//	}
	//}
	//for (int i = Nx + delta + 1; i < Nx + 2 * delta + 1; i++)//������ ������������� ������� PML
	//{
	//	for (int j = 1; j < delta + 1; j++)
	//	{
	//		Sigma[i][j].sigmaE_x = var_Sigma_max_x * pow(distanceE_x<ftype>(Nx, delta, i), n);
	//		Sigma[i][j].sigmaH_x = var_Sigma_max_x * pow(distanceH_x<ftype>(Nx, delta, i), n);

	//		Sigma[i][j].sigmaE_y = var_Sigma_max_y * pow(distanceE_y<ftype>(Ny, delta, j), n);
	//		Sigma[i][j].sigmaH_y = var_Sigma_max_y * pow(distanceH_y<ftype>(Ny, delta, j), n);
	//	}
	//	for (int j = delta + 1; j < Ny + delta + 1; j++)
	//	{
	//		Sigma[i][j].sigmaE_x = var_Sigma_max_x * pow(distanceE_x<ftype>(Nx, delta, i), n);
	//		Sigma[i][j].sigmaH_x = var_Sigma_max_x * pow(distanceH_x<ftype>(Nx, delta, i), n);
	//	}
	//	for (int j = Ny + delta + 1; j < Ny + 2 * delta + 1; j++)
	//	{
	//		Sigma[i][j].sigmaE_x = var_Sigma_max_x * pow(distanceE_x<ftype>(Nx, delta, i), n);
	//		Sigma[i][j].sigmaH_x = var_Sigma_max_x * pow(distanceH_x<ftype>(Nx, delta, i), n);

	//		Sigma[i][j].sigmaE_y = var_Sigma_max_y * pow(distanceE_y<ftype>(Ny, delta, j), n);
	//		Sigma[i][j].sigmaH_y = var_Sigma_max_y * pow(distanceH_y<ftype>(Ny, delta, j), n);
	//	}
	//}
}

template <class ftype>
void Initializing_Coeff_two_dimen_PML(vector<vector<COEFF<ftype>>>& Coeff, vector<vector<SIGMA<ftype>>>& Sigma, int Nx, int Ny, int delta, ftype dt)
{
	for (int i = 1; i < Nx + 2 * delta + 1; i++)
		for (int j = 1; j < Ny + 2 * delta + 1; j++)
		{
			if ((i >= delta + 1) && (i < Nx + delta + 1) && (j >= delta + 1) && (j < Ny + delta + 1))
			{
			}
			else {
				Coeff[i][j].Exy1 = (ftype)exp(-dt * Sigma[i][j].sigmaE_y);
				Coeff[i][j].Exz1 = (ftype)exp(-dt * Sigma[i][j].sigmaE_z);
				Coeff[i][j].Eyx1 = (ftype)exp(-dt * Sigma[i][j].sigmaE_x);
				Coeff[i][j].Eyz1 = (ftype)exp(-dt * Sigma[i][j].sigmaE_z);
				Coeff[i][j].Ezx1 = (ftype)exp(-dt * Sigma[i][j].sigmaE_x);
				Coeff[i][j].Ezy1 = (ftype)exp(-dt * Sigma[i][j].sigmaE_y);

				Coeff[i][j].Bxy1 = (ftype)exp(-dt * Sigma[i][j].sigmaH_y);
				Coeff[i][j].Bxz1 = (ftype)exp(-dt * Sigma[i][j].sigmaH_z);
				Coeff[i][j].Byx1 = (ftype)exp(-dt * Sigma[i][j].sigmaH_x);
				Coeff[i][j].Byz1 = (ftype)exp(-dt * Sigma[i][j].sigmaH_z);
				Coeff[i][j].Bzx1 = (ftype)exp(-dt * Sigma[i][j].sigmaH_x);
				Coeff[i][j].Bzy1 = (ftype)exp(-dt * Sigma[i][j].sigmaH_y);

				if (Sigma[i][j].sigmaE_x != (ftype)0.0)
				{
					Coeff[i][j].Eyx2 = (ftype)1.0 / Sigma[i][j].sigmaE_x * ((ftype)1.0 - Coeff[i][j].Eyx1);
					Coeff[i][j].Ezx2 = (ftype)1.0 / Sigma[i][j].sigmaE_x * ((ftype)1.0 - Coeff[i][j].Ezx1);
				}
				else {
					Coeff[i][j].Eyx2 = dt;
					Coeff[i][j].Ezx2 = dt;
				}
				if (Sigma[i][j].sigmaE_y != (ftype)0.0)
				{
					Coeff[i][j].Exy2 = (ftype)1.0 / Sigma[i][j].sigmaE_y * ((ftype)1.0 - Coeff[i][j].Exy1);
					Coeff[i][j].Ezy2 = (ftype)1.0 / Sigma[i][j].sigmaE_y * ((ftype)1.0 - Coeff[i][j].Ezy1);
				}
				else {
					Coeff[i][j].Exy2 = dt;
					Coeff[i][j].Ezy2 = dt;
				}
				if (Sigma[i][j].sigmaE_z != (ftype)0.0)
				{
					Coeff[i][j].Exz2 = (ftype)1.0 / Sigma[i][j].sigmaE_z * ((ftype)1.0 - Coeff[i][j].Exz1);
					Coeff[i][j].Eyz2 = (ftype)1.0 / Sigma[i][j].sigmaE_z * ((ftype)1.0 - Coeff[i][j].Eyz1);
				}
				else {
					Coeff[i][j].Exz2 = dt;
					Coeff[i][j].Eyz2 = dt;
				}
				if (Sigma[i][j].sigmaH_x != (ftype)0.0)
				{
					Coeff[i][j].Byx2 = (ftype)1.0 / Sigma[i][j].sigmaH_x * ((ftype)1.0 - Coeff[i][j].Byx1);
					Coeff[i][j].Bzx2 = (ftype)1.0 / Sigma[i][j].sigmaH_x * ((ftype)1.0 - Coeff[i][j].Bzx1);
				}
				else {
					Coeff[i][j].Byx2 = dt;
					Coeff[i][j].Bzx2 = dt;
				}
				if (Sigma[i][j].sigmaH_y != (ftype)0.0)
				{
					Coeff[i][j].Bxy2 = (ftype)1.0 / Sigma[i][j].sigmaH_y * ((ftype)1.0 - Coeff[i][j].Bxy1);
					Coeff[i][j].Bzy2 = (ftype)1.0 / Sigma[i][j].sigmaH_y * ((ftype)1.0 - Coeff[i][j].Bzy1);
				}
				else {
					Coeff[i][j].Bxy2 = dt;
					Coeff[i][j].Bzy2 = dt;
				}
				if (Sigma[i][j].sigmaH_z != (ftype)0.0)
				{
					Coeff[i][j].Bxz2 = (ftype)1.0 / Sigma[i][j].sigmaH_z * ((ftype)1.0 - Coeff[i][j].Bxz1);
					Coeff[i][j].Byz2 = (ftype)1.0 / Sigma[i][j].sigmaH_z * ((ftype)1.0 - Coeff[i][j].Byz1);
				}
				else {
					Coeff[i][j].Bxz2 = dt;
					Coeff[i][j].Byz2 = dt;
				}
			}
		}
}

template <class ftype>
void Update_electric_field_two_dimen_PML(vector<vector<Component<ftype>>>& cube,
	vector<vector<ComponentSplit<ftype>>>& cube_split, vector<vector<COEFF<ftype>>>& Coeff,
	ftype _1dx, ftype _1dy, int i, int j)
{
	ftype tExy, tExz, tEyx, tEyz, tEzx, tEzy;
	// ����� ���������� �� �������� ����, �� ������������� ���� ������ ������������
	// ftype -> ftypePML  (static cast)

	tExy = cube_split[i][j].Exy * Coeff[i][j].Exy1 + (cube[i][j + 1].Bz - cube[i][j].Bz) * Coeff[i][j].Exy2 * (_1dy);

	tExz = cube_split[i][j].Exz * Coeff[i][j].Exz1;

	tEyx = cube_split[i][j].Eyx * Coeff[i][j].Eyx1 - (cube[i + 1][j].Bz - cube[i][j].Bz) * Coeff[i][j].Eyx2 * (_1dx);

	tEyz = cube_split[i][j].Eyz * Coeff[i][j].Eyz1;

	tEzx = cube_split[i][j].Ezx * Coeff[i][j].Ezx1 + (cube[i + 1][j].By - cube[i][j].By) * Coeff[i][j].Ezx2 * (_1dx);

	tEzy = cube_split[i][j].Ezy * Coeff[i][j].Ezy1 - (cube[i][j + 1].Bx - cube[i][j].Bx) * Coeff[i][j].Ezy2 * (_1dy);

	cube_split[i][j].Exy = tExy;
	cube_split[i][j].Exz = tExz;
	cube_split[i][j].Eyx = tEyx;
	cube_split[i][j].Eyz = tEyz;
	cube_split[i][j].Ezx = tEzx;
	cube_split[i][j].Ezy = tEzy;
	cube[i][j].Ex = tExy + tExz;
	cube[i][j].Ey = tEyx + tEyz;
	cube[i][j].Ez = tEzx + tEzy;
}

template <class ftype>
void Update_magnetic_field_two_dimen_PML(vector<vector<Component<ftype>>>& cube,
	vector<vector<ComponentSplit<ftype>>>& cube_split, vector<vector<COEFF<ftype>>>& Coeff,
	ftype _1dx, ftype _1dy, int i, int j)
{
	ftype tBxy, tBxz, tByx, tByz, tBzx, tBzy;

	tBxy = cube_split[i][j].Bxy * Coeff[i][j].Bxy1 - (cube[i][j].Ez - cube[i][j - 1].Ez) * Coeff[i][j].Bxy2 * (_1dy);

	tBxz = cube_split[i][j].Bxz * Coeff[i][j].Bxz1;

	tByx = cube_split[i][j].Byx * Coeff[i][j].Byx1 + (cube[i][j].Ez - cube[i - 1][j].Ez) * Coeff[i][j].Byx2 * (_1dx);

	tByz = cube_split[i][j].Byz * Coeff[i][j].Byz1;

	tBzx = cube_split[i][j].Bzx * Coeff[i][j].Bzx1 - (cube[i][j].Ey - cube[i - 1][j].Ey) * Coeff[i][j].Bzx2 * (_1dx);

	tBzy = cube_split[i][j].Bzy * Coeff[i][j].Bzy1 + (cube[i][j].Ex - cube[i][j - 1].Ex) * Coeff[i][j].Bzy2 * (_1dy);
	
	cube_split[i][j].Bxy = tBxy;
	cube_split[i][j].Bxz = tBxz;
	cube_split[i][j].Byx = tByx;
	cube_split[i][j].Byz = tByz;
	cube_split[i][j].Bzx = tBzx;
	cube_split[i][j].Bzy = tBzy;
	cube[i][j].Bx = tBxy + tBxz;
	cube[i][j].By = tByx + tByz;
	cube[i][j].Bz = tBzx + tBzy;
}

template <class ftype>
void Add_Currents_2_dimen(vector<vector<Component<ftype>>>& cube,
	int Nx, int Ny, int delta, ftype dx, ftype dy, const ftype dt, int it, int index_start, Direction direction,
	pair <ftype, ftype> ab, pair <ftype, ftype> cd)
{
	ftype x1, x2, y1, y2, t1, t2;

	ftype ax_transverse = ab.second * 3. / 4.; // поперечное
	ftype ay_transverse = cd.second * 3. / 4.;

	ftype tp_x = ab.second / 6.;
	ftype tp_y = cd.second / 6.;

	double lambda_ = 2. * M_PI / 4.; // wvelength               //?
	double k_ = 2. * M_PI / lambda_; // wavenumber = 2 pi/ wavelength
	double omega_ = 2. * M_PI / lambda_; // 2 pi c /wavelength 
	float A = 1.; //amplitude

	t1 = dt * (ftype)it;
	t2 = t1 + (ftype)0.5 * dt;

	float x0 = (ab.second - ab.first) / 2.;
	float y0 = (cd.second - cd.first) / 2.;

	float t0_x = 3. * tp_x;
	float t0_y = 3. * tp_y;

	int offset = 5;
	x1 = (ftype)(ab.first + dx * (ftype)offset + (ftype)0.5 * dx);
	x2 = (ftype)(ab.first + dx * (ftype)offset);

	for (int j = 0; j < Ny; j++)
	{
		y1 = (ftype)(cd.first + dy * (ftype)j + (ftype)0.5 * dy);
		y2 = (ftype)(cd.first + dy * (ftype)j);

		// Ez and By has the same y coordinate
		cube[offset + index_start][j + index_start].Ey += A * exp(-(t1 - t0_x) * (t1 - t0_x) / (tp_x * tp_x))
			* exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse)) *
			sin(omega_ * t1 - k_ * x1); // sin(phase), phase = omega * t - k * x
		cube[offset + index_start][j + index_start].Bz += A * exp(-(t2 - t0_x) * (t2 - t0_x) / (tp_x * tp_x))
			* exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse)) *
			sin(omega_ * t2 - k_ * x2); // sin(phase), phase = omega * t - k * x
		}
}

template <class ftype>
double FDTD_two_dimen_with_PML(int Nx, int Ny, ftype T, ftype dt, int delta, double R, string type_sum, string name_file)
{
	setlocale(LC_ALL, "Russian");

	int Nt = T / dt; int n = 4;
	// delta  =  width boundary layer

	vector<vector<Component<ftype>>> cube(Nx + 2 * delta + 2, vector<Component<ftype>>(Ny + 2 * delta + 2));
	vector<vector<ComponentSplit<ftype>>> cube_split(Nx + 2 * delta + 2, vector<ComponentSplit<ftype>>(Ny + 2 * delta + 2));
	vector<vector<SIGMA<ftype>>> Sigma(Nx + 2 * delta + 2, vector<SIGMA<ftype>>(Ny + 2 * delta + 2));
	vector<vector<COEFF<ftype>>> Coeff(Nx + 2 * delta + 2, vector<COEFF<ftype>>(Ny + 2 * delta + 2));

	pair<ftype, ftype> ab(0.0, 4.0 * M_PI), cd(0.0, 16.0 * M_PI);
	ftype dx = (ab.second - ab.first) / (ftype)Nx, dy = (cd.second - cd.first) / (ftype)Ny;

	Direction direction = x_comp_yz;

	ftype dt_x = dt / (dx), dt_y = dt / (dy);
	ftype _1dx = (ftype)1.0 / dx, _1dy = (ftype)1.0 / dy;
	double abs_err_t = 0.0, rel_err_t = 0.0, smape_err_t = 0.0;
	vector<double> vec_energy(Nt + 1);

	// Initializing_FDTD_two_dimen_Gauss_PML<ftype>(cube, cube_split, Nx, Ny, delta, dx, dy, dt, delta + 1, direction, ab, cd, abs_err_t, rel_err_t, smape_err_t);
	Initializing_cube_split_two_dimen_PML<ftype>(cube_split, Nx, Ny, delta);
	Initializing_Sigma_two_dimen_PML<ftype>(Sigma, Nx, Ny, dx, dy, delta, n, R);
	Initializing_Coeff_two_dimen_PML<ftype>(Coeff, Sigma, Nx, Ny, delta, dt);


	if (type_sum == "Kahan")
	{
	}
	else
	{
		vec_energy[0] = CalculateEnergyPML_2dimen(cube, Nx, Ny, delta);
		for (int it = 0; it < Nt; it++)
		{
			if (it % 1000 == 0)
				cout << it << endl;
			
			Add_Currents_2_dimen(cube, Nx, Ny, delta, dx, dy, dt, it, delta + 1, direction, ab, cd);

#pragma omp parallel for collapse(2)
			for (int i = 1; i < Nx + 2 * delta + 1; i++)
				for (int j = 1; j < Ny + 2 * delta + 1; j++)
					{
						if ((i >= delta + 1) && (i < Nx + delta + 1)&& (j >= delta + 1) && (j < Ny + delta + 1))
						{
							Update_electric_field_two_dimen<ftype>(cube, dt_x, dt_y, i, j);
							//cout << "Update_electric_field_two_dimen  " << i << " " << j << endl;
						}
						else {
							Update_electric_field_two_dimen_PML<ftype>(cube, cube_split, Coeff, _1dx, _1dy, i, j);
							// cout << "Update_electric_field_two_dimen PML  " << i << " " << j << endl;
						}
					}

#pragma omp parallel for collapse(2)
			for (int i = 1; i < Nx + 2 * delta + 1; i++)
				for (int j = 1; j < Ny + 2 * delta + 1; j++)
					{
						if ((i >= delta + 1) && (i < Nx + delta + 1) && (j >= delta + 1) && (j < Ny + delta + 1))
						{
							Update_magnetic_field_two_dimen<ftype>(cube, dt_x, dt_y, i, j);
							//cout << "Update_magnetic_field_two_dimen  " << i << " " << j << endl;
						}
						else {
							Update_magnetic_field_two_dimen_PML<ftype>(cube, cube_split, Coeff, _1dx, _1dy, i, j);
							//cout << "Update_magnetic_field_two_dimen PML " << i << " " << j << endl;
						}
					}
			// Border conditions
			//for (int i = 0; i < Ny + 2 * delta + 2; i++) {
			//	cube[0][i].Ex = cube[Nx + 2 * delta][i].Ex;
			//	cube[0][i].Ey = cube[Nx + 2 * delta][i].Ey;
			//	cube[0][i].Ez = cube[Nx + 2 * delta][i].Ez;
			//	cube[0][i].Bx = cube[Nx + 2 * delta][i].Bx;
			//	cube[0][i].By = cube[Nx + 2 * delta][i].By;
			//	cube[0][i].Bz = cube[Nx + 2 * delta][i].Bz;
			//	cube[Nx + 2 * delta + 1][i].Ex = cube[1][i].Ex;
			//	cube[Nx + 2 * delta + 1][i].Ey = cube[1][i].Ey;
			//	cube[Nx + 2 * delta + 1][i].Ez = cube[1][i].Ez;
			//	cube[Nx + 2 * delta + 1][i].Bx = cube[1][i].Bx;
			//	cube[Nx + 2 * delta + 1][i].By = cube[1][i].By;
			//	cube[Nx + 2 * delta + 1][i].Bz = cube[1][i].Bz;
			//}
			//for (int i = 0; i < Nx + 2 * delta + 2; i++) {
			//	cube[i][0].Ex = cube[i][Ny + 2 * delta].Ex;
			//	cube[i][0].Ey = cube[i][Ny + 2 * delta].Ey;
			//	cube[i][0].Ez = cube[i][Ny + 2 * delta].Ez;
			//	cube[i][0].Bx = cube[i][Ny + 2 * delta].Bx;
			//	cube[i][0].By = cube[i][Ny + 2 * delta].By;
			//	cube[i][0].Bz = cube[i][Ny + 2 * delta].Bz;
			//	cube[i][Ny + 2 * delta + 1].Ex = cube[i][1].Ex;
			//	cube[i][Ny + 2 * delta + 1].Ey = cube[i][1].Ey;
			//	cube[i][Ny + 2 * delta + 1].Ez = cube[i][1].Ez;
			//	cube[i][Ny + 2 * delta + 1].Bx = cube[i][1].Bx;
			//	cube[i][Ny + 2 * delta + 1].By = cube[i][1].By;
			//	cube[i][Ny + 2 * delta + 1].Bz = cube[i][1].Bz;
			//}
			vec_energy[it + 1] = CalculateEnergyPML_2dimen(cube, Nx, Ny, delta);
		}
	}
	cout << endl << vec_energy[Nt] / vec_energy[0] << endl;
	Graph_Solution_in_two_planes(cube, Nx, Ny, delta, dx, dy, direction);

	//����� ������� ������� � ����
	ofstream numb_energy("Energy_graph_with_PML.csv");
	for (int s = 0; s < Nt; s++)
	{
		numb_energy << dt * (double)(s + 1) << ";" << vec_energy[s] << endl;
	}
	numb_energy << endl << endl;
	numb_energy.close();

	return 0.0;
}

template <class ftype>
double CalculateEnergyPML_2dimen(vector<vector<Component<ftype>>>& cube, int Nx, int Ny, int delta)
{
	double energy = 0.0;

	for (int i =  1; i < Nx +  2 * delta + 1; i++)
		for (int j =  1; j < Ny + 2 * delta + 1; j++)
			{
				energy += cube[i][j].Ex * cube[i][j].Ex + cube[i][j].Ey * cube[i][j].Ey + cube[i][j].Ez * cube[i][j].Ez;
				energy += cube[i][j].Bx * cube[i][j].Bx + cube[i][j].By * cube[i][j].By + cube[i][j].Bz * cube[i][j].Bz;
			}
	return energy;
}

template <class ftype>
void Graph_for_Sigma(vector<vector<SIGMA<ftype>>>& Sigma, int Nx, int Ny, int delta, ftype dx, ftype dy, int fix_index)
{
	//����� ������� ���� � ����
	ofstream numb_sigma("Sigma_graph_1.csv");
	numb_sigma << "Graph for Ox" << endl;
	numb_sigma << "value X" << "value sigma_x" << endl;
	for (int s = 1; s < Nx + 2 * delta + 1; s++)
	{
		numb_sigma << dx * (ftype)(s - 1) << ";" << Sigma[s][fix_index].sigmaH_x << endl;
		numb_sigma << dx * ((ftype)(s - 1) + 0.5) << ";" << Sigma[s][fix_index].sigmaE_x << endl;
	}
	numb_sigma << endl << endl;

	numb_sigma << "Graph for Oy" << endl;
	numb_sigma << "value Y" << "value sigma_y" << endl;
	for (int s = 1; s < Ny + 2 * delta + 1; s++)
	{
		numb_sigma << dy * (ftype)(s - 1) << ";" << Sigma[fix_index][s].sigmaH_y << endl;
		numb_sigma << dy * ((ftype)(s - 1) + 0.5) << ";" << Sigma[fix_index][s].sigmaE_y << endl;
	}
	numb_sigma << endl << endl;
	numb_sigma.close();
}

template <class ftype>
void Graph_Solution_in_two_planes(vector<vector<Component<ftype>>>& cube, int Nx, int Ny, int delta, ftype dx, ftype dy, Direction direction)
{
	ofstream numbE("graph_solution_E_x.csv"), numbB("graph_solution_B_x.csv");
	numbE << direction << ";" << "y" << endl;
	numbE << "x" << ";" << ";";
	numbB << direction << ";" << "y" << endl;
	numbB << "x" << ";" << ";";
	for (int i = 1; i < Nx + 2 * delta + 1; i++)
	{
		numbE << dx * (ftype)i << ";";
		numbB << dx * (ftype)i << ";";
	}
	numbE << endl; 		numbB << endl;
	
	switch (direction)
	{
	case x_comp_yz:
	{
		for (int i = 1; i < Ny + 2 * delta + 1; i++)
		{
			ftype y = dy * (ftype)i;
			numbE << ";" << y << ";"; 			numbB << ";" << y << ";";
			for (int j = 1; j < Nx + 2 * delta + 1; j++)
			{
				numbE << cube[j][i].Ey << ";";
				numbB << cube[j][i].Bz << ";";
			}
			numbE << endl; 			numbB << endl;
		}
		break;
	}
	case y_comp_zx:
	{
		for (int i = 1; i < Ny + 2 * delta + 1; i++)
		{
			ftype y = dy * (ftype)i;
			numbE << ";" << y << ";"; 			numbB << ";" << y << ";";
			for (int j = 1; j < Nx + 2 * delta + 1; j++)
			{
				numbE << cube[j][i].Ez << ";";
				numbB << cube[j][i].Bx << ";";
			}
			numbE << endl; 			numbB << endl;
		}
		break;
	}
	}
	numbE.close();
	numbB.close();
}

template <class ftype>
void Graph_Solution_Coeff_in_two_planes(vector<vector<COEFF<ftype>>>& Coeff, int Nx, int Ny, int delta, ftype dx, ftype dy)
{
	ofstream numb2("graph_COeff.csv");
	numb2 << ";" << "y" << endl;
	numb2 << "x" << ";" << ";";
	for (int i = 1; i < Nx + 2 * delta + 1; i++)
	{
		numb2 << dx * (ftype)i << ";";
	}
	numb2 << endl;
	for (int i = 1; i < Ny + 2 * delta + 1; i++)
	{
		ftype y = dy * (ftype)i;
		numb2 << ";" << y << ";";
		for (int j = 1; j < Nx + 2 * delta + 1; j++)
		{
			numb2 << Coeff[j][i].Exy2 << ";";
		}
		numb2 << endl;
	}

	numb2.close();
}