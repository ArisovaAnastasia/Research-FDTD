#define _USE_MATH_DEFINES
#include <algorithm>
#include <iostream>
#include <fstream>
#include <clocale>
#include <chrono>
#include <vector>
#include <cmath>
#include <omp.h>
#include "half.hpp"

using namespace std;


template <class ftype>
void Initializing_FDTD_three_dimen_task1(vector<vector<vector<Component<ftype>>>>& cube, int Nx, int Ny, int Nz, ftype dx, const ftype dt,
	pair <ftype, ftype> ab, double& abs_err_t, double& rel_err_t, double& smape_err_t)
{
	ftype x1, x2, t;
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			for (int k = 0; k < Nz; k++)
			{
				x1 = (ftype)(ab.first + dx * (ftype)i + (ftype)0.5 * dx);
				x2 = (ftype)(ab.first + dx * (ftype)i);
				t = (ftype)0.5 * dt;

				cube[i][j][k].Ex = (ftype)0.0;
				cube[i][j][k].Ey = sin(x1);
				cube[i][j][k].Ez = (ftype)0.0;
				cube[i][j][k].Bx = (ftype)0.0;
				cube[i][j][k].By = (ftype)0.0;
				cube[i][j][k].Bz = sin(x2 - t);

				Check_error_current_step_electric_task1(cube, dt, dx, ab.first, 0, i, j, k, abs_err_t, rel_err_t, smape_err_t);
				Check_error_current_step_magnetic_task1(cube, dt, dx, ab.first, 0, i, j, k, abs_err_t, rel_err_t, smape_err_t);
			}

	//Заполняем отдельно границу куба
	for (int i = 0; i < Ny + 2; i++)
		for (int j = 0; j < Nz + 2; j++)
		{
			cube[Nx][i][j].Ex = cube[0][i][j].Ex;
			cube[Nx][i][j].Ey = cube[0][i][j].Ey;
			cube[Nx][i][j].Ez = cube[0][i][j].Ez;
			cube[Nx][i][j].Bx = cube[0][i][j].Bx;
			cube[Nx][i][j].By = cube[0][i][j].By;
			cube[Nx][i][j].Bz = cube[0][i][j].Bz;

			cube[Nx + 1][i][j].Ex = cube[1][i][j].Ex;
			cube[Nx + 1][i][j].Ey = cube[1][i][j].Ey;
			cube[Nx + 1][i][j].Ez = cube[1][i][j].Ez;
			cube[Nx + 1][i][j].Bx = cube[1][i][j].Bx;
			cube[Nx + 1][i][j].By = cube[1][i][j].By;
			cube[Nx + 1][i][j].Bz = cube[1][i][j].Bz;
		}
	for (int i = 0; i < Nx + 2; i++)
		for (int j = 0; j < Nz + 2; j++)
		{
			cube[i][Ny][j].Ex = cube[i][0][j].Ex;
			cube[i][Ny][j].Ey = cube[i][0][j].Ey;
			cube[i][Ny][j].Ez = cube[i][0][j].Ez;
			cube[i][Ny][j].Bx = cube[i][0][j].Bx;
			cube[i][Ny][j].By = cube[i][0][j].By;
			cube[i][Ny][j].Bz = cube[i][0][j].Bz;

			cube[i][Ny + 1][j].Ex = cube[i][1][j].Ex;
			cube[i][Ny + 1][j].Ey = cube[i][1][j].Ey;
			cube[i][Ny + 1][j].Ez = cube[i][1][j].Ez;
			cube[i][Ny + 1][j].Bx = cube[i][1][j].Bx;
			cube[i][Ny + 1][j].By = cube[i][1][j].By;
			cube[i][Ny + 1][j].Bz = cube[i][1][j].Bz;

		}
	for (int i = 0; i < Nx + 2; i++)
		for (int j = 0; j < Ny + 2; j++)
		{
			cube[i][j][Nz].Ex = cube[i][j][0].Ex;
			cube[i][j][Nz].Ey = cube[i][j][0].Ey;
			cube[i][j][Nz].Ez = cube[i][j][0].Ez;
			cube[i][j][Nz].Bx = cube[i][j][0].Bx;
			cube[i][j][Nz].By = cube[i][j][0].By;
			cube[i][j][Nz].Bz = cube[i][j][0].Bz;

			cube[i][j][Nz + 1].Ex = cube[i][j][1].Ex;
			cube[i][j][Nz + 1].Ey = cube[i][j][1].Ey;
			cube[i][j][Nz + 1].Ez = cube[i][j][1].Ez;
			cube[i][j][Nz + 1].Bx = cube[i][j][1].Bx;
			cube[i][j][Nz + 1].By = cube[i][j][1].By;
			cube[i][j][Nz + 1].Bz = cube[i][j][1].Bz;
		}
} 

template <class ftype>
void Initializing_FDTD_three_dimen_task2(vector<vector<vector<Component<ftype>>>>& cube, int Nx, int Ny, int Nz, ftype dy, const ftype dt,
	pair <ftype, ftype> cd, double& abs_err_t, double& rel_err_t, double& smape_err_t)
{
	ftype y1, y2, t;
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			for (int k = 0; k < Nz; k++)
			{
				y1 = (ftype)(cd.first + dy * (ftype)j + (ftype)0.5 * dy);
				y2 = (ftype)(cd.first + dy * (ftype)j);
				t = (ftype)0.5 * dt;

				cube[i][j][k].Ex = (ftype)0.0;
				cube[i][j][k].Ey = (ftype)0.0;
				cube[i][j][k].Ez = sin(y1);
				cube[i][j][k].Bx = sin(y2 - t);
				cube[i][j][k].By = (ftype)0.0;
				cube[i][j][k].Bz = (ftype)0.0;
				Check_error_current_step_electric_task2(cube, dt, dy, cd.first, 0, i, j, k, abs_err_t, rel_err_t, smape_err_t);
				Check_error_current_step_magnetic_task2(cube, dt, dy, cd.first, 0, i, j, k, abs_err_t, rel_err_t, smape_err_t);
			}

	//Заполняем отдельно границу куба
	for (int i = 0; i < Ny + 2; i++)
		for (int j = 0; j < Nz + 2; j++)
		{
			cube[Nx][i][j].Ex = cube[0][i][j].Ex;
			cube[Nx][i][j].Ey = cube[0][i][j].Ey;
			cube[Nx][i][j].Ez = cube[0][i][j].Ez;
			cube[Nx][i][j].Bx = cube[0][i][j].Bx;
			cube[Nx][i][j].By = cube[0][i][j].By;
			cube[Nx][i][j].Bz = cube[0][i][j].Bz;

			cube[Nx + 1][i][j].Ex = cube[1][i][j].Ex;
			cube[Nx + 1][i][j].Ey = cube[1][i][j].Ey;
			cube[Nx + 1][i][j].Ez = cube[1][i][j].Ez;
			cube[Nx + 1][i][j].Bx = cube[1][i][j].Bx;
			cube[Nx + 1][i][j].By = cube[1][i][j].By;
			cube[Nx + 1][i][j].Bz = cube[1][i][j].Bz;
		}
	for (int i = 0; i < Nx + 2; i++)
		for (int j = 0; j < Nz + 2; j++)
		{
			cube[i][Ny][j].Ex = cube[i][0][j].Ex;
			cube[i][Ny][j].Ey = cube[i][0][j].Ey;
			cube[i][Ny][j].Ez = cube[i][0][j].Ez;
			cube[i][Ny][j].Bx = cube[i][0][j].Bx;
			cube[i][Ny][j].By = cube[i][0][j].By;
			cube[i][Ny][j].Bz = cube[i][0][j].Bz;

			cube[i][Ny + 1][j].Ex = cube[i][1][j].Ex;
			cube[i][Ny + 1][j].Ey = cube[i][1][j].Ey;
			cube[i][Ny + 1][j].Ez = cube[i][1][j].Ez;
			cube[i][Ny + 1][j].Bx = cube[i][1][j].Bx;
			cube[i][Ny + 1][j].By = cube[i][1][j].By;
			cube[i][Ny + 1][j].Bz = cube[i][1][j].Bz;

		}
	for (int i = 0; i < Nx + 2; i++)
		for (int j = 0; j < Ny + 2; j++)
		{
			cube[i][j][Nz].Ex = cube[i][j][0].Ex;
			cube[i][j][Nz].Ey = cube[i][j][0].Ey;
			cube[i][j][Nz].Ez = cube[i][j][0].Ez;
			cube[i][j][Nz].Bx = cube[i][j][0].Bx;
			cube[i][j][Nz].By = cube[i][j][0].By;
			cube[i][j][Nz].Bz = cube[i][j][0].Bz;

			cube[i][j][Nz + 1].Ex = cube[i][j][1].Ex;
			cube[i][j][Nz + 1].Ey = cube[i][j][1].Ey;
			cube[i][j][Nz + 1].Ez = cube[i][j][1].Ez;
			cube[i][j][Nz + 1].Bx = cube[i][j][1].Bx;
			cube[i][j][Nz + 1].By = cube[i][j][1].By;
			cube[i][j][Nz + 1].Bz = cube[i][j][1].Bz;
		}

}

template <class ftype>
void Initializing_FDTD_three_dimen_task3(vector<vector<vector<Component<ftype>>>>& cube, int Nx, int Ny, int Nz, ftype dz, const ftype dt,
	pair <ftype, ftype> fg, double& abs_err_t, double& rel_err_t, double& smape_err_t)
{
	ftype z1, z2, t;
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			for (int k = 0; k < Nz; k++)
			{
				z1 = (ftype)(fg.first + dz * (ftype)k + (ftype)0.5 * dz);
				z2 = (ftype)(fg.first + dz * (ftype)k);
				t = (ftype)0.5 * dt;

				cube[i][j][k].Ex = sin(z1);
				cube[i][j][k].Ey = (ftype)0.0;
				cube[i][j][k].Ez = (ftype)0.0;
				cube[i][j][k].Bx = (ftype)0.0;
				cube[i][j][k].By = sin(z2 - t);
				cube[i][j][k].Bz = (ftype)0.0;
				Check_error_current_step_electric_task3(cube, dt, dz, fg.first, 0, i, j, k, abs_err_t, rel_err_t, smape_err_t);
				Check_error_current_step_magnetic_task3(cube, dt, dz, fg.first, 0, i, j, k, abs_err_t, rel_err_t, smape_err_t);
			}

	//Заполняем отдельно границу куба
	for (int i = 0; i < Ny + 2; i++)
		for (int j = 0; j < Nz + 2; j++)
		{
			cube[Nx][i][j].Ex = cube[0][i][j].Ex;
			cube[Nx][i][j].Ey = cube[0][i][j].Ey;
			cube[Nx][i][j].Ez = cube[0][i][j].Ez;
			cube[Nx][i][j].Bx = cube[0][i][j].Bx;
			cube[Nx][i][j].By = cube[0][i][j].By;
			cube[Nx][i][j].Bz = cube[0][i][j].Bz;

			cube[Nx + 1][i][j].Ex = cube[1][i][j].Ex;
			cube[Nx + 1][i][j].Ey = cube[1][i][j].Ey;
			cube[Nx + 1][i][j].Ez = cube[1][i][j].Ez;
			cube[Nx + 1][i][j].Bx = cube[1][i][j].Bx;
			cube[Nx + 1][i][j].By = cube[1][i][j].By;
			cube[Nx + 1][i][j].Bz = cube[1][i][j].Bz;
		}
	for (int i = 0; i < Nx + 2; i++)
		for (int j = 0; j < Nz + 2; j++)
		{
			cube[i][Ny][j].Ex = cube[i][0][j].Ex;
			cube[i][Ny][j].Ey = cube[i][0][j].Ey;
			cube[i][Ny][j].Ez = cube[i][0][j].Ez;
			cube[i][Ny][j].Bx = cube[i][0][j].Bx;
			cube[i][Ny][j].By = cube[i][0][j].By;
			cube[i][Ny][j].Bz = cube[i][0][j].Bz;

			cube[i][Ny + 1][j].Ex = cube[i][1][j].Ex;
			cube[i][Ny + 1][j].Ey = cube[i][1][j].Ey;
			cube[i][Ny + 1][j].Ez = cube[i][1][j].Ez;
			cube[i][Ny + 1][j].Bx = cube[i][1][j].Bx;
			cube[i][Ny + 1][j].By = cube[i][1][j].By;
			cube[i][Ny + 1][j].Bz = cube[i][1][j].Bz;

		}
	for (int i = 0; i < Nx + 2; i++)
		for (int j = 0; j < Ny + 2; j++)
		{
			cube[i][j][Nz].Ex = cube[i][j][0].Ex;
			cube[i][j][Nz].Ey = cube[i][j][0].Ey;
			cube[i][j][Nz].Ez = cube[i][j][0].Ez;
			cube[i][j][Nz].Bx = cube[i][j][0].Bx;
			cube[i][j][Nz].By = cube[i][j][0].By;
			cube[i][j][Nz].Bz = cube[i][j][0].Bz;

			cube[i][j][Nz + 1].Ex = cube[i][j][1].Ex;
			cube[i][j][Nz + 1].Ey = cube[i][j][1].Ey;
			cube[i][j][Nz + 1].Ez = cube[i][j][1].Ez;
			cube[i][j][Nz + 1].Bx = cube[i][j][1].Bx;
			cube[i][j][Nz + 1].By = cube[i][j][1].By;
			cube[i][j][Nz + 1].Bz = cube[i][j][1].Bz;
		}

}

//template <class ftype>
//void Update_electric_field_three_dimen(vector<vector<vector<Component<ftype>>>>& cube, vector<vector<vector<Component<ftype>>>>& cube2,
//	ftype dt_x, ftype dt_y, ftype dt_z, ftype dt, ftype dx, ftype ab_first, int it,
//	int i, int j, int k)
//{
//	ftype tEx, tEy, tEz;
//
//	tEx = (cube[i][j][k].Bz - cube[i][j - 1][k].Bz) * (dt_y)
//		- (cube[i][j][k].By - cube[i][j][k - 1].By) * (dt_z) + cube[i][j][k].Ex;
//	tEy = (cube[i][j][k].Bx - cube[i][j][k - 1].Bx) * (dt_z)
//		- (cube[i][j][k].Bz - cube[i - 1][j][k].Bz) * (dt_x) + cube[i][j][k].Ey;
//	tEz = (cube[i][j][k].By - cube[i - 1][j][k].By) * (dt_x)
//		- (cube[i][j][k].Bx - cube[i][j - 1][k].Bx) * (dt_y) + cube[i][j][k].Ez;
//
//	cube2[i][j][k].Ex = tEx;
//	cube2[i][j][k].Ey = tEy;
//	cube2[i][j][k].Ez = tEz;
//}
template <class ftype>
void Update_electric_field_three_dimen(vector<vector<vector<Component<ftype>>>>& cube,
	ftype dt_x, ftype dt_y, ftype dt_z, int i, int j, int k)
{
	ftype tEx, tEy, tEz;

	tEx = (cube[i][j + 1][k].Bz - cube[i][j][k].Bz) * (dt_y)
		-(cube[i][j][k + 1].By - cube[i][j][k].By) * (dt_z)+cube[i][j][k].Ex;
	tEy = (cube[i][j][k + 1].Bx - cube[i][j][k].Bx) * (dt_z)
		-(cube[i + 1][j][k].Bz - cube[i][j][k].Bz) * (dt_x)+cube[i][j][k].Ey;
	tEz = (cube[i + 1][j][k].By - cube[i][j][k].By) * (dt_x)
		-(cube[i][j + 1][k].Bx - cube[i][j][k].Bx) * (dt_y)+cube[i][j][k].Ez;

	cube[i][j][k].Ex = tEx;
	cube[i][j][k].Ey = tEy;
	cube[i][j][k].Ez = tEz;
}

//template <class ftype>
//void Update_magnetic_field_three_dimen(vector<vector<vector<Component<ftype>>>>& cube, vector<vector<vector<Component<ftype>>>>& cube2,
//	ftype dt_x, ftype dt_y, ftype dt_z, ftype dt, ftype dx, ftype ab_first, int it,
//	int i, int j, int k)
//{
//	ftype tBx, tBy, tBz;
//
//	tBx = (cube[i][j][k].Ez - cube[i][j + 1][k].Ez) * (dt_y)
//		+(cube[i][j][k + 1].Ey - cube[i][j][k].Ey) * (dt_z)+cube[i][j][k].Bx;
//	tBy = (cube[i][j][k].Ex - cube[i][j][k + 1].Ex) * (dt_z)
//		+(cube[i + 1][j][k].Ez - cube[i][j][k].Ez) * (dt_x)+cube[i][j][k].By;
//	tBz = (cube[i][j][k].Ey - cube[i + 1][j][k].Ey) * (dt_x)
//		+(cube[i][j + 1][k].Ex - cube[i][j][k].Ex) * (dt_y)+cube[i][j][k].Bz;
//
//	cube2[i][j][k].Bx = tBx;
//	cube2[i][j][k].By = tBy;
//	cube2[i][j][k].Bz = tBz;
//}
template <class ftype>
void Update_magnetic_field_three_dimen(vector<vector<vector<Component<ftype>>>>& cube,
	ftype dt_x, ftype dt_y, ftype dt_z,	int i, int j, int k)
{
	ftype tBx, tBy, tBz;

	tBx = -(cube[i][j][k].Ez - cube[i][j - 1][k].Ez) * (dt_y)
		+ (cube[i][j][k].Ey - cube[i][j][k - 1].Ey) * (dt_z) + cube[i][j][k].Bx;

	tBy = -(cube[i][j][k].Ex - cube[i][j][k - 1].Ex) * (dt_z)
		+ (cube[i][j][k].Ez - cube[i - 1][j][k].Ez) * (dt_x) + cube[i][j][k].By;
	tBz = -(cube[i][j][k].Ey - cube[i - 1][j][k].Ey) * (dt_x)
		+ (cube[i][j][k].Ex - cube[i][j - 1][k].Ex) * (dt_y) + cube[i][j][k].Bz;

	cube[i][j][k].Bx = tBx;
	cube[i][j][k].By = tBy;
	cube[i][j][k].Bz = tBz;
}

//template <class ftype>
//void Update_electric_field_three_dimen_calc_error_task1(vector<vector<vector<Component<ftype>>>>& cube, vector<vector<vector<Component<ftype>>>>& cube2,
//	ftype dt_x, ftype dt_y, ftype dt_z, ftype dt, ftype dx, ftype ab_first, int it,
//	int i, int j, int k, double& abs_err_t, double& rel_err_t, double& smape_err_t)
//{
//	ftype tEx, tEy, tEz;
//
//	tEx = (cube[i][j + 1][k].Bz - cube[i][j][k].Bz) * (dt_y)
//		-(cube[i][j][k + 1].By - cube[i][j][k].By) * (dt_z)+cube[i][j][k].Ex;
//	tEy = (cube[i][j][k + 1].Bx - cube[i][j][k].Bx) * (dt_z)
//		-(cube[i + 1][j][k].Bz - cube[i][j][k].Bz) * (dt_x)+cube[i][j][k].Ey;
//	tEz = (cube[i + 1][j][k].By - cube[i][j][k].By) * (dt_x)
//		-(cube[i][j + 1][k].Bx - cube[i][j][k].Bx) * (dt_y)+cube[i][j][k].Ez;
//
//	cube2[i][j][k].Ex = tEx;
//	cube2[i][j][k].Ey = tEy;
//	cube2[i][j][k].Ez = tEz;
//
//	double x = (double)ab_first + (double)i * (double)dx + (ftype)0.5 * dx, t = (double)it * (double)dt;
//
//	//Считаем текущую абсолютную ошибку в данном узле
//	abs_err_t = max(abs_err_t, abs(((double)sin(x - t)) - (double)cube2[i][j][k].Ey));
//
//	//Считаем текущую относительную ошибку в данном узле
//	if (t != 0)
//	{
//		rel_err_t = max(rel_err_t, abs(sin(x - t) - (double)cube2[i][j][k].Ey) / abs(sin(x - t)));
//	}
//
//	//Считаем текущую ошибку smape в данном узле
//	double temp;
//	temp = (double)2 * abs((double)sin(x - t) - (double)cube2[i][j][k].Ey) /
//		(abs((double)sin(x - t)) + abs((double)cube2[i][j][k].Ey));
//	smape_err_t += temp;
//}
//
//template <class ftype>
//void Update_magnetic_field_three_dimen_calc_error_task1(vector<vector<vector<Component<ftype>>>>& cube, vector<vector<vector<Component<ftype>>>>& cube2,
//	ftype dt_x, ftype dt_y, ftype dt_z, ftype dt, ftype dx, ftype ab_first, int it,
//	int i, int j, int k, double& abs_err_t, double& rel_err_t, double& smape_err_t)
//{
//	ftype tBx, tBy, tBz;
//
//	tBx = (cube[i][j - 1][k].Ez - cube[i][j][k].Ez) * (dt_y)
//		+(cube[i][j][k].Ey - cube[i][j][k - 1].Ey) * (dt_z)+cube[i][j][k].Bx;
//	tBy = (cube[i][j][k - 1].Ex - cube[i][j][k].Ex) * (dt_z)
//		+(cube[i][j][k].Ez - cube[i - 1][j][k].Ez) * (dt_x)+cube[i][j][k].By;
//	tBz = (cube[i - 1][j][k].Ey - cube[i][j][k].Ey) * (dt_x)
//		+(cube[i][j][k].Ex - cube[i][j - 1][k].Ex) * (dt_y)+cube[i][j][k].Bz;
//
//	cube2[i][j][k].Bx = tBx;
//	cube2[i][j][k].By = tBy;
//	cube2[i][j][k].Bz = tBz;
//
//	double x = ab_first + (double)i * dx, t = (double)it * dt + (double)0.5 * dt;
//	
//	//Считаем текущую абсолютную ошибку в данном узле
//	abs_err_t = max(abs_err_t, abs(((double)sin(x - t)) - (double)cube2[i][j][k].Bz));
//
//	//Считаем текущую относительную ошибку в данном узле
//	if (t != 0)
//	{
//		rel_err_t = max(rel_err_t, abs(sin(x - t) - (double)cube2[i][j][k].Bz) / abs(sin(x - t)));
//	}
//
//	//Считаем текущую ошибку smape в данном узле
//	double temp;
//	temp = (double)2 * abs((double)sin(x - t) - (double)cube2[i][j][k].Bz) /
//		(abs((double)sin(x - t)) + abs((double)cube2[i][j][k].Bz));
//	smape_err_t += temp;
//}
//
//template <class ftype>
//void Update_electric_field_three_dimen_calc_error_task2(vector<vector<vector<Component<ftype>>>>& cube, vector<vector<vector<Component<ftype>>>>& cube2,
//	ftype dt_x, ftype dt_y, ftype dt_z, ftype dt, ftype dy, ftype cd_first, int it,
//	int i, int j, int k, double& abs_err_t, double& rel_err_t, double& smape_err_t)
//{
//	ftype tEx, tEy, tEz;
//
//	tEx = (cube[i][j + 1][k].Bz - cube[i][j][k].Bz) * (dt_y)
//		-(cube[i][j][k + 1].By - cube[i][j][k].By) * (dt_z)+cube[i][j][k].Ex;
//	tEy = (cube[i][j][k + 1].Bx - cube[i][j][k].Bx) * (dt_z)
//		-(cube[i + 1][j][k].Bz - cube[i][j][k].Bz) * (dt_x)+cube[i][j][k].Ey;
//	tEz = (cube[i + 1][j][k].By - cube[i][j][k].By) * (dt_x)
//		-(cube[i][j + 1][k].Bx - cube[i][j][k].Bx) * (dt_y)+cube[i][j][k].Ez;
//
//	cube2[i][j][k].Ex = tEx;
//	cube2[i][j][k].Ey = tEy;
//	cube2[i][j][k].Ez = tEz;
//}
//
//template <class ftype>
//void Update_magnetic_field_three_dimen_calc_error_task2(vector<vector<vector<Component<ftype>>>>& cube, vector<vector<vector<Component<ftype>>>>& cube2,
//	ftype dt_x, ftype dt_y, ftype dt_z, ftype dt, ftype dy, ftype cd_first, int it,
//	int i, int j, int k, double& abs_err_t, double& rel_err_t, double& smape_err_t)
//{
//	ftype tBx, tBy, tBz;
//
//	tBx = (cube[i][j - 1][k].Ez - cube[i][j][k].Ez) * (dt_y)
//		+(cube[i][j][k].Ey - cube[i][j][k - 1].Ey) * (dt_z)+cube[i][j][k].Bx;
//	tBy = (cube[i][j][k - 1].Ex - cube[i][j][k].Ex) * (dt_z)
//		+(cube[i][j][k].Ez - cube[i - 1][j][k].Ez) * (dt_x)+cube[i][j][k].By;
//	tBz = (cube[i - 1][j][k].Ey - cube[i][j][k].Ey) * (dt_x)
//		+(cube[i][j][k].Ex - cube[i][j - 1][k].Ex) * (dt_y)+cube[i][j][k].Bz;
//
//	cube2[i][j][k].Bx = tBx;
//	cube2[i][j][k].By = tBy;
//	cube2[i][j][k].Bz = tBz;
//
//	double y = cd_first + (double)j * dy, t = (double)it * dt + (double)0.5 * dt;
//
//	//Считаем текущую абсолютную ошибку в данном узле
//	abs_err_t = max(abs_err_t, abs(((double)sin(y - t)) - (double)cube2[i][j][k].Bx));
//
//	//Считаем текущую относительную ошибку в данном узле
//	if (t != 0)
//	{
//		rel_err_t = max(rel_err_t, abs(sin(y - t) - (double)cube2[i][j][k].Bx) / abs(sin(y - t)));
//	}
//
//	//Считаем текущую ошибку smape в данном узле
//	double temp;
//	temp = (double)2 * abs((double)sin(y - t) - (double)cube2[i][j][k].Bx) /
//		(abs((double)sin(y - t)) + abs((double)cube2[i][j][k].Bx));
//	smape_err_t += temp;
//}

template <class ftype>
void Check_error_current_step_electric_task1(vector<vector<vector<Component<ftype>>>>& cube2, ftype dt, ftype dx, ftype ab_first, int it,
	int i, int j, int k, double& abs_err_t, double& rel_err_t, double& smape_err_t)
{
	double x1 = (double)(ab_first + dx * (double)i + (double)0.5 * dx), t1 = (double)it * dt;

	//Считаем текущую абсолютную ошибку в данном узле
	abs_err_t = max(abs_err_t, abs(((double)sin((x1 - t1))) - (double)cube2[i][j][k].Ey));

	////Считаем текущую относительную ошибку в данном узле
	//if (t != 0)
	//{
	//	rel_err_t = max(rel_err_t, abs(sin(y2 - t) - (double)cube2[i][j][k].Bx) / abs(sin(y - t)));
	//	rel_err_t = max(rel_err_t, abs(sin(y2 - t) - (double)cube2[i][j][k].Bx) / abs(sin(y - t)));
	//}

	//Считаем текущую ошибку smape в данном узле
	double temp;
	temp = (double)2 * abs((double)sin(x1 - t1) - (double)cube2[i][j][k].Ey) /
		(abs((double)sin(x1 - t1)) + abs((double)cube2[i][j][k].Ey));
	smape_err_t += temp;
}

template <class ftype>
void Check_error_current_step_magnetic_task1(vector<vector<vector<Component<ftype>>>>& cube2, ftype dt, ftype dx, ftype ab_first, int it,
	int i, int j, int k, double& abs_err_t, double& rel_err_t, double& smape_err_t)
{
	double x2 = (double)(ab_first + dx * (double)i), t2 = (double)it * dt + (double)0.5 * dt;

	//Считаем текущую абсолютную ошибку в данном узле
	abs_err_t = max(abs_err_t, abs(((double)sin(x2 - t2)) - (double)cube2[i][j][k].Bz));

	////Считаем текущую относительную ошибку в данном узле
	//if (t != 0)
	//{
	//	rel_err_t = max(rel_err_t, abs(sin(y2 - t) - (double)cube2[i][j][k].Bx) / abs(sin(y - t)));
	//	rel_err_t = max(rel_err_t, abs(sin(y2 - t) - (double)cube2[i][j][k].Bx) / abs(sin(y - t)));
	//}

	//Считаем текущую ошибку smape в данном узле
	double temp;
	temp = (double)2 * abs((double)sin(x2 - t2) - (double)cube2[i][j][k].Bz) /
		(abs((double)sin(x2 - t2)) + abs((double)cube2[i][j][k].Bz));
	smape_err_t += temp;
}

template <class ftype>
void Check_error_current_step_electric_task2(vector<vector<vector<Component<ftype>>>>& cube2, ftype dt, ftype dy, ftype cd_first, int it,
	int i, int j, int k, double& abs_err_t, double& rel_err_t, double& smape_err_t)
{
	double y1 = (double)(cd_first + dy * (double)j + (double)0.5 * dy), t1 = (double)it * dt;

	//Считаем текущую абсолютную ошибку в данном узле
	abs_err_t = max(abs_err_t, abs(((double)sin(y1 - t1)) - (double)cube2[i][j][k].Ez));

	////Считаем текущую относительную ошибку в данном узле
	//if (t != 0)
	//{
	//	rel_err_t = max(rel_err_t, abs(sin(y2 - t) - (double)cube2[i][j][k].Bx) / abs(sin(y - t)));
	//	rel_err_t = max(rel_err_t, abs(sin(y2 - t) - (double)cube2[i][j][k].Bx) / abs(sin(y - t)));
	//}

	////Считаем текущую ошибку smape в данном узле
	//double temp;
	//temp = (double)2 * abs((double)sin(y2 - t) - (double)cube2[i][j][k].Bx) /
	//	(abs((double)sin(y2 - t)) + abs((double)cube2[i][j][k].Bx));
	//smape_err_t += temp;
}

template <class ftype>
void Check_error_current_step_magnetic_task2(vector<vector<vector<Component<ftype>>>>& cube2, ftype dt, ftype dy, ftype cd_first, int it,
	int i, int j, int k, double& abs_err_t, double& rel_err_t, double& smape_err_t)
{
	double y2 = (double)(cd_first + dy * (double)j), t2 = (double)it * dt + (double)0.5 * dt;

	//Считаем текущую абсолютную ошибку в данном узле
	abs_err_t = max(abs_err_t, abs(((double)sin(y2 - t2)) - (double)cube2[i][j][k].Bx));

	////Считаем текущую относительную ошибку в данном узле
	//if (t != 0)
	//{
	//	rel_err_t = max(rel_err_t, abs(sin(y2 - t) - (double)cube2[i][j][k].Bx) / abs(sin(y - t)));
	//	rel_err_t = max(rel_err_t, abs(sin(y2 - t) - (double)cube2[i][j][k].Bx) / abs(sin(y - t)));
	//}

	////Считаем текущую ошибку smape в данном узле
	//double temp;
	//temp = (double)2 * abs((double)sin(y2 - t) - (double)cube2[i][j][k].Bx) /
	//	(abs((double)sin(y2 - t)) + abs((double)cube2[i][j][k].Bx));
	//smape_err_t += temp;
}

template <class ftype>
void Check_error_current_step_electric_task3(vector<vector<vector<Component<ftype>>>>& cube2, ftype dt, ftype dz, ftype fg_first, int it,
	int i, int j, int k, double& abs_err_t, double& rel_err_t, double& smape_err_t)
{
	double z1 = (double)(fg_first + dz * (double)k + (double)0.5 * dz), t1 = (double)it * dt;

	//Считаем текущую абсолютную ошибку в данном узле
	abs_err_t = max(abs_err_t, abs(((double)sin(z1 - t1)) - (double)cube2[i][j][k].Ex));

	////Считаем текущую относительную ошибку в данном узле
	//if (t != 0)
	//{
	//	rel_err_t = max(rel_err_t, abs(sin(y2 - t) - (double)cube2[i][j][k].Bx) / abs(sin(y - t)));
	//	rel_err_t = max(rel_err_t, abs(sin(y2 - t) - (double)cube2[i][j][k].Bx) / abs(sin(y - t)));
	//}

	////Считаем текущую ошибку smape в данном узле
	//double temp;
	//temp = (double)2 * abs((double)sin(y2 - t) - (double)cube2[i][j][k].Bx) /
	//	(abs((double)sin(y2 - t)) + abs((double)cube2[i][j][k].Bx));
	//smape_err_t += temp;
}

template <class ftype>
void Check_error_current_step_magnetic_task3(vector<vector<vector<Component<ftype>>>>& cube2, ftype dt, ftype dz, ftype fg_first, int it,
	int i, int j, int k, double& abs_err_t, double& rel_err_t, double& smape_err_t)
{
	double z2 = (double)(fg_first + dz * (double)k), t2 = (double)it * dt + (double)0.5 * dt;

	//Считаем текущую абсолютную ошибку в данном узле
	abs_err_t = max(abs_err_t, abs(((double)sin(z2 - t2)) - (double)cube2[i][j][k].By));

	////Считаем текущую относительную ошибку в данном узле
	//if (t != 0)
	//{
	//	rel_err_t = max(rel_err_t, abs(sin(y2 - t) - (double)cube2[i][j][k].Bx) / abs(sin(y - t)));
	//	rel_err_t = max(rel_err_t, abs(sin(y2 - t) - (double)cube2[i][j][k].Bx) / abs(sin(y - t)));
	//}

	////Считаем текущую ошибку smape в данном узле
	//double temp;
	//temp = (double)2 * abs((double)sin(y2 - t) - (double)cube2[i][j][k].Bx) /
	//	(abs((double)sin(y2 - t)) + abs((double)cube2[i][j][k].Bx));
	//smape_err_t += temp;
}
//template <class ftype>
//double FDTD_three_dimen(int Nx, int Ny, int Nz, ftype T, ftype dt, string type_sum)
//{
//	setlocale(LC_ALL, "Russian");
//	int Nt = T / dt;
//
//	vector<vector<vector<Component<ftype>>>> cube(Nx + 2, vector<vector<Component<ftype>>>(Ny + 2, vector<Component<ftype>>(Nz + 2)));
//	vector<vector<vector<Component<ftype>>>> cube2(Nx + 2, vector<vector<Component<ftype>>>(Ny + 2, vector<Component<ftype>>(Nz + 2)));
//	/*vector<vector<vector<Component<ftype>>>> compensator(N_3 + 2, vector<vector<Component<ftype>>>(N_3 + 2, vector<Component<ftype>>(N_3 + 2)));
//	vector<vector<vector<Component<ftype>>>> compensator2(N_3 + 2, vector<vector<Component<ftype>>>(N_3 + 2, vector<Component<ftype>>(N_3 + 2)));*/
//
//	pair<ftype, ftype> ab(0.0, 2.0 * M_PI), cd(0.0, 2.0 * M_PI), fg(0.0, 2.0 * M_PI);
//	ftype dx = (ab.second - ab.first) / (ftype)Nx, dy = (cd.second - cd.first) / (ftype)Ny,
//		dz = (fg.second - fg.first) / (ftype)Nz;
//
//	ftype dt_x = dt / (dx);
//	ftype dt_y = dt / (dy);
//	ftype dt_z = dt / (dz);
//
//	Initializing_FDTD_three_dimen<ftype>(cube, Nx, Ny, Nz, dx, dt, ab);
//	cout <<endl << "Ey = "<<cube[0][0][0].Ey << "  Bz = " << cube[0][0][0].Bz << endl;
//	if (type_sum == "Kahan")
//	{
//	}
//	else
//	{
//		for (int it = 0; it < Nt; it++)
//		{
//			if (it % 100 == 0)
//				cout << it << endl;
//
//#pragma omp parallel for collapse(3)
//			for (int i = 0; i < Nx + 1; i++)
//				for (int j = 0; j < Ny + 1; j++)
//					for (int k = 0; k < Nz + 1; k++)
//						Update_electric_field_three_dimen<ftype>(cube, cube2, dt_x, dt_y, dt_z, dt, dx, ab.first, it, i, j, k);
//
//			for (int i = 0; i < Nx + 1; i++)
//				for (int j = 0; j < Ny + 1; j++)
//					for (int k = 0; k < Nz + 1; k++)
//					{
//						cube[i][j][k].Ex = cube2[i][j][k].Ex;
//						cube[i][j][k].Ey = cube2[i][j][k].Ey;
//						cube[i][j][k].Ez = cube2[i][j][k].Ez;
//					}
//
//#pragma omp parallel for collapse(3)
//			for (int i = 1; i < Nx + 1; i++)
//				for (int j = 1; j < Ny + 1; j++)
//					for (int k = 1; k < Nz + 1; k++)
//						Update_magnetic_field_three_dimen<ftype>(cube, cube2, dt_x, dt_y, dt_z, dt, dx, ab.first, it, i, j, k);
//
//			for (int i = 1; i < Nx + 1; i++)
//				for (int j = 1; j < Ny + 1; j++)
//					for (int k = 1; k < Nz + 1; k++)
//					{
//						cube[i][j][k].Bx = cube2[i][j][k].Bx;
//						cube[i][j][k].By = cube2[i][j][k].By;
//						cube[i][j][k].Bz = cube2[i][j][k].Bz;
//					}
//
//			// Изменяем значения на границе
//			for (int i = 0; i < Ny + 2; i++)
//				for (int j = 0; j < Nz + 2; j++)
//				{
//					cube[Nx + 1][i][j].Ex = cube[1][i][j].Ex;
//					cube[Nx + 1][i][j].Ey = cube[1][i][j].Ey;
//					cube[Nx + 1][i][j].Ez = cube[1][i][j].Ez;
//					cube[Nx + 1][i][j].Bx = cube[1][i][j].Bx;
//					cube[Nx + 1][i][j].By = cube[1][i][j].By;
//					cube[Nx + 1][i][j].Bz = cube[1][i][j].Bz;
//
//					cube[0][i][j].Ex = cube[Nx][i][j].Ex;
//					cube[0][i][j].Ey = cube[Nx][i][j].Ey;
//					cube[0][i][j].Ez = cube[Nx][i][j].Ez;
//					cube[0][i][j].Bx = cube[Nx][i][j].Bx;
//					cube[0][i][j].By = cube[Nx][i][j].By;
//					cube[0][i][j].Bz = cube[Nx][i][j].Bz;
//				}
//			for (int i = 0; i < Nx + 2; i++)
//				for (int j = 0; j < Nz + 2; j++)
//				{
//					cube[i][Ny + 1][j].Ex = cube[i][1][j].Ex;
//					cube[i][Ny + 1][j].Ey = cube[i][1][j].Ey;
//					cube[i][Ny + 1][j].Ez = cube[i][1][j].Ez;
//					cube[i][Ny + 1][j].Bx = cube[i][1][j].Bx;
//					cube[i][Ny + 1][j].By = cube[i][1][j].By;
//					cube[i][Ny + 1][j].Bz = cube[i][1][j].Bz;
//
//					cube[i][0][j].Ex = cube[i][Ny][j].Ex;
//					cube[i][0][j].Ey = cube[i][Ny][j].Ey;
//					cube[i][0][j].Ez = cube[i][Ny][j].Ez;
//					cube[i][0][j].Bx = cube[i][Ny][j].Bx;
//					cube[i][0][j].By = cube[i][Ny][j].By;
//					cube[i][0][j].Bz = cube[i][Ny][j].Bz;
//				}
//			for (int i = 0; i < Nx + 2; i++)
//				for (int j = 0; j < Ny + 2; j++)
//				{
//					cube[i][j][Nz + 1].Ex = cube[i][j][1].Ex;
//					cube[i][j][Nz + 1].Ey = cube[i][j][1].Ey;
//					cube[i][j][Nz + 1].Ez = cube[i][j][1].Ez;
//					cube[i][j][Nz + 1].Bx = cube[i][j][1].Bx;
//					cube[i][j][Nz + 1].By = cube[i][j][1].By;
//					cube[i][j][Nz + 1].Bz = cube[i][j][1].Bz;
//
//					cube[i][j][0].Ex = cube[i][j][Nz].Ex;
//					cube[i][j][0].Ey = cube[i][j][Nz].Ey;
//					cube[i][j][0].Ez = cube[i][j][Nz].Ez;
//					cube[i][j][0].Bx = cube[i][j][Nz].Bx;
//					cube[i][j][0].By = cube[i][j][Nz].By;
//					cube[i][j][0].Bz = cube[i][j][Nz].Bz;
//				}
//		}
//	}
//
//	vector<vector<vector<Component<ftype>>>> analitical_cube(Nx + 1, vector<vector<Component<ftype>>>(Ny + 1, vector<Component<ftype>>(Nz + 1)));
//
//	double t1 = (double)Nt * dt;
//	double t2 = (double)Nt * dt + 0.5 * dt;
//	double x1, x2;
//
//	for (int i = 0; i < Nx + 1; i++)
//		for (int j = 0; j < Ny + 1; j++)
//			for (int k = 0; k < Nz + 1; k++)
//			{
//				x1 = ab.first + dx * (double)i + 0.5 * dx;
//				x2 = ab.first + dx * (double)i;
//
//				analitical_cube[i][j][k].Ey = sin(x1 - t1);
//				analitical_cube[i][j][k].Bz = sin(x2 - t2);
//			}
//
//	// СЧИТАЕМ ОШИБКУ
//
//	double abs_err = 0.0, rel_err = 0.0, rel_err_new = 0.0, smape_err = 0.0;
//	double temp, temp1 = 0.0;
//	int i_max = 0, j_max = 0, k_max = 0;
//	int i_max_rel = 0, j_max_rel = 0, k_max_rel = 0;
//	double critical_val;
//
//	if (sizeof(cube[0][0][0].Ey) == 2)
//		critical_val = 0.1;
//	else critical_val = 0.0001;
//
//	for (int i = 0; i < Nx; i++)
//		for (int j = 0; j < Ny; j++)
//			for (int k = 0; k < Nz; k++)
//			{
//				temp1 = abs((double)analitical_cube[i][j][k].Bz - (double)cube[i][j][k].Bz);
//				if (temp1 > abs_err)
//				{
//					abs_err = temp1;
//					i_max = i; j_max = j; k_max = k;
//				}
//				temp1 = abs((double)analitical_cube[i][j][k].Ey - (double)cube[i][j][k].Ey);
//				if (temp1 > abs_err)
//				{
//					abs_err = temp1;
//					i_max = i; j_max = j; k_max = k;
//				}
//
//				if ((double)analitical_cube[i][j][k].Bz > critical_val)
//					rel_err = max(rel_err, abs((double)analitical_cube[i][j][k].Ey - (double)cube[i][j][k].Ey) / (double)analitical_cube[i][j][k].Ey);
//				if ((double)analitical_cube[i][j][k].Bz > critical_val)
//					rel_err = max(rel_err, abs((double)analitical_cube[i][j][k].Bz - (double)cube[i][j][k].Bz) / (double)analitical_cube[i][j][k].Bz);
//
//				temp = (double)2 * abs((double)analitical_cube[i][j][k].Bz - (double)cube[i][j][k].Bz) /
//					(abs((double)analitical_cube[i][j][k].Bz) + abs((double)cube[i][j][k].Bz));
//				smape_err += temp;
//				temp = (double)2 * abs((double)analitical_cube[i][j][k].Ey - (double)cube[i][j][k].Ey) /
//					(abs((double)analitical_cube[i][j][k].Ey) + abs((double)cube[i][j][k].Ey));
//				smape_err += temp;
//
//			}
//
//	// НОВАЯ ВЕРСИЯ ПОДСЧЕТА ОТНОСИТЕЛЬНОЙ ОШИБКИ
//	rel_err_new = max(rel_err_new, abs((double)analitical_cube[i_max][j_max][k_max].Ey - (double)cube[i_max][j_max][k_max].Ey / (double)analitical_cube[i_max][j_max][k_max].Ey));
//	rel_err_new = max(rel_err_new, abs((double)analitical_cube[i_max][j_max][k_max].Bz - (double)cube[i_max][j_max][k_max].Bz / (double)analitical_cube[i_max][j_max][k_max].Bz));
//
//
//	cout << endl << "i_max = " << i_max << "  j_max = " << j_max << " k_max = " << k_max << endl;
//	cout << "Absolute error " << abs_err << endl;
//	cout << "Relative error " << rel_err * (double)100 << endl;
//	cout << "Relative error new version " << rel_err_new * (double)100 << endl;
//	cout << "Smape error " << smape_err / (double)(6 * Nx * Ny * Nz) << endl;
//
//	cout << endl << "Ey = " << cube[0][0][0].Ey << "  Bz = " << cube[0][0][0].Bz << endl;
//	return abs_err;
//
//
//}

template <class ftype>
double FDTD_three_dimen_print_to_file_task1(int Nx, int Ny, int Nz, ftype T, ftype dt, string type_sum, string name_file)
{
	setlocale(LC_ALL, "Russian");

	int Nt = T / dt;

	vector<vector<vector<Component<ftype>>>> cube(Nx + 2, vector<vector<Component<ftype>>>(Ny + 2, vector<Component<ftype>>(Nz + 2)));

	vector<double> vec_abs_err(Nt + 1), vec_rel_err(Nt), vec_smape_err(Nt);
	vector<double> vec_energy(Nt), vec_value(Nt + 1);
	pair<ftype, ftype> ab(0.0, 2.0 * M_PI), cd(0.0, 2.0 * M_PI), fg(0.0, 2.0 * M_PI);
	ftype dx = (ab.second - ab.first) / (ftype)Nx, dy = (cd.second - cd.first) / (ftype)Ny,
		dz = (fg.second - fg.first) / (ftype)Nz;

	ftype dt_x = dt / (dx), dt_y = dt / (dy), dt_z = dt / (dz);
	double abs_err_t = 0.0, rel_err_t = 0.0, smape_err_t = 0.0;

	Initializing_FDTD_three_dimen_task1<ftype>(cube, Nx, Ny, Nz, dx, dt, ab, abs_err_t, rel_err_t, smape_err_t);
	vec_abs_err[0] = abs_err_t;
	vec_smape_err[0] = abs_err_t;
	vec_value[0] = cube[0][0][0].Ey;
	if (type_sum == "Kahan")
	{
	}
	else
	{
		for (int it = 0; it < Nt; it++)
		{
			if (it % 1000 == 0)
				cout << it << endl;
			abs_err_t = 0.0; rel_err_t = 0.0; smape_err_t = 0.0;

#pragma omp parallel for collapse(3)
			for (int i = 0; i < Nx + 1; i++)
				for (int j = 0; j < Ny + 1; j++)
					for (int k = 0; k < Nz + 1; k++)
					{
						Update_electric_field_three_dimen<ftype>(cube, dt_x, dt_y, dt_z, i, j, k);
						Check_error_current_step_electric_task1(cube, dt, dx, ab.first, it, i, j, k, abs_err_t, rel_err_t, smape_err_t);
					}

#pragma omp parallel for collapse(3)
			for (int i = 1; i < Nx + 1; i++)
				for (int j = 1; j < Ny + 1; j++)
					for (int k = 1; k < Nz + 1; k++)
					{
						Update_magnetic_field_three_dimen<ftype>(cube, dt_x, dt_y, dt_z, i, j, k);
						Check_error_current_step_magnetic_task1(cube, dt, dx, ab.first, it, i, j, k, abs_err_t, rel_err_t, smape_err_t);
					}

			// Добавляем значение в вектор
			vec_abs_err[it + 1] = (abs_err_t);
			vec_rel_err[it] = (rel_err_t * (double)100);
			vec_smape_err[it] = (smape_err_t / (double)(6 * Nx * Ny * Nz));
			vec_energy[it] = CalculateEnergy(cube, Nx, Ny, Nz);
			vec_value[it + 1] = cube[0][0][0].Ey;

			// Изменяем значения на границе
			for (int i = 0; i < Ny + 2; i++)
				for (int j = 0; j < Nz + 2; j++)
				{
					cube[Nx + 1][i][j].Ex = cube[1][i][j].Ex;
					cube[Nx + 1][i][j].Ey = cube[1][i][j].Ey;
					cube[Nx + 1][i][j].Ez = cube[1][i][j].Ez;
					cube[Nx + 1][i][j].Bx = cube[1][i][j].Bx;
					cube[Nx + 1][i][j].By = cube[1][i][j].By;
					cube[Nx + 1][i][j].Bz = cube[1][i][j].Bz;

					cube[0][i][j].Ex = cube[Nx][i][j].Ex;
					cube[0][i][j].Ey = cube[Nx][i][j].Ey;
					cube[0][i][j].Ez = cube[Nx][i][j].Ez;
					cube[0][i][j].Bx = cube[Nx][i][j].Bx;
					cube[0][i][j].By = cube[Nx][i][j].By;
					cube[0][i][j].Bz = cube[Nx][i][j].Bz;
				}
			for (int i = 0; i < Nx + 2; i++)
				for (int j = 0; j < Nz + 2; j++)
				{
					cube[i][Ny + 1][j].Ex = cube[i][1][j].Ex;
					cube[i][Ny + 1][j].Ey = cube[i][1][j].Ey;
					cube[i][Ny + 1][j].Ez = cube[i][1][j].Ez;
					cube[i][Ny + 1][j].Bx = cube[i][1][j].Bx;
					cube[i][Ny + 1][j].By = cube[i][1][j].By;
					cube[i][Ny + 1][j].Bz = cube[i][1][j].Bz;

					cube[i][0][j].Ex = cube[i][Ny][j].Ex;
					cube[i][0][j].Ey = cube[i][Ny][j].Ey;
					cube[i][0][j].Ez = cube[i][Ny][j].Ez;
					cube[i][0][j].Bx = cube[i][Ny][j].Bx;
					cube[i][0][j].By = cube[i][Ny][j].By;
					cube[i][0][j].Bz = cube[i][Ny][j].Bz;
				}
			for (int i = 0; i < Nx + 2; i++)
				for (int j = 0; j < Ny + 2; j++)
				{

					cube[i][j][Nz + 1].Ex = cube[i][j][1].Ex;
					cube[i][j][Nz + 1].Ey = cube[i][j][1].Ey;
					cube[i][j][Nz + 1].Ez = cube[i][j][1].Ez;
					cube[i][j][Nz + 1].Bx = cube[i][j][1].Bx;
					cube[i][j][Nz + 1].By = cube[i][j][1].By;
					cube[i][j][Nz + 1].Bz = cube[i][j][1].Bz;

					cube[i][j][0].Ex = cube[i][j][Nz].Ex;
					cube[i][j][0].Ey = cube[i][j][Nz].Ey;
					cube[i][j][0].Ez = cube[i][j][Nz].Ez;
					cube[i][j][0].Bx = cube[i][j][Nz].Bx;
					cube[i][j][0].By = cube[i][j][Nz].By;
					cube[i][j][0].Bz = cube[i][j][Nz].Bz;
				}
		}
	}
	//вывод векторов ошибок в файл
	ofstream numb1(name_file);
	if (sizeof(cube[0][0][0].Ey) == 2)		numb1 << "Тип данных" << "half" << endl;
	if (sizeof(cube[0][0][0].Ey) == 4)		numb1 << "Тип данных" << "float" << endl;
	if (sizeof(cube[0][0][0].Ey) == 8)		numb1 << "Тип данных" << "double" << endl;
	numb1 << "FDTD" << endl;
	numb1 << "Тип суммирования " << type_sum << endl;
	numb1 << "T =" << T << endl;
	numb1 << "Nx = " << Nx << "Ny = " << Ny << "Nz = " << Nz << endl;
	numb1 << "dt = " << dt << endl;
	numb1 << "Nt = " << (int)(T / dt) << endl << endl;
	for (int s = 0; s < Nt + 1; s++)
	{
		numb1 << dt * (double)(s) << ";" << vec_abs_err[s] << ";" << ";" << dt * (double)(s) << ";" << vec_smape_err[s] << endl;
	}
	numb1 << endl << endl;
	numb1.close();

	////вывод вектора значений в файл
	//ofstream numb_value("Value task1 Ey 000.csv");
	//for (int s = 0; s < Nt; s++)
	//{
	//	numb_value << dt * (double)(s + 1) << ";" << vec_value[s] << endl;
	//}
	//numb_value << endl << endl;
	//numb_value.close();

	//вывод вектора энергии в файл
	ofstream numb_energy("Energy_graph_without_PML.csv");
	for (int s = 0; s < Nt; s++)
	{
		numb_energy << dt * (double)(s + 1) << ";" << vec_energy[s]  << endl;
	}
	numb_energy << endl << endl;
	numb_energy.close();

	////График решений в разных плоскостях
	//ofstream numb2("graph_solution.csv");
	//numb2 << ";" << "z" << endl;
	//numb2 << "y" << ";" << ";";
	//for (int i = 0; i < Ny + 1; i++)
	//{
	//	numb2 << dy * (ftype)i << ";";
	//}
	//numb2 << endl;
	//for (int i = 0; i < Nz + 1; i++)
	//{
	//	ftype z = dz * (ftype)i;
	//	numb2 << ";" << z << ";";
	//	for (int j = 0; j < Ny + 1; j++)
	//	{
	//		numb2 << cube[16][j][i].Bz << ";";
	//	}
	//	numb2 << endl;
	//}
	//numb2 << dx * 16.0;
	//numb2.close();

	////График решений в разных плоскостях одномерный
	//ofstream numb3("graph_solution_Bx_.csv");
	//for(int i = 0; i < Ny + 1; i++)
	//{ 
	//	numb3 << dy * (ftype)i << ";"<<cube[0][i][0].Bx << endl;
	//}
	//numb3.close();

	vector<vector<vector<Component<ftype>>>> analitical_cube(Nx + 1, vector<vector<Component<ftype>>>(Ny + 1, vector<Component<ftype>>(Nz + 1)));

	double t1 = (double)Nt * dt;
	double t2 = (double)Nt * dt + 0.5 * dt;
	double x1, x2;

	for (int i = 0; i < Nx + 1; i++)
		for (int j = 0; j < Ny + 1; j++)
			for (int k = 0; k < Nz + 1; k++)
			{
				x1 = ab.first + dx * (double)i + 0.5 * dx;
				x2 = ab.first + dx * (double)i;

				analitical_cube[i][j][k].Ey = sin(x1 - t1);
				analitical_cube[i][j][k].Bz = sin(x2 - t2);
			}

	// СЧИТАЕМ ОШИБКУ

	double abs_err = 0.0, rel_err = 0.0, rel_err_new = 0.0, smape_err = 0.0;
	double temp, temp1 = 0.0;
	int i_max = 0, j_max = 0, k_max = 0;
	int i_max_rel = 0, j_max_rel = 0, k_max_rel = 0;
	double critical_val;

	if (sizeof(cube[0][0][0].Ey) == 2)
		critical_val = 0.1;
	else critical_val = 0.0001;

	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			for (int k = 0; k < Nz; k++)
			{
				temp1 = abs((double)analitical_cube[i][j][k].Bz - (double)cube[i][j][k].Bz);
				if (temp1 > abs_err)
				{
					abs_err = temp1;
					i_max = i; j_max = j; k_max = k;
				}
				temp1 = abs((double)analitical_cube[i][j][k].Ey - (double)cube[i][j][k].Ey);
				if (temp1 > abs_err)
				{
					abs_err = temp1;
					i_max = i; j_max = j; k_max = k;
				}

				if ((double)analitical_cube[i][j][k].Bz > critical_val)
					rel_err = max(rel_err, abs((double)analitical_cube[i][j][k].Ey - (double)cube[i][j][k].Ey) / (double)analitical_cube[i][j][k].Ey);
				if ((double)analitical_cube[i][j][k].Bz > critical_val)
					rel_err = max(rel_err, abs((double)analitical_cube[i][j][k].Bz - (double)cube[i][j][k].Bz) / (double)analitical_cube[i][j][k].Bz);

				temp = (double)2 * abs((double)analitical_cube[i][j][k].Bz - (double)cube[i][j][k].Bz) /
					(abs((double)analitical_cube[i][j][k].Bz) + abs((double)cube[i][j][k].Bz));
				smape_err += temp;
				temp = (double)2 * abs((double)analitical_cube[i][j][k].Ey - (double)cube[i][j][k].Ey) /
					(abs((double)analitical_cube[i][j][k].Ey) + abs((double)cube[i][j][k].Ey));
				smape_err += temp;

			}

	// НОВАЯ ВЕРСИЯ ПОДСЧЕТА ОТНОСИТЕЛЬНОЙ ОШИБКИ
	rel_err_new = max(rel_err_new, abs((double)analitical_cube[i_max][j_max][k_max].Ey - (double)cube[i_max][j_max][k_max].Ey / (double)analitical_cube[i_max][j_max][k_max].Ey));
	rel_err_new = max(rel_err_new, abs((double)analitical_cube[i_max][j_max][k_max].Bz - (double)cube[i_max][j_max][k_max].Bz / (double)analitical_cube[i_max][j_max][k_max].Bz));


	cout << endl << "i_max = " << i_max << "  j_max = " << j_max << " k_max = " << k_max << endl;
	cout << "Absolute error " << abs_err << endl;
	cout << "Relative error " << rel_err * (double)100 << endl;
	cout << "Relative error new version " << rel_err_new * (double)100 << endl;
	cout << "Smape error " << smape_err / (double)(6 * Nx * Ny * Nz) << endl;
	return 0.0;
}

template <class ftype>
double FDTD_three_dimen_print_to_file_task2(int Nx, int Ny, int Nz, ftype T, ftype dt, string type_sum, string name_file)
{
	setlocale(LC_ALL, "Russian");

	int Nt = T / dt;

	vector<vector<vector<Component<ftype>>>> cube(Nx + 2, vector<vector<Component<ftype>>>(Ny + 2, vector<Component<ftype>>(Nz + 2)));

	vector<double> vec_abs_err(Nt + 1), vec_rel_err(Nt), vec_smape_err(Nt);
	vector<double> vec_energy(Nt), vec_value(Nt + 1);
	pair<ftype, ftype> ab(0.0, 2.0 * M_PI), cd(0.0, 2.0 * M_PI), fg(0.0, 2.0 * M_PI);
	ftype dx = (ab.second - ab.first) / (ftype)Nx, dy = (cd.second - cd.first) / (ftype)Ny,
		dz = (fg.second - fg.first) / (ftype)Nz;

	ftype dt_x = dt / (dx), dt_y = dt / (dy), dt_z = dt / (dz);
	double abs_err_t = 0.0, rel_err_t = 0.0, smape_err_t = 0.0;

	Initializing_FDTD_three_dimen_task2<ftype>(cube, Nx, Ny, Nz, dy, dt, cd, abs_err_t, rel_err_t, smape_err_t);
	vec_abs_err[0] = abs_err_t;
	vec_value[0] = cube[0][0][0].Ez;
	if (type_sum == "Kahan")
	{
	}
	else
	{
		for (int it = 0; it < Nt; it++)
		{
			if (it % 100 == 0)
				cout << it << endl;
			abs_err_t = 0.0; rel_err_t = 0.0; smape_err_t = 0.0;

#pragma omp parallel for collapse(3)
			for (int i = 0; i < Nx + 1; i++)
				for (int j = 0; j < Ny + 1; j++)
					for (int k = 0; k < Nz + 1; k++)
					{
						Update_electric_field_three_dimen<ftype>(cube, dt_x, dt_y, dt_z, i, j, k);
						Check_error_current_step_electric_task2(cube, dt, dy, cd.first, it, i, j, k, abs_err_t, rel_err_t, smape_err_t);
					}

#pragma omp parallel for collapse(3)
			for (int i = 1; i < Nx + 1; i++)
				for (int j = 1; j < Ny + 1; j++)
					for (int k = 1; k < Nz + 1; k++)
					{
						Update_magnetic_field_three_dimen<ftype>(cube, dt_x, dt_y, dt_z, i, j, k);
						Check_error_current_step_magnetic_task2(cube, dt, dy, cd.first, it, i, j, k, abs_err_t, rel_err_t, smape_err_t);
					}
			// Добавляем значение в вектор
			vec_abs_err[it + 1] = (abs_err_t);
			vec_rel_err[it] = (rel_err_t * (double)100);
			vec_smape_err[it] = (smape_err_t / (double)(6 * Nx * Ny * Nz));
			vec_energy[it] = CalculateEnergy(cube, Nx, Ny, Nz);
			vec_value[it + 1] = cube[0][0][0].Ez;

			// Изменяем значения на границе
			for (int i = 0; i < Ny + 2; i++)
				for (int j = 0; j < Nz + 2; j++)
				{
					cube[Nx + 1][i][j].Ex = cube[1][i][j].Ex;
					cube[Nx + 1][i][j].Ey = cube[1][i][j].Ey;
					cube[Nx + 1][i][j].Ez = cube[1][i][j].Ez;
					cube[Nx + 1][i][j].Bx = cube[1][i][j].Bx;
					cube[Nx + 1][i][j].By = cube[1][i][j].By;
					cube[Nx + 1][i][j].Bz = cube[1][i][j].Bz;

					cube[0][i][j].Ex = cube[Nx][i][j].Ex;
					cube[0][i][j].Ey = cube[Nx][i][j].Ey;
					cube[0][i][j].Ez = cube[Nx][i][j].Ez;
					cube[0][i][j].Bx = cube[Nx][i][j].Bx;
					cube[0][i][j].By = cube[Nx][i][j].By;
					cube[0][i][j].Bz = cube[Nx][i][j].Bz;
				}
			for (int i = 0; i < Nx + 2; i++)
				for (int j = 0; j < Nz + 2; j++)
				{
					cube[i][Ny + 1][j].Ex = cube[i][1][j].Ex;
					cube[i][Ny + 1][j].Ey = cube[i][1][j].Ey;
					cube[i][Ny + 1][j].Ez = cube[i][1][j].Ez;
					cube[i][Ny + 1][j].Bx = cube[i][1][j].Bx;
					cube[i][Ny + 1][j].By = cube[i][1][j].By;
					cube[i][Ny + 1][j].Bz = cube[i][1][j].Bz;

					cube[i][0][j].Ex = cube[i][Ny][j].Ex;
					cube[i][0][j].Ey = cube[i][Ny][j].Ey;
					cube[i][0][j].Ez = cube[i][Ny][j].Ez;
					cube[i][0][j].Bx = cube[i][Ny][j].Bx;
					cube[i][0][j].By = cube[i][Ny][j].By;
					cube[i][0][j].Bz = cube[i][Ny][j].Bz;
				}
			for (int i = 0; i < Nx + 2; i++)
				for (int j = 0; j < Ny + 2; j++)
				{

					cube[i][j][Nz + 1].Ex = cube[i][j][1].Ex;
					cube[i][j][Nz + 1].Ey = cube[i][j][1].Ey;
					cube[i][j][Nz + 1].Ez = cube[i][j][1].Ez;
					cube[i][j][Nz + 1].Bx = cube[i][j][1].Bx;
					cube[i][j][Nz + 1].By = cube[i][j][1].By;
					cube[i][j][Nz + 1].Bz = cube[i][j][1].Bz;

					cube[i][j][0].Ex = cube[i][j][Nz].Ex;
					cube[i][j][0].Ey = cube[i][j][Nz].Ey;
					cube[i][j][0].Ez = cube[i][j][Nz].Ez;
					cube[i][j][0].Bx = cube[i][j][Nz].Bx;
					cube[i][j][0].By = cube[i][j][Nz].By;
					cube[i][j][0].Bz = cube[i][j][Nz].Bz;
				}
		}
	}
	//вывод векторов ошибок в файл
	ofstream numb1(name_file);
	if (sizeof(cube[0][0][0].Ey) == 2)		numb1 << "Тип данных" << "half" << endl;
	if (sizeof(cube[0][0][0].Ey) == 4)		numb1 << "Тип данных" << "float" << endl;
	if (sizeof(cube[0][0][0].Ey) == 8)		numb1 << "Тип данных" << "double" << endl;
	numb1 << "FDTD" << endl;
	numb1 << "Тип суммирования " << type_sum << endl;
	numb1 << "T =" << T << endl;
	numb1 << "Nx = " << Nx << "Ny = " << Ny << "Nz = " << Nz << endl;
	numb1 << "dt = " << dt << endl;
	numb1 << "Nt = " << (int)(T / dt) << endl << endl;
	for (int s = 0; s < Nt + 1; s++)
	{
		numb1 << dt * (double)(s) << ";" << vec_abs_err[s] << endl;
	}
	numb1 << endl << endl;
	numb1.close();

	////вывод вектора значений в файл
	//ofstream numb_value("Value task2 Ez 000.csv");
	//for (int s = 0; s < Nt + 1; s++)
	//{
	//	numb_value << dt * (double)(s) << ";" << vec_value[s] << endl;
	//}
	//numb_value << endl << endl;
	//numb_value.close();

	//вывод вектора энергии в файл
	ofstream numb_energy("Energy_graph_2.csv");
	for (int s = 0; s < Nt; s++)
	{
		numb_energy << dt * (double)(s + 1) << ";" << vec_energy[s]  << endl;
	}
	numb_energy << endl << endl;
	numb_energy.close();

	////График решений в разных плоскостях
	//ofstream numb2("graph_solution.csv");
	//numb2 << ";" << "z" << endl;
	//numb2 << "y" << ";" << ";";
	//for (int i = 0; i < Ny + 1; i++)
	//{
	//	numb2 << dy * (ftype)i << ";";
	//}
	//numb2 << endl;
	//for (int i = 0; i < Nz + 1; i++)
	//{
	//	ftype z = dz * (ftype)i;
	//	numb2 << ";" << z << ";";
	//	for (int j = 0; j < Ny + 1; j++)
	//	{
	//		numb2 << cube[16][j][i].Bz << ";";
	//	}
	//	numb2 << endl;
	//}
	//numb2 << dx * 16.0;
	//numb2.close();

	////График решений в разных плоскостях одномерный
	//ofstream numb3("graph_solution_Bx_.csv");
	//for(int i = 0; i < Ny + 1; i++)
	//{ 
	//	numb3 << dy * (ftype)i << ";"<<cube[0][i][0].Bx << endl;
	//}
	//numb3.close();

	vector<vector<vector<Component<ftype>>>> analitical_cube(Nx + 1, vector<vector<Component<ftype>>>(Ny + 1, vector<Component<ftype>>(Nz + 1)));

	double t1 = (double)Nt * dt;
	double t2 = (double)Nt * dt + 0.5 * dt;
	double y1, y2;

	for (int i = 0; i < Nx + 1; i++)
		for (int j = 0; j < Ny + 1; j++)
			for (int k = 0; k < Nz + 1; k++)
			{
				y1 = cd.first + dy * (double)j + 0.5 * dy;
				y2 = cd.first + dy * (double)j;

				analitical_cube[i][j][k].Ez = sin(y1 - t1);
				analitical_cube[i][j][k].Bx = sin(y2 - t2);
			}

	// СЧИТАЕМ ОШИБКУ

	double abs_err = 0.0, rel_err = 0.0, rel_err_new = 0.0, smape_err = 0.0;
	double temp, temp1 = 0.0;
	int i_max = 0, j_max = 0, k_max = 0;
	int i_max_rel = 0, j_max_rel = 0, k_max_rel = 0;
	double critical_val;

	if (sizeof(cube[0][0][0].Ey) == 2)
		critical_val = 0.1;
	else critical_val = 0.0001;

	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			for (int k = 0; k < Nz; k++)
			{
				temp1 = abs((double)analitical_cube[i][j][k].Bz - (double)cube[i][j][k].Bz);
				if (temp1 > abs_err)
				{
					abs_err = temp1;
					i_max = i; j_max = j; k_max = k;
				}
				temp1 = abs((double)analitical_cube[i][j][k].Ey - (double)cube[i][j][k].Ey);
				if (temp1 > abs_err)
				{
					abs_err = temp1;
					i_max = i; j_max = j; k_max = k;
				}

				if ((double)analitical_cube[i][j][k].Bz > critical_val)
					rel_err = max(rel_err, abs((double)analitical_cube[i][j][k].Ey - (double)cube[i][j][k].Ey) / (double)analitical_cube[i][j][k].Ey);
				if ((double)analitical_cube[i][j][k].Bz > critical_val)
					rel_err = max(rel_err, abs((double)analitical_cube[i][j][k].Bz - (double)cube[i][j][k].Bz) / (double)analitical_cube[i][j][k].Bz);

				temp = (double)2 * abs((double)analitical_cube[i][j][k].Bz - (double)cube[i][j][k].Bz) /
					(abs((double)analitical_cube[i][j][k].Bz) + abs((double)cube[i][j][k].Bz));
				smape_err += temp;
				temp = (double)2 * abs((double)analitical_cube[i][j][k].Ey - (double)cube[i][j][k].Ey) /
					(abs((double)analitical_cube[i][j][k].Ey) + abs((double)cube[i][j][k].Ey));
				smape_err += temp;

			}

	// НОВАЯ ВЕРСИЯ ПОДСЧЕТА ОТНОСИТЕЛЬНОЙ ОШИБКИ
	rel_err_new = max(rel_err_new, abs((double)analitical_cube[i_max][j_max][k_max].Ey - (double)cube[i_max][j_max][k_max].Ey / (double)analitical_cube[i_max][j_max][k_max].Ey));
	rel_err_new = max(rel_err_new, abs((double)analitical_cube[i_max][j_max][k_max].Bz - (double)cube[i_max][j_max][k_max].Bz / (double)analitical_cube[i_max][j_max][k_max].Bz));


	cout << endl << "i_max = " << i_max << "  j_max = " << j_max << " k_max = " << k_max << endl;
	cout << "Absolute error " << abs_err << endl;
	cout << "Relative error " << rel_err * (double)100 << endl;
	cout << "Relative error new version " << rel_err_new * (double)100 << endl;
	cout << "Smape error " << smape_err / (double)(6 * Nx * Ny * Nz) << endl;
	return 0.0;
}

template <class ftype>
double FDTD_three_dimen_print_to_file_task3(int Nx, int Ny, int Nz, ftype T, ftype dt, string type_sum, string name_file)
{
	setlocale(LC_ALL, "Russian");

	int Nt = T / dt;

	vector<vector<vector<Component<ftype>>>> cube(Nx + 2, vector<vector<Component<ftype>>>(Ny + 2, vector<Component<ftype>>(Nz + 2)));

	vector<double> vec_abs_err(Nt + 1), vec_rel_err(Nt), vec_smape_err(Nt);
	vector<double> vec_energy(Nt), vec_value(Nt);
	pair<ftype, ftype> ab(0.0, 2.0 * M_PI), cd(0.0, 2.0 * M_PI), fg(0.0, 2.0 * M_PI);
	ftype dx = (ab.second - ab.first) / (ftype)Nx, dy = (cd.second - cd.first) / (ftype)Ny,
		dz = (fg.second - fg.first) / (ftype)Nz;

	ftype dt_x = dt / (dx), dt_y = dt / (dy), dt_z = dt / (dz);
	double abs_err_t = 0.0, rel_err_t = 0.0, smape_err_t = 0.0;

	Initializing_FDTD_three_dimen_task3<ftype>(cube, Nx, Ny, Nz, dz, dt, fg, abs_err_t, rel_err_t, smape_err_t);
	vec_abs_err[0] = abs_err_t;

	if (type_sum == "Kahan")
	{
	}
	else
	{
		for (int it = 0; it < Nt; it++)
		{
			if (it % 100 == 0)
				cout << it << endl;
			abs_err_t = 0.0; rel_err_t = 0.0; smape_err_t = 0.0;

#pragma omp parallel for collapse(3)
			for (int i = 0; i < Nx + 1; i++)
				for (int j = 0; j < Ny + 1; j++)
					for (int k = 0; k < Nz + 1; k++)
					{
						Update_electric_field_three_dimen<ftype>(cube, dt_x, dt_y, dt_z, i, j, k);
						Check_error_current_step_electric_task3(cube, dt, dz, fg.first, it, i, j, k, abs_err_t, rel_err_t, smape_err_t);
					}

#pragma omp parallel for collapse(3)
			for (int i = 1; i < Nx + 1; i++)
				for (int j = 1; j < Ny + 1; j++)
					for (int k = 1; k < Nz + 1; k++)
					{
						Update_magnetic_field_three_dimen<ftype>(cube, dt_x, dt_y, dt_z, i, j, k);
						Check_error_current_step_magnetic_task3(cube, dt, dz, fg.first, it, i, j, k, abs_err_t, rel_err_t, smape_err_t);
					}
			// Добавляем значение в вектор
			vec_abs_err[it + 1] = (abs_err_t);
			vec_rel_err[it] = (rel_err_t * (double)100);
			vec_smape_err[it] = (smape_err_t / (double)(6 * Nx * Ny * Nz));
			vec_energy[it] = CalculateEnergy(cube, Nx, Ny, Nz);
			vec_value[it] = cube[0][0][0].Bx;

			// Изменяем значения на границе
			for (int i = 0; i < Ny + 2; i++)
				for (int j = 0; j < Nz + 2; j++)
				{
					cube[Nx + 1][i][j].Ex = cube[1][i][j].Ex;
					cube[Nx + 1][i][j].Ey = cube[1][i][j].Ey;
					cube[Nx + 1][i][j].Ez = cube[1][i][j].Ez;
					cube[Nx + 1][i][j].Bx = cube[1][i][j].Bx;
					cube[Nx + 1][i][j].By = cube[1][i][j].By;
					cube[Nx + 1][i][j].Bz = cube[1][i][j].Bz;

					cube[0][i][j].Ex = cube[Nx][i][j].Ex;
					cube[0][i][j].Ey = cube[Nx][i][j].Ey;
					cube[0][i][j].Ez = cube[Nx][i][j].Ez;
					cube[0][i][j].Bx = cube[Nx][i][j].Bx;
					cube[0][i][j].By = cube[Nx][i][j].By;
					cube[0][i][j].Bz = cube[Nx][i][j].Bz;
				}
			for (int i = 0; i < Nx + 2; i++)
				for (int j = 0; j < Nz + 2; j++)
				{
					cube[i][Ny + 1][j].Ex = cube[i][1][j].Ex;
					cube[i][Ny + 1][j].Ey = cube[i][1][j].Ey;
					cube[i][Ny + 1][j].Ez = cube[i][1][j].Ez;
					cube[i][Ny + 1][j].Bx = cube[i][1][j].Bx;
					cube[i][Ny + 1][j].By = cube[i][1][j].By;
					cube[i][Ny + 1][j].Bz = cube[i][1][j].Bz;

					cube[i][0][j].Ex = cube[i][Ny][j].Ex;
					cube[i][0][j].Ey = cube[i][Ny][j].Ey;
					cube[i][0][j].Ez = cube[i][Ny][j].Ez;
					cube[i][0][j].Bx = cube[i][Ny][j].Bx;
					cube[i][0][j].By = cube[i][Ny][j].By;
					cube[i][0][j].Bz = cube[i][Ny][j].Bz;
				}
			for (int i = 0; i < Nx + 2; i++)
				for (int j = 0; j < Ny + 2; j++)
				{

					cube[i][j][Nz + 1].Ex = cube[i][j][1].Ex;
					cube[i][j][Nz + 1].Ey = cube[i][j][1].Ey;
					cube[i][j][Nz + 1].Ez = cube[i][j][1].Ez;
					cube[i][j][Nz + 1].Bx = cube[i][j][1].Bx;
					cube[i][j][Nz + 1].By = cube[i][j][1].By;
					cube[i][j][Nz + 1].Bz = cube[i][j][1].Bz;

					cube[i][j][0].Ex = cube[i][j][Nz].Ex;
					cube[i][j][0].Ey = cube[i][j][Nz].Ey;
					cube[i][j][0].Ez = cube[i][j][Nz].Ez;
					cube[i][j][0].Bx = cube[i][j][Nz].Bx;
					cube[i][j][0].By = cube[i][j][Nz].By;
					cube[i][j][0].Bz = cube[i][j][Nz].Bz;
				}
		}
	}
	//вывод векторов ошибок в файл
	ofstream numb1(name_file);
	if (sizeof(cube[0][0][0].Ey) == 2)		numb1 << "Тип данных" << "half" << endl;
	if (sizeof(cube[0][0][0].Ey) == 4)		numb1 << "Тип данных" << "float" << endl;
	if (sizeof(cube[0][0][0].Ey) == 8)		numb1 << "Тип данных" << "double" << endl;
	numb1 << "FDTD" << endl;
	numb1 << "Тип суммирования " << type_sum << endl;
	numb1 << "T =" << T << endl;
	numb1 << "Nx = " << Nx << "Ny = " << Ny << "Nz = " << Nz << endl;
	numb1 << "dt = " << dt << endl;
	numb1 << "Nt = " << (int)(T / dt) << endl << endl;
	for (int s = 0; s < Nt + 1; s++)
	{
		numb1 << dt * (double)(s) << ";" << vec_abs_err[s] << endl;
	}
	numb1 << endl << endl;
	numb1.close();

	////вывод вектора значений в файл
	//ofstream numb_value("Value.csv");
	//for (int s = 0; s < Nt; s++)
	//{
	//	numb_value << dt * (double)(s + 1) << ";" << vec_value[s] << endl;
	//}
	//numb_value << endl << endl;
	//numb_value.close();

	//вывод вектора энергии в файл
	ofstream numb_energy("Energy_graph_3.csv");
	for (int s = 0; s < Nt; s++)
	{
		numb_energy << dt * (double)(s + 1) << ";" << vec_energy[s]  << endl;
	}
	numb_energy << endl << endl;
	numb_energy.close();


	vector<vector<vector<Component<ftype>>>> analitical_cube(Nx + 1, vector<vector<Component<ftype>>>(Ny + 1, vector<Component<ftype>>(Nz + 1)));

	double t1 = (double)Nt * dt;
	double t2 = (double)Nt * dt + 0.5 * dt;
	double z1, z2;

	for (int i = 0; i < Nx + 1; i++)
		for (int j = 0; j < Ny + 1; j++)
			for (int k = 0; k < Nz + 1; k++)
			{
				z1 = fg.first + dz * (double)k + 0.5 * dz;
				z2 = fg.first + dz * (double)k;

				analitical_cube[i][j][k].Ex = sin(z1 - t1);
				analitical_cube[i][j][k].By = sin(z2 - t2);
			}

	// СЧИТАЕМ ОШИБКУ

	double abs_err = 0.0, rel_err = 0.0, rel_err_new = 0.0, smape_err = 0.0;
	double temp, temp1 = 0.0;
	int i_max = 0, j_max = 0, k_max = 0;
	int i_max_rel = 0, j_max_rel = 0, k_max_rel = 0;
	double critical_val;

	if (sizeof(cube[0][0][0].Ey) == 2)
		critical_val = 0.1;
	else critical_val = 0.0001;

	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			for (int k = 0; k < Nz; k++)
			{
				temp1 = abs((double)analitical_cube[i][j][k].Bz - (double)cube[i][j][k].Bz);
				if (temp1 > abs_err)
				{
					abs_err = temp1;
					i_max = i; j_max = j; k_max = k;
				}
				temp1 = abs((double)analitical_cube[i][j][k].Ey - (double)cube[i][j][k].Ey);
				if (temp1 > abs_err)
				{
					abs_err = temp1;
					i_max = i; j_max = j; k_max = k;
				}

				if ((double)analitical_cube[i][j][k].Bz > critical_val)
					rel_err = max(rel_err, abs((double)analitical_cube[i][j][k].Ey - (double)cube[i][j][k].Ey) / (double)analitical_cube[i][j][k].Ey);
				if ((double)analitical_cube[i][j][k].Bz > critical_val)
					rel_err = max(rel_err, abs((double)analitical_cube[i][j][k].Bz - (double)cube[i][j][k].Bz) / (double)analitical_cube[i][j][k].Bz);

				temp = (double)2 * abs((double)analitical_cube[i][j][k].Bz - (double)cube[i][j][k].Bz) /
					(abs((double)analitical_cube[i][j][k].Bz) + abs((double)cube[i][j][k].Bz));
				smape_err += temp;
				temp = (double)2 * abs((double)analitical_cube[i][j][k].Ey - (double)cube[i][j][k].Ey) /
					(abs((double)analitical_cube[i][j][k].Ey) + abs((double)cube[i][j][k].Ey));
				smape_err += temp;

			}

	// НОВАЯ ВЕРСИЯ ПОДСЧЕТА ОТНОСИТЕЛЬНОЙ ОШИБКИ
	rel_err_new = max(rel_err_new, abs((double)analitical_cube[i_max][j_max][k_max].Ey - (double)cube[i_max][j_max][k_max].Ey / (double)analitical_cube[i_max][j_max][k_max].Ey));
	rel_err_new = max(rel_err_new, abs((double)analitical_cube[i_max][j_max][k_max].Bz - (double)cube[i_max][j_max][k_max].Bz / (double)analitical_cube[i_max][j_max][k_max].Bz));


	cout << endl << "i_max = " << i_max << "  j_max = " << j_max << " k_max = " << k_max << endl;
	cout << "Absolute error " << abs_err << endl;
	cout << "Relative error " << rel_err * (double)100 << endl;
	cout << "Relative error new version " << rel_err_new * (double)100 << endl;
	cout << "Smape error " << smape_err / (double)(6 * Nx * Ny * Nz) << endl;
	return 0.0;
}

template <class ftype>
double CalculateEnergy(vector<vector<vector<Component<ftype>>>>& cube, int Nx, int Ny, int Nz)
{
	double energy = 0.0;

	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			for (int k = 0; k < Nz; k++)
			{
				energy += cube[i][j][k].Ex * cube[i][j][k].Ex + cube[i][j][k].Ey * cube[i][j][k].Ey + cube[i][j][k].Ez * cube[i][j][k].Ez;
				energy += cube[i][j][k].Bx * cube[i][j][k].Bx + cube[i][j][k].By * cube[i][j][k].By + cube[i][j][k].Bz * cube[i][j][k].Bz;
			}
	return energy;
}

template <class ftype>
void Check_Curant(ftype dx, ftype dy, ftype dz, ftype dt)
{
	double temp = dt * sqrt(1.0/(dx*dx) + 1.0/(dy*dy) + 1.0/(dz*dz)) ;
	cout << temp << " < 1  ";
	if (temp < 1)
		cout << "TRUE" << endl;
	else cout << "FALSE" << endl;
}