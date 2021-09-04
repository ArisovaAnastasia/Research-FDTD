#define _USE_MATH_DEFINES
#include <algorithm>
#include <iostream>
#include <fstream>
#include <clocale>
#include <chrono>
#include <vector>
#include <cmath>
#include <omp.h>
#include "type data.h"

using namespace std;


template <class ftype>
void Initializing_the_cube_central(vector<vector<vector<Component<ftype>>>>& cube, long long int _N, ftype _dx,
	pair <ftype, ftype> ab, pair <ftype, ftype> cd, pair <ftype, ftype> fg)
{
	for (long long int i = 0; i < _N; i++)
	{
		for (long long int j = 0; j < _N; j++)
		{
			for (long long int k = 0; k < _N; k++)
			{
				//Аналитическое решение
				//Ex = 0	Ey = sin(x)		Ez=0
				//Bx = 0	By = 0			Bz = sin(x)

				ftype x = (ftype)(ab.first + (ftype)i * _dx);

				cube[i][j][k].Ex = (ftype)0.0;
				cube[i][j][k].Ey = sin(x);
				cube[i][j][k].Ez = (ftype)0.0;
				cube[i][j][k].Bx = (ftype)0.0;
				cube[i][j][k].By = (ftype)0.0;
				cube[i][j][k].Bz = sin(x);
			}
		}
	}
	//Заполняем отдельно границу куба
	for (long long int i = 0; i < _N + 2; i++) {
		for (long long int j = 0; j < _N + 2; j++)
		{
			cube[_N][i][j].Ex = cube[0][i][j].Ex;
			cube[_N][i][j].Ey = cube[0][i][j].Ey;
			cube[_N][i][j].Ez = cube[0][i][j].Ez;
			cube[_N][i][j].Bx = cube[0][i][j].Bx;
			cube[_N][i][j].By = cube[0][i][j].By;
			cube[_N][i][j].Bz = cube[0][i][j].Bz;

			cube[_N + 1][i][j].Ex = cube[1][i][j].Ex;
			cube[_N + 1][i][j].Ey = cube[1][i][j].Ey;
			cube[_N + 1][i][j].Ez = cube[1][i][j].Ez;
			cube[_N + 1][i][j].Bx = cube[1][i][j].Bx;
			cube[_N + 1][i][j].By = cube[1][i][j].By;
			cube[_N + 1][i][j].Bz = cube[1][i][j].Bz;

		}
	}
	for (long long int i = 0; i < _N + 2; i++) {
		for (long long int j = 0; j < _N + 2; j++)
		{
			cube[i][_N][j].Ex = cube[i][0][j].Ex;
			cube[i][_N][j].Ey = cube[i][0][j].Ey;
			cube[i][_N][j].Ez = cube[i][0][j].Ez;
			cube[i][_N][j].Bx = cube[i][0][j].Bx;
			cube[i][_N][j].By = cube[i][0][j].By;
			cube[i][_N][j].Bz = cube[i][0][j].Bz;

			cube[i][_N + 1][j].Ex = cube[i][1][j].Ex;
			cube[i][_N + 1][j].Ey = cube[i][1][j].Ey;
			cube[i][_N + 1][j].Ez = cube[i][1][j].Ez;
			cube[i][_N + 1][j].Bx = cube[i][1][j].Bx;
			cube[i][_N + 1][j].By = cube[i][1][j].By;
			cube[i][_N + 1][j].Bz = cube[i][1][j].Bz;

		}
	}

	for (long long int i = 0; i < _N + 2; i++) {
		for (long long int j = 0; j < _N + 2; j++)
		{
			cube[i][j][_N].Ex = cube[i][j][0].Ex;
			cube[i][j][_N].Ey = cube[i][j][0].Ey;
			cube[i][j][_N].Ez = cube[i][j][0].Ez;
			cube[i][j][_N].Bx = cube[i][j][0].Bx;
			cube[i][j][_N].By = cube[i][j][0].By;
			cube[i][j][_N].Bz = cube[i][j][0].Bz;

			cube[i][j][_N + 1].Ex = cube[i][j][1].Ex;
			cube[i][j][_N + 1].Ey = cube[i][j][1].Ey;
			cube[i][j][_N + 1].Ez = cube[i][j][1].Ez;
			cube[i][j][_N + 1].Bx = cube[i][j][1].Bx;
			cube[i][j][_N + 1].By = cube[i][j][1].By;
			cube[i][j][_N + 1].Bz = cube[i][j][1].Bz;

		}
	}
}

template <class ftype>
void Update_cell_central(vector<vector<vector<Component<ftype>>>>& cube, vector<vector<vector<Component<ftype>>>>& cube2,
	ftype _dt_x, ftype _dt_y, ftype _dt_z, ftype dt, ftype dx, ftype ab_first, long long int it,
	long long int i, long long int j, long long int k)
{
	ftype tEx, tEy, tEz, tBx, tBy, tBz;
	ftype  J = (ftype)0;

	tBx = (cube[i][j - 1][k].Ez - cube[i][j + 1][k].Ez) * (_dt_y) / (ftype)2.0
		+ (cube[i][j][k + 1].Ey - cube[i][j][k - 1].Ey) * (_dt_z) / (ftype)2.0 + cube[i][j][k].Bx;
	tBy = (cube[i][j][k - 1].Ex - cube[i][j][k + 1].Ex) * (_dt_z) / (ftype)2.0
		+ (cube[i + 1][j][k].Ez - cube[i - 1][j][k].Ez) * (_dt_x) / (ftype)2.0 + cube[i][j][k].By;
	tBz = (cube[i - 1][j][k].Ey - cube[i + 1][j][k].Ey) * (_dt_x) / (ftype)2.0
		+ (cube[i][j + 1][k].Ex - cube[i][j - 1][k].Ex) * (_dt_y) / (ftype)2.0 + cube[i][j][k].Bz;

	cube2[i][j][k].Bx = tBx;
	cube2[i][j][k].By = tBy;
	cube2[i][j][k].Bz = tBz;

	tEx = (cube[i][j + 1][k].Bz - cube[i][j - 1][k].Bz) * (_dt_y) / (ftype)2.0
		- (cube[i][j][k + 1].By - cube[i][j][k - 1].By) * (_dt_z) / (ftype)2.0 + cube[i][j][k].Ex - J * dt;
	tEy = (cube[i][j][k + 1].Bx - cube[i][j][k - 1].Bx) * (_dt_z) / (ftype)2.0
		-(cube[i + 1][j][k].Bz - cube[i - 1][j][k].Bz) * (_dt_x) / (ftype)2.0 + cube[i][j][k].Ey - J * dt;
	tEz = (cube[i + 1][j][k].By - cube[i - 1][j][k].By) * (_dt_x) / (ftype)2.0
		- (cube[i][j + 1][k].Bx - cube[i][j - 1][k].Bx) * (_dt_y) / (ftype)2.0 + cube[i][j][k].Ez - J * dt;

	cube2[i][j][k].Ex = tEx;
	cube2[i][j][k].Ey = tEy;
	cube2[i][j][k].Ez = tEz;
}

template <class ftype>
void Update_cell_central_with_KahanSum(vector<vector<vector<Component<ftype>>>>& cube, vector<vector<vector<Component<ftype>>>>& cube2,
	vector<vector<vector<Component<ftype>>>>& compensator, vector<vector<vector<Component<ftype>>>>& compensator2, ftype _dt_x, ftype _dt_y, ftype _dt_z, ftype dt, ftype dx, ftype ab_first,
	long long int it, long long int i, long long int j, long long int k)
{
	ftype tEx, tEy, tEz, tBx, tBy, tBz;
	ftype  J = (ftype)0;
	ftype sum, y;

	sum = cube[i][j][k].Bx;
	y = (cube[i][j - 1][k].Ez - cube[i][j + 1][k].Ez - (compensator[i][j - 1][k].Ez - compensator[i][j + 1][k].Ez)) * (_dt_y) / (ftype)2.0
		+ (cube[i][j][k + 1].Ey - cube[i][j][k - 1].Ey - (compensator[i][j][k + 1].Ey - compensator[i][j][k - 1].Ey)) * (_dt_z) / (ftype)2.0
		-compensator[i][j][k].Bx;
	tBx = sum + y;
	compensator2[i][j][k].Bx = (tBx - sum) - y;

	sum = cube[i][j][k].By;
	y = (cube[i][j][k - 1].Ex - cube[i][j][k + 1].Ex - (compensator[i][j][k - 1].Ex - compensator[i][j][k + 1].Ex)) * (_dt_z) / (ftype)2.0
		+(cube[i + 1][j][k].Ez - cube[i - 1][j][k].Ez - (compensator[i + 1][j][k].Ez - compensator[i - 1][j][k].Ez)) * (_dt_x) / (ftype)2.0
		-compensator[i][j][k].By;
	tBy = sum + y;
	compensator2[i][j][k].By = (tBy - sum) - y;

	sum = cube[i][j][k].Bz;
	y = (cube[i - 1][j][k].Ey - cube[i + 1][j][k].Ey - (compensator[i - 1][j][k].Ey - compensator[i + 1][j][k].Ey)) * (_dt_x) / (ftype)2.0
		+(cube[i][j + 1][k].Ex - cube[i][j - 1][k].Ex - (compensator[i][j + 1][k].Ex - compensator[i][j - 1][k].Ex)) * (_dt_y) / (ftype)2.0
		-compensator[i][j][k].Bz;
	tBz = sum + y;
	compensator2[i][j][k].Bz = (tBz - sum) - y;

	cube2[i][j][k].Bx = tBx;
	cube2[i][j][k].By = tBy;
	cube2[i][j][k].Bz = tBz;

	sum = cube[i][j][k].Ex;
	y = (cube[i][j + 1][k].Bz - cube[i][j - 1][k].Bz - (compensator[i][j + 1][k].Bz - compensator[i][j - 1][k].Bz)) * (_dt_y) / (ftype)2.0
		-(cube[i][j][k + 1].By - cube[i][j][k - 1].By - (compensator[i][j][k + 1].By - compensator[i][j][k - 1].By)) * (_dt_z) / (ftype)2.0
		-J * (dt)-compensator[i][j][k].Ex;
	tEx = sum + y;
	compensator2[i][j][k].Ex = (tEx - sum) - y;

	sum = cube[i][j][k].Ey;
	y = (cube[i][j][k + 1].Bx - cube[i][j][k - 1].Bx - (compensator[i][j][k + 1].Bx - compensator[i][j][k - 1].Bx)) * (_dt_z) / (ftype)2.0
		-(cube[i + 1][j][k].Bz - cube[i - 1][j][k].Bz - (compensator[i + 1][j][k].Bz - compensator[i - 1][j][k].Bz)) * (_dt_x) / (ftype)2.0
		-J * dt - compensator[i][j][k].Ey;
	tEy = sum + y;
	compensator2[i][j][k].Ey = (tEy - sum) - y;

	sum = cube[i][j][k].Ez;
	y = (cube[i + 1][j][k].By - cube[i - 1][j][k].By - (compensator[i + 1][j][k].By - compensator[i - 1][j][k].By)) * (_dt_x) / (ftype)2.0
		-(cube[i][j + 1][k].Bx - cube[i][j - 1][k].Bx - (compensator[i][j + 1][k].Bx - compensator[i][j - 1][k].Bx)) * (_dt_y) / (ftype)2.0
		-J * (dt)-compensator[i][j][k].Ez;
	tEz = sum + y;
	compensator2[i][j][k].Ez = (tEz - sum) - y;

	cube2[i][j][k].Ex = tEx;
	cube2[i][j][k].Ey = tEy;
	cube2[i][j][k].Ez = tEz;
}

template <class ftype>
void Update_cell_central_with_error_calculation(vector<vector<vector<Component<ftype>>>>& cube, vector<vector<vector<Component<ftype>>>>& cube2,
	ftype _dt_x, ftype _dt_y, ftype _dt_z, ftype dt, ftype dx, ftype ab_first, long long int it,
	long long int i, long long int j, long long int k, double& abs_err_t, double& rel_err_t, double& smape_err_t)
{
	ftype tEx, tEy, tEz, tBx, tBy, tBz;
	ftype  J = (ftype)0;

	tBx = (cube[i][j - 1][k].Ez - cube[i][j + 1][k].Ez) * (_dt_y) / (ftype)2.0
		+ (cube[i][j][k + 1].Ey - cube[i][j][k - 1].Ey) * (_dt_z) / (ftype)2.0 + cube[i][j][k].Bx;
	tBy = (cube[i][j][k - 1].Ex - cube[i][j][k + 1].Ex) * (_dt_z) / (ftype)2.0
		+ (cube[i + 1][j][k].Ez - cube[i - 1][j][k].Ez) * (_dt_x) / (ftype)2.0 + cube[i][j][k].By;
	tBz = (cube[i - 1][j][k].Ey - cube[i + 1][j][k].Ey) * (_dt_x) / (ftype)2.0
		+ (cube[i][j + 1][k].Ex - cube[i][j - 1][k].Ex) * (_dt_y) / (ftype)2.0 + cube[i][j][k].Bz;

	cube2[i][j][k].Bx = tBx;
	cube2[i][j][k].By = tBy;
	cube2[i][j][k].Bz = tBz;

	tEx = (cube[i][j + 1][k].Bz - cube[i][j - 1][k].Bz) * (_dt_y) / (ftype)2.0
		- (cube[i][j][k + 1].By - cube[i][j][k - 1].By) * (_dt_z) / (ftype)2.0 + cube[i][j][k].Ex - J * dt;
	tEy = (cube[i][j][k + 1].Bx - cube[i][j][k - 1].Bx) * (_dt_z) / (ftype)2.0
		-(cube[i + 1][j][k].Bz - cube[i - 1][j][k].Bz) * (_dt_x) / (ftype)2.0 + cube[i][j][k].Ey - J * dt;
	tEz = (cube[i + 1][j][k].By - cube[i - 1][j][k].By) * (_dt_x) / (ftype)2.0
		- (cube[i][j + 1][k].Bx - cube[i][j - 1][k].Bx) * (_dt_y) / (ftype)2.0 + cube[i][j][k].Ez - J * dt;

	cube2[i][j][k].Ex = tEx;
	cube2[i][j][k].Ey = tEy;
	cube2[i][j][k].Ez = tEz;

	//Считаем текущую абсолютную ошибку в данном узле
	double x = (double)ab_first + (double)i * (double)dx, t = (double)it * (double)dt;

	abs_err_t = max(abs_err_t, abs(((double)sin(x - t)) - (double)cube2[i][j][k].Bz));
	abs_err_t = max(abs_err_t, abs(((double)sin(x - t)) - (double)cube2[i][j][k].Ey));

	//Считаем текущую относительную ошибку в данном узле
	if (t != 0)
	{
		rel_err_t = max(rel_err_t, abs(sin(x - t) - (double)cube2[i][j][k].Bz) / abs(sin(x - t)));
		rel_err_t = max(rel_err_t, abs(sin(x - t) - (double)cube2[i][j][k].Ey) / abs(sin(x - t)));
	}

	//Считаем текущую ошибку smape в данном узле
	double temp;
	temp = (double)2 * abs((double)sin(x - t) - (double)cube2[i][j][k].Bz) /
		(abs((double)sin(x - t)) + abs((double)cube2[i][j][k].Bz));
	smape_err_t += temp;
	temp = (double)2 * abs((double)sin(x - t) - (double)cube2[i][j][k].Ey) /
		(abs((double)sin(x - t)) + abs((double)cube2[i][j][k].Ey));
	smape_err_t += temp;
}

template <class ftype>
void Update_cell_central_with_KahanSum_and_error_calculation(vector<vector<vector<Component<ftype>>>>& cube, vector<vector<vector<Component<ftype>>>>& cube2,
	vector<vector<vector<Component<ftype>>>>& compensator, vector<vector<vector<Component<ftype>>>>& compensator2,  ftype _dt_x, ftype _dt_y, ftype _dt_z, ftype dt, ftype dx, ftype ab_first,
	long long int it, long long int i, long long int j, long long int k, double& abs_err_t, double& rel_err_t, double& smape_err_t)
{
	ftype tEx, tEy, tEz, tBx, tBy, tBz;
	ftype  J = (ftype)0;
	ftype sum, y;

	sum = cube[i][j][k].Bx;
	y = (cube[i][j - 1][k].Ez - cube[i][j + 1][k].Ez - (compensator[i][j - 1][k].Ez - compensator[i][j + 1][k].Ez)) * (_dt_y) / (ftype)2.0
		+ (cube[i][j][k + 1].Ey - cube[i][j][k - 1].Ey - (compensator[i][j][k + 1].Ey - compensator[i][j][k - 1].Ey)) * (_dt_z) / (ftype)2.0
		- compensator[i][j][k].Bx;
	tBx = sum + y;
	compensator2[i][j][k].Bx = (tBx - sum) - y;

	sum = cube[i][j][k].By;
	y = (cube[i][j][k - 1].Ex - cube[i][j][k + 1].Ex - (compensator[i][j][k - 1].Ex - compensator[i][j][k + 1].Ex)) * (_dt_z) / (ftype)2.0
		+ (cube[i + 1][j][k].Ez - cube[i - 1][j][k].Ez - (compensator[i + 1][j][k].Ez - compensator[i - 1][j][k].Ez)) * (_dt_x) / (ftype)2.0
		- compensator[i][j][k].By;
	tBy = sum + y;
	compensator2[i][j][k].By = (tBy - sum) - y;

	sum = cube[i][j][k].Bz;
	y = (cube[i - 1][j][k].Ey - cube[i + 1][j][k].Ey - (compensator[i - 1][j][k].Ey - compensator[i + 1][j][k].Ey)) * (_dt_x) / (ftype)2.0
		+ (cube[i][j + 1][k].Ex - cube[i][j - 1][k].Ex - (compensator[i][j + 1][k].Ex - compensator[i][j - 1][k].Ex)) * (_dt_y) / (ftype)2.0
		- compensator[i][j][k].Bz;
	tBz = sum + y;
	compensator2[i][j][k].Bz = (tBz - sum) - y;

	cube2[i][j][k].Bx = tBx;
	cube2[i][j][k].By = tBy;
	cube2[i][j][k].Bz = tBz;

	sum = cube[i][j][k].Ex;
	y = (cube[i][j + 1][k].Bz - cube[i][j - 1][k].Bz - (compensator[i][j + 1][k].Bz - compensator[i][j - 1][k].Bz)) * (_dt_y) / (ftype)2.0
		- (cube[i][j][k + 1].By - cube[i][j][k - 1].By - (compensator[i][j][k + 1].By - compensator[i][j][k - 1].By)) * (_dt_z) / (ftype)2.0
		- J * (dt)-compensator[i][j][k].Ex;
	tEx = sum + y;
	compensator2[i][j][k].Ex = (tEx - sum) - y;

	sum = cube[i][j][k].Ey;
	y = (cube[i][j][k + 1].Bx - cube[i][j][k - 1].Bx - (compensator[i][j][k + 1].Bx - compensator[i][j][k - 1].Bx)) * (_dt_z) / (ftype)2.0
		- (cube[i + 1][j][k].Bz - cube[i - 1][j][k].Bz - (compensator[i + 1][j][k].Bz - compensator[i - 1][j][k].Bz)) * (_dt_x) / (ftype)2.0
		- J * dt - compensator[i][j][k].Ey;
	tEy = sum + y;
	compensator2[i][j][k].Ey = (tEy - sum) - y;

	sum = cube[i][j][k].Ez;
	y = (cube[i + 1][j][k].By - cube[i - 1][j][k].By - (compensator[i + 1][j][k].By - compensator[i - 1][j][k].By)) * (_dt_x) / (ftype)2.0
		- (cube[i][j + 1][k].Bx - cube[i][j - 1][k].Bx - (compensator[i][j + 1][k].Bx - compensator[i][j - 1][k].Bx)) * (_dt_y) / (ftype)2.0
		- J * (dt)-compensator[i][j][k].Ez;
	tEz = sum + y;
	compensator2[i][j][k].Ez = (tEz - sum) - y;

	cube2[i][j][k].Ex = tEx;
	cube2[i][j][k].Ey = tEy;
	cube2[i][j][k].Ez = tEz;

	//Считаем текущую абсолютную ошибку в данном узле
	double x = (double)ab_first + (double)i * (double)dx, t = (double)it * (double)dt;

	abs_err_t = max(abs_err_t, abs(((double)sin(x - t)) - (double)cube2[i][j][k].Bz));
	abs_err_t = max(abs_err_t, abs(((double)sin(x - t)) - (double)cube2[i][j][k].Ey));

	//Считаем текущую относительную ошибку в данном узле
	if (t != 0)
	{
		rel_err_t = max(rel_err_t, abs(sin(x - t) - (double)cube2[i][j][k].Bz) / abs(sin(x - t)));
		rel_err_t = max(rel_err_t, abs(sin(x - t) - (double)cube2[i][j][k].Ey) / abs(sin(x - t)));
	}

	//Считаем текущую ошибку smape в данном узле
	double temp;
	temp = (double)2 * abs((double)sin(x - t) - (double)cube2[i][j][k].Bz) /
		(abs((double)sin(x - t)) + abs((double)cube2[i][j][k].Bz));
	smape_err_t += temp;
	temp = (double)2 * abs((double)sin(x - t) - (double)cube2[i][j][k].Ey) /
		(abs((double)sin(x - t)) + abs((double)cube2[i][j][k].Ey));
	smape_err_t += temp;
}


template <class ftype>
double Maxwells_central_solution_for_measuring_time(long long int _N, ftype T, ftype _dt, string type_sum)
{

	setlocale(LC_ALL, "Russian");

	long long int N_3 = _N;
	long long int Nt = T / _dt;

	vector<vector<vector<Component<ftype>>>> cube(N_3 + 2, vector<vector<Component<ftype>>>(N_3 + 2, vector<Component<ftype>>(N_3 + 2)));
	vector<vector<vector<Component<ftype>>>> cube2(N_3 + 2, vector<vector<Component<ftype>>>(N_3 + 2, vector<Component<ftype>>(N_3 + 2)));
	vector<vector<vector<Component<ftype>>>> compensator(N_3 + 2, vector<vector<Component<ftype>>>(N_3 + 2, vector<Component<ftype>>(N_3 + 2)));
	vector<vector<vector<Component<ftype>>>> compensator2(N_3 + 2, vector<vector<Component<ftype>>>(N_3 + 2, vector<Component<ftype>>(N_3 + 2)));

	pair<ftype, ftype> ab(0.0, 2.0 * M_PI), cd(0.0, 2.0 * M_PI), fg(0.0, 2.0 * M_PI);
	ftype dt = _dt, dx = (ab.second - ab.first) / (ftype)N_3, dy = (cd.second - cd.first) / (ftype)N_3,
		dz = (fg.second - fg.first) / (ftype)N_3;

	ftype dt_x = dt / (dx);
	ftype dt_y = dt / (dy);
	ftype dt_z = dt / (dz);

	Initializing_the_cube_central<ftype>(cube, N_3, dx, ab, cd, fg);
	if (type_sum == "Kahan")
	{
		auto begin = chrono::high_resolution_clock::now();

		for (long long int it = 0; it < Nt; it++)
		{

//#pragma omp parallel for collapse(3)
			for (long long int i = 1; i < N_3 + 1; i++)
				for (long long int j = 1; j < N_3 + 1; j++)
					for (long long int k = 1; k < N_3 + 1; k++)
					{
						Update_cell_central_with_KahanSum<ftype>(cube, cube2, compensator, compensator2, dt_x, dt_y, dt_z, dt, dx, ab.first, it, i, j, k);
					}
//#pragma omp parallel for collapse(3)
			for (long long int i = 1; i < N_3 + 1; i++)
				for (long long int j = 1; j < N_3 + 1; j++)
					for (long long int k = 1; k < N_3 + 1; k++)
					{
						cube[i][j][k].Ex = cube2[i][j][k].Ex;
						cube[i][j][k].Ey = cube2[i][j][k].Ey;
						cube[i][j][k].Ez = cube2[i][j][k].Ez;
						cube[i][j][k].Bx = cube2[i][j][k].Bx;
						cube[i][j][k].By = cube2[i][j][k].By;
						cube[i][j][k].Bz = cube2[i][j][k].Bz;
					}
//#pragma omp parallel for collapse(3)
			for (long long int i = 0; i < N_3 + 2; i++)
				for (long long int j = 0; j < N_3 + 2; j++)
					for (long long int k = 0; k < N_3 + 2; k++)
					{
						compensator[i][j][k].Ex = compensator2[i][j][k].Ex;
						compensator[i][j][k].Ey = compensator2[i][j][k].Ey;
						compensator[i][j][k].Ez = compensator2[i][j][k].Ez;
						compensator[i][j][k].Bx = compensator2[i][j][k].Bx;
						compensator[i][j][k].By = compensator2[i][j][k].By;
						compensator[i][j][k].Bz = compensator2[i][j][k].Bz;
					}
			// Изменяем значения на границе
//#pragma omp parallel for collapse(2)
			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{
					cube[N_3 + 1][i][j].Ex = cube[1][i][j].Ex;
					cube[N_3 + 1][i][j].Ey = cube[1][i][j].Ey;
					cube[N_3 + 1][i][j].Ez = cube[1][i][j].Ez;
					cube[N_3 + 1][i][j].Bx = cube[1][i][j].Bx;
					cube[N_3 + 1][i][j].By = cube[1][i][j].By;
					cube[N_3 + 1][i][j].Bz = cube[1][i][j].Bz;

					cube[0][i][j].Ex = cube[N_3][i][j].Ex;
					cube[0][i][j].Ey = cube[N_3][i][j].Ey;
					cube[0][i][j].Ez = cube[N_3][i][j].Ez;
					cube[0][i][j].Bx = cube[N_3][i][j].Bx;
					cube[0][i][j].By = cube[N_3][i][j].By;
					cube[0][i][j].Bz = cube[N_3][i][j].Bz;
				}
			}
//#pragma omp parallel for collapse(2)
			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{
					cube[i][N_3 + 1][j].Ex = cube[i][1][j].Ex;
					cube[i][N_3 + 1][j].Ey = cube[i][1][j].Ey;
					cube[i][N_3 + 1][j].Ez = cube[i][1][j].Ez;
					cube[i][N_3 + 1][j].Bx = cube[i][1][j].Bx;
					cube[i][N_3 + 1][j].By = cube[i][1][j].By;
					cube[i][N_3 + 1][j].Bz = cube[i][1][j].Bz;

					cube[i][0][j].Ex = cube[i][N_3][j].Ex;
					cube[i][0][j].Ey = cube[i][N_3][j].Ey;
					cube[i][0][j].Ez = cube[i][N_3][j].Ez;
					cube[i][0][j].Bx = cube[i][N_3][j].Bx;
					cube[i][0][j].By = cube[i][N_3][j].By;
					cube[i][0][j].Bz = cube[i][N_3][j].Bz;
				}
			}
//#pragma omp parallel for collapse(2)
			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{

					cube[i][j][N_3 + 1].Ex = cube[i][j][1].Ex;
					cube[i][j][N_3 + 1].Ey = cube[i][j][1].Ey;
					cube[i][j][N_3 + 1].Ez = cube[i][j][1].Ez;
					cube[i][j][N_3 + 1].Bx = cube[i][j][1].Bx;
					cube[i][j][N_3 + 1].By = cube[i][j][1].By;
					cube[i][j][N_3 + 1].Bz = cube[i][j][1].Bz;

					cube[i][j][0].Ex = cube[i][j][N_3].Ex;
					cube[i][j][0].Ey = cube[i][j][N_3].Ey;
					cube[i][j][0].Ez = cube[i][j][N_3].Ez;
					cube[i][j][0].Bx = cube[i][j][N_3].Bx;
					cube[i][j][0].By = cube[i][j][N_3].By;
					cube[i][j][0].Bz = cube[i][j][N_3].Bz;
				}
			}
		}

		auto end = chrono::high_resolution_clock::now();
		return chrono::duration_cast<chrono::milliseconds>(end - begin).count();
	}
	else
	{
		auto begin = chrono::high_resolution_clock::now();

		for (long long int it = 0; it < Nt; it++)
		{

//#pragma omp parallel for collapse(3)
			for (long long int i = 1; i < N_3 + 1; i++)
				for (long long int j = 1; j < N_3 + 1; j++)
					for (long long int k = 1; k < N_3 + 1; k++)
					{
						Update_cell_central<ftype>(cube, cube2, dt_x, dt_y, dt_z, dt, dx, ab.first, it, i, j, k);
					}
//#pragma omp parallel for collapse(3)
			for (long long int i = 1; i < N_3 + 1; i++)
				for (long long int j = 1; j < N_3 + 1; j++)
					for (long long int k = 1; k < N_3 + 1; k++)
					{
						cube[i][j][k].Ex = cube2[i][j][k].Ex;
						cube[i][j][k].Ey = cube2[i][j][k].Ey;
						cube[i][j][k].Ez = cube2[i][j][k].Ez;
						cube[i][j][k].Bx = cube2[i][j][k].Bx;
						cube[i][j][k].By = cube2[i][j][k].By;
						cube[i][j][k].Bz = cube2[i][j][k].Bz;
					}

			// Изменяем значения на границе
//#pragma omp parallel for collapse(2)
			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{
					cube[N_3 + 1][i][j].Ex = cube[1][i][j].Ex;
					cube[N_3 + 1][i][j].Ey = cube[1][i][j].Ey;
					cube[N_3 + 1][i][j].Ez = cube[1][i][j].Ez;
					cube[N_3 + 1][i][j].Bx = cube[1][i][j].Bx;
					cube[N_3 + 1][i][j].By = cube[1][i][j].By;
					cube[N_3 + 1][i][j].Bz = cube[1][i][j].Bz;

					cube[0][i][j].Ex = cube[N_3][i][j].Ex;
					cube[0][i][j].Ey = cube[N_3][i][j].Ey;
					cube[0][i][j].Ez = cube[N_3][i][j].Ez;
					cube[0][i][j].Bx = cube[N_3][i][j].Bx;
					cube[0][i][j].By = cube[N_3][i][j].By;
					cube[0][i][j].Bz = cube[N_3][i][j].Bz;
				}
			}
//#pragma omp parallel for collapse(2)
			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{
					cube[i][N_3 + 1][j].Ex = cube[i][1][j].Ex;
					cube[i][N_3 + 1][j].Ey = cube[i][1][j].Ey;
					cube[i][N_3 + 1][j].Ez = cube[i][1][j].Ez;
					cube[i][N_3 + 1][j].Bx = cube[i][1][j].Bx;
					cube[i][N_3 + 1][j].By = cube[i][1][j].By;
					cube[i][N_3 + 1][j].Bz = cube[i][1][j].Bz;

					cube[i][0][j].Ex = cube[i][N_3][j].Ex;
					cube[i][0][j].Ey = cube[i][N_3][j].Ey;
					cube[i][0][j].Ez = cube[i][N_3][j].Ez;
					cube[i][0][j].Bx = cube[i][N_3][j].Bx;
					cube[i][0][j].By = cube[i][N_3][j].By;
					cube[i][0][j].Bz = cube[i][N_3][j].Bz;
				}
			}
//#pragma omp parallel for collapse(2)
			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{

					cube[i][j][N_3 + 1].Ex = cube[i][j][1].Ex;
					cube[i][j][N_3 + 1].Ey = cube[i][j][1].Ey;
					cube[i][j][N_3 + 1].Ez = cube[i][j][1].Ez;
					cube[i][j][N_3 + 1].Bx = cube[i][j][1].Bx;
					cube[i][j][N_3 + 1].By = cube[i][j][1].By;
					cube[i][j][N_3 + 1].Bz = cube[i][j][1].Bz;

					cube[i][j][0].Ex = cube[i][j][N_3].Ex;
					cube[i][j][0].Ey = cube[i][j][N_3].Ey;
					cube[i][j][0].Ez = cube[i][j][N_3].Ez;
					cube[i][j][0].Bx = cube[i][j][N_3].Bx;
					cube[i][j][0].By = cube[i][j][N_3].By;
					cube[i][j][0].Bz = cube[i][j][N_3].Bz;
				}
			}
		}

		auto end = chrono::high_resolution_clock::now();
		return chrono::duration_cast<chrono::milliseconds>(end - begin).count();
	}
}

template <class ftype>
double Maxwells_central_solution_print_error_to_console(long long int _N, ftype T, ftype _dt, string type_sum)
{

	setlocale(LC_ALL, "Russian");

	long long int N_3 = _N;
	long long int Nt = T / _dt;

	vector<vector<vector<Component<ftype>>>> cube(N_3 + 2, vector<vector<Component<ftype>>>(N_3 + 2, vector<Component<ftype>>(N_3 + 2)));
	vector<vector<vector<Component<ftype>>>> cube2(N_3 + 2, vector<vector<Component<ftype>>>(N_3 + 2, vector<Component<ftype>>(N_3 + 2)));
	vector<vector<vector<Component<ftype>>>> compensator(N_3 + 2, vector<vector<Component<ftype>>>(N_3 + 2, vector<Component<ftype>>(N_3 + 2)));
	vector<vector<vector<Component<ftype>>>> compensator2(N_3 + 2, vector<vector<Component<ftype>>>(N_3 + 2, vector<Component<ftype>>(N_3 + 2)));
	
	pair<ftype, ftype> ab(0.0, 2.0 * M_PI), cd(0.0, 2.0 * M_PI), fg(0.0, 2.0 * M_PI);
	ftype dt = _dt, dx = (ab.second - ab.first) / (ftype)N_3, dy = (cd.second - cd.first) / (ftype)N_3,
		dz = (fg.second - fg.first) / (ftype)N_3;

	ftype dt_x = dt / (dx);
	ftype dt_y = dt / (dy);
	ftype dt_z = dt / (dz);

	Initializing_the_cube_central<ftype> (cube, N_3, dx, ab, cd, fg);

	if (type_sum == "Kahan")
	{
		for (long long int it = 0; it < Nt; it++)
		{
			if (it % 100 == 0)
				cout << it << endl;
			#pragma omp parallel for collapse(3)
			for (long long int i = 1; i < N_3 + 1; i++)
				for (long long int j = 1; j < N_3 + 1; j++)
					for (long long int k = 1; k < N_3 + 1; k++)
					{
						Update_cell_central_with_KahanSum<ftype>(cube, cube2, compensator, compensator2, dt_x, dt_y, dt_z, dt, dx, ab.first, it, i, j, k);
					}

			for (long long int i = 1; i < N_3 + 1; i++)
				for (long long int j = 1; j < N_3 + 1; j++)
					for (long long int k = 1; k < N_3 + 1; k++)
					{
						cube[i][j][k].Ex = cube2[i][j][k].Ex;
						cube[i][j][k].Ey = cube2[i][j][k].Ey;
						cube[i][j][k].Ez = cube2[i][j][k].Ez;
						cube[i][j][k].Bx = cube2[i][j][k].Bx;
						cube[i][j][k].By = cube2[i][j][k].By;
						cube[i][j][k].Bz = cube2[i][j][k].Bz;
					}
			for (long long int i = 0; i < N_3 + 2; i++)
				for (long long int j = 0; j < N_3 + 2; j++)
					for (long long int k = 0; k < N_3 + 2; k++)
					{
						compensator[i][j][k].Ex = compensator2[i][j][k].Ex;
						compensator[i][j][k].Ey = compensator2[i][j][k].Ey;
						compensator[i][j][k].Ez = compensator2[i][j][k].Ez;
						compensator[i][j][k].Bx = compensator2[i][j][k].Bx;
						compensator[i][j][k].By = compensator2[i][j][k].By;
						compensator[i][j][k].Bz = compensator2[i][j][k].Bz;
					}
			// Изменяем значения на границе
			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{
					cube[N_3 + 1][i][j].Ex = cube[1][i][j].Ex;
					cube[N_3 + 1][i][j].Ey = cube[1][i][j].Ey;
					cube[N_3 + 1][i][j].Ez = cube[1][i][j].Ez;
					cube[N_3 + 1][i][j].Bx = cube[1][i][j].Bx;
					cube[N_3 + 1][i][j].By = cube[1][i][j].By;
					cube[N_3 + 1][i][j].Bz = cube[1][i][j].Bz;

					cube[0][i][j].Ex = cube[N_3][i][j].Ex;
					cube[0][i][j].Ey = cube[N_3][i][j].Ey;
					cube[0][i][j].Ez = cube[N_3][i][j].Ez;
					cube[0][i][j].Bx = cube[N_3][i][j].Bx;
					cube[0][i][j].By = cube[N_3][i][j].By;
					cube[0][i][j].Bz = cube[N_3][i][j].Bz;
				}
			}

			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{
					cube[i][N_3 + 1][j].Ex = cube[i][1][j].Ex;
					cube[i][N_3 + 1][j].Ey = cube[i][1][j].Ey;
					cube[i][N_3 + 1][j].Ez = cube[i][1][j].Ez;
					cube[i][N_3 + 1][j].Bx = cube[i][1][j].Bx;
					cube[i][N_3 + 1][j].By = cube[i][1][j].By;
					cube[i][N_3 + 1][j].Bz = cube[i][1][j].Bz;

					cube[i][0][j].Ex = cube[i][N_3][j].Ex;
					cube[i][0][j].Ey = cube[i][N_3][j].Ey;
					cube[i][0][j].Ez = cube[i][N_3][j].Ez;
					cube[i][0][j].Bx = cube[i][N_3][j].Bx;
					cube[i][0][j].By = cube[i][N_3][j].By;
					cube[i][0][j].Bz = cube[i][N_3][j].Bz;
				}
			}

			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{

					cube[i][j][N_3 + 1].Ex = cube[i][j][1].Ex;
					cube[i][j][N_3 + 1].Ey = cube[i][j][1].Ey;
					cube[i][j][N_3 + 1].Ez = cube[i][j][1].Ez;
					cube[i][j][N_3 + 1].Bx = cube[i][j][1].Bx;
					cube[i][j][N_3 + 1].By = cube[i][j][1].By;
					cube[i][j][N_3 + 1].Bz = cube[i][j][1].Bz;

					cube[i][j][0].Ex = cube[i][j][N_3].Ex;
					cube[i][j][0].Ey = cube[i][j][N_3].Ey;
					cube[i][j][0].Ez = cube[i][j][N_3].Ez;
					cube[i][j][0].Bx = cube[i][j][N_3].Bx;
					cube[i][j][0].By = cube[i][j][N_3].By;
					cube[i][j][0].Bz = cube[i][j][N_3].Bz;
				}
			}
		}
	}
	else 
	{
		for (long long int it = 0; it < Nt; it++)
		{
			if (it % 100 == 0)
				cout << it << endl;
			#pragma omp parallel for collapse(3)
			for (long long int i = 1; i < N_3 + 1; i++)
				for (long long int j = 1; j < N_3 + 1; j++)
					for (long long int k = 1; k < N_3 + 1; k++)
					{
						Update_cell_central<ftype>(cube, cube2, dt_x, dt_y, dt_z, dt, dx, ab.first, it, i, j, k);
					}

			for (long long int i = 1; i < N_3 + 1; i++)
				for (long long int j = 1; j < N_3 + 1; j++)
					for (long long int k = 1; k < N_3 + 1; k++)
					{
						cube[i][j][k].Ex = cube2[i][j][k].Ex;
						cube[i][j][k].Ey = cube2[i][j][k].Ey;
						cube[i][j][k].Ez = cube2[i][j][k].Ez;
						cube[i][j][k].Bx = cube2[i][j][k].Bx;
						cube[i][j][k].By = cube2[i][j][k].By;
						cube[i][j][k].Bz = cube2[i][j][k].Bz;
					}

			// Изменяем значения на границе
			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{
					cube[N_3 + 1][i][j].Ex = cube[1][i][j].Ex;
					cube[N_3 + 1][i][j].Ey = cube[1][i][j].Ey;
					cube[N_3 + 1][i][j].Ez = cube[1][i][j].Ez;
					cube[N_3 + 1][i][j].Bx = cube[1][i][j].Bx;
					cube[N_3 + 1][i][j].By = cube[1][i][j].By;
					cube[N_3 + 1][i][j].Bz = cube[1][i][j].Bz;

					cube[0][i][j].Ex = cube[N_3][i][j].Ex;
					cube[0][i][j].Ey = cube[N_3][i][j].Ey;
					cube[0][i][j].Ez = cube[N_3][i][j].Ez;
					cube[0][i][j].Bx = cube[N_3][i][j].Bx;
					cube[0][i][j].By = cube[N_3][i][j].By;
					cube[0][i][j].Bz = cube[N_3][i][j].Bz;
				}
			}

			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{
					cube[i][N_3 + 1][j].Ex = cube[i][1][j].Ex;
					cube[i][N_3 + 1][j].Ey = cube[i][1][j].Ey;
					cube[i][N_3 + 1][j].Ez = cube[i][1][j].Ez;
					cube[i][N_3 + 1][j].Bx = cube[i][1][j].Bx;
					cube[i][N_3 + 1][j].By = cube[i][1][j].By;
					cube[i][N_3 + 1][j].Bz = cube[i][1][j].Bz;

					cube[i][0][j].Ex = cube[i][N_3][j].Ex;
					cube[i][0][j].Ey = cube[i][N_3][j].Ey;
					cube[i][0][j].Ez = cube[i][N_3][j].Ez;
					cube[i][0][j].Bx = cube[i][N_3][j].Bx;
					cube[i][0][j].By = cube[i][N_3][j].By;
					cube[i][0][j].Bz = cube[i][N_3][j].Bz;
				}
			}

			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{

					cube[i][j][N_3 + 1].Ex = cube[i][j][1].Ex;
					cube[i][j][N_3 + 1].Ey = cube[i][j][1].Ey;
					cube[i][j][N_3 + 1].Ez = cube[i][j][1].Ez;
					cube[i][j][N_3 + 1].Bx = cube[i][j][1].Bx;
					cube[i][j][N_3 + 1].By = cube[i][j][1].By;
					cube[i][j][N_3 + 1].Bz = cube[i][j][1].Bz;

					cube[i][j][0].Ex = cube[i][j][N_3].Ex;
					cube[i][j][0].Ey = cube[i][j][N_3].Ey;
					cube[i][j][0].Ez = cube[i][j][N_3].Ez;
					cube[i][j][0].Bx = cube[i][j][N_3].Bx;
					cube[i][j][0].By = cube[i][j][N_3].By;
					cube[i][j][0].Bz = cube[i][j][N_3].Bz;
				}
			}
		}
	}

	vector<vector<vector<Component<ftype>>>> analitical_cube(N_3 + 1, vector<vector<Component<ftype>>>(N_3 + 1, vector<Component<ftype>>(N_3 + 1)));

	double t = (double)Nt * dt;
	for (long long int i = 0; i < N_3 + 1; i++)
	{
		for (long long int j = 0; j < N_3 + 1; j++)
		{
			for (long long int k = 0; k < N_3 + 1; k++)
			{
				ftype x = (ftype)(ab.first + (ftype)i * dx);

				analitical_cube[i][j][k].Bz = sin(x - t);
				analitical_cube[i][j][k].Ey = sin(x - t);
			}
		}
	}

	// СЧИТАЕМ ОШИБКУ

	double abs_err = 0.0, rel_err = 0.0, rel_err_new = 0.0;
	double smape_err = 0.0;
	double temp, temp1 = 0.0;
	int i_max = 0, j_max = 0, k_max = 0;
	int i_max_rel = 0, j_max_rel = 0, k_max_rel = 0;
	double critical_val;

	if (sizeof(cube[0][0][0].Ey) == 2)
		critical_val = 0.1;
	else critical_val = 0.0001;

	for (long long int i = 0; i < N_3; i++)
	{
		for (long long int j = 0; j < N_3; j++)
		{
			for (long long int k = 0; k < N_3; k++)
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
		}
	}
	// НОВАЯ ВЕРСИЯ ПОДСЧЕТА ОТНОСИТЕЛЬНОЙ ОШИБКИ
	rel_err_new = max(rel_err_new, abs((double)analitical_cube[i_max][j_max][k_max].Ey - (double)cube[i_max][j_max][k_max].Ey / (double)analitical_cube[i_max][j_max][k_max].Ey));
	rel_err_new = max(rel_err_new, abs((double)analitical_cube[i_max][j_max][k_max].Bz - (double)cube[i_max][j_max][k_max].Bz / (double)analitical_cube[i_max][j_max][k_max].Bz));


	cout << endl << "i_max = " << i_max << "  j_max = " << j_max << " k_max = " << k_max << endl;
	cout << "Absolute error " << abs_err << endl;
	cout << "Relative error " << rel_err * (double)100 << endl;
	cout << "Relative error new version " << rel_err_new * (double)100 << endl;
	cout << "Smape error " << smape_err / (double)(6 * N_3 * N_3 * N_3) << endl;
	return 0.0;
}

template <class ftype>
double Maxwells_central_solution_print_error_to_file(long long int _N, ftype T, ftype _dt, string type_sum, string name_file)
{
	setlocale(LC_ALL, "Russian");

	long long int N_3 = _N;
	long long int Nt = T / _dt;
	vector<double> vec_abs_err(Nt), vec_rel_err(Nt), vec_smape_err(Nt);
	vector<double> comp_Ey(Nt), comp_analit(Nt);

	vector<vector<vector<Component<ftype>>>> cube(N_3 + 2, vector<vector<Component<ftype>>>(N_3 + 2, vector<Component<ftype>>(N_3 + 2)));
	vector<vector<vector<Component<ftype>>>> cube2(N_3 + 2, vector<vector<Component<ftype>>>(N_3 + 2, vector<Component<ftype>>(N_3 + 2)));
	vector<vector<vector<Component<ftype>>>> compensator(N_3 + 2, vector<vector<Component<ftype>>>(N_3 + 2, vector<Component<ftype>>(N_3 + 2)));
	vector<vector<vector<Component<ftype>>>> compensator2(N_3 + 2, vector<vector<Component<ftype>>>(N_3 + 2, vector<Component<ftype>>(N_3 + 2)));

	pair<ftype, ftype> ab(0.0, 2.0 * M_PI), cd(0.0, 2.0 * M_PI), fg(0.0, 2.0 * M_PI);
	ftype dt = _dt, dx = (ab.second - ab.first) / (ftype)N_3, dy = (cd.second - cd.first) / (ftype)N_3,
		dz = (fg.second - fg.first) / (ftype)N_3;

	ftype dt_x = dt / (dx);
	ftype dt_y = dt / (dy);
	ftype dt_z = dt / (dz);

	Initializing_the_cube_central<ftype>(cube, N_3, dx, ab, cd, fg);

	double abs_err_t, rel_err_t, smape_err_t;

	ofstream numb2("graph_solution.csv");
	numb2 << ";" << "t" << endl;
	numb2 << "x" << ";" << ";";
	for (int i = 0; i < N_3 + 1; i++)
	{
		numb2 << dx * (ftype)i << ";";
	}
	numb2 << endl;

	if (type_sum == "Kahan")
	{
		for (long long int it = 0; it < Nt; it++)
		{
			if (it % 100 == 0)
				cout << it << endl;
			abs_err_t = 0.0; rel_err_t = 0.0; smape_err_t = 0.0;

			#pragma omp parallel for collapse(3)
			for (long long int i = 1; i < N_3 + 1; i++)
				for (long long int j = 1; j < N_3 + 1; j++)
					for (long long int k = 1; k < N_3 + 1; k++)
					{
						Update_cell_central_with_KahanSum_and_error_calculation<ftype>(cube, cube2, compensator, compensator2, dt_x, dt_y, dt_z, dt, dx, ab.first, it, i, j, k, abs_err_t, rel_err_t, smape_err_t);
					}

			comp_Ey[it] = cube2[1][1][1].Ey;
			comp_analit[it] = sin(-it * dt);

			for (long long int i = 1; i < N_3 + 1; i++)
				for (long long int j = 1; j < N_3 + 1; j++)
					for (long long int k = 1; k < N_3 + 1; k++)
					{
						cube[i][j][k].Ex = cube2[i][j][k].Ex;
						cube[i][j][k].Ey = cube2[i][j][k].Ey;
						cube[i][j][k].Ez = cube2[i][j][k].Ez;
						cube[i][j][k].Bx = cube2[i][j][k].Bx;
						cube[i][j][k].By = cube2[i][j][k].By;
						cube[i][j][k].Bz = cube2[i][j][k].Bz;
					}

			 //Добавляем значение в вектор

			vec_abs_err[it] = (abs_err_t);

			vec_rel_err[it] = (rel_err_t * (double)100);

			vec_smape_err[it] = (smape_err_t / (double)(6 * N_3 * N_3 * N_3));

			for (long long int i = 0; i < N_3 + 2; i++)
				for (long long int j = 0; j < N_3 + 2; j++)
					for (long long int k = 0; k < N_3 + 2; k++)
					{
						compensator[i][j][k].Ex = compensator2[i][j][k].Ex;
						compensator[i][j][k].Ey = compensator2[i][j][k].Ey;
						compensator[i][j][k].Ez = compensator2[i][j][k].Ez;
						compensator[i][j][k].Bx = compensator2[i][j][k].Bx;
						compensator[i][j][k].By = compensator2[i][j][k].By;
						compensator[i][j][k].Bz = compensator2[i][j][k].Bz;
					}
			// Изменяем значения на границе
			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{
					cube[N_3 + 1][i][j].Ex = cube[1][i][j].Ex;
					cube[N_3 + 1][i][j].Ey = cube[1][i][j].Ey;
					cube[N_3 + 1][i][j].Ez = cube[1][i][j].Ez;
					cube[N_3 + 1][i][j].Bx = cube[1][i][j].Bx;
					cube[N_3 + 1][i][j].By = cube[1][i][j].By;
					cube[N_3 + 1][i][j].Bz = cube[1][i][j].Bz;

					cube[0][i][j].Ex = cube[N_3][i][j].Ex;
					cube[0][i][j].Ey = cube[N_3][i][j].Ey;
					cube[0][i][j].Ez = cube[N_3][i][j].Ez;
					cube[0][i][j].Bx = cube[N_3][i][j].Bx;
					cube[0][i][j].By = cube[N_3][i][j].By;
					cube[0][i][j].Bz = cube[N_3][i][j].Bz;
				}
			}

			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{
					cube[i][N_3 + 1][j].Ex = cube[i][1][j].Ex;
					cube[i][N_3 + 1][j].Ey = cube[i][1][j].Ey;
					cube[i][N_3 + 1][j].Ez = cube[i][1][j].Ez;
					cube[i][N_3 + 1][j].Bx = cube[i][1][j].Bx;
					cube[i][N_3 + 1][j].By = cube[i][1][j].By;
					cube[i][N_3 + 1][j].Bz = cube[i][1][j].Bz;

					cube[i][0][j].Ex = cube[i][N_3][j].Ex;
					cube[i][0][j].Ey = cube[i][N_3][j].Ey;
					cube[i][0][j].Ez = cube[i][N_3][j].Ez;
					cube[i][0][j].Bx = cube[i][N_3][j].Bx;
					cube[i][0][j].By = cube[i][N_3][j].By;
					cube[i][0][j].Bz = cube[i][N_3][j].Bz;
				}
			}

			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{

					cube[i][j][N_3 + 1].Ex = cube[i][j][1].Ex;
					cube[i][j][N_3 + 1].Ey = cube[i][j][1].Ey;
					cube[i][j][N_3 + 1].Ez = cube[i][j][1].Ez;
					cube[i][j][N_3 + 1].Bx = cube[i][j][1].Bx;
					cube[i][j][N_3 + 1].By = cube[i][j][1].By;
					cube[i][j][N_3 + 1].Bz = cube[i][j][1].Bz;

					cube[i][j][0].Ex = cube[i][j][N_3].Ex;
					cube[i][j][0].Ey = cube[i][j][N_3].Ey;
					cube[i][j][0].Ez = cube[i][j][N_3].Ez;
					cube[i][j][0].Bx = cube[i][j][N_3].Bx;
					cube[i][j][0].By = cube[i][j][N_3].By;
					cube[i][j][0].Bz = cube[i][j][N_3].Bz;
				}
			}
		}
	}
	else
	{
		for (long long int it = 0; it < Nt; it++)
		{
			if (it % 100 == 0)
				cout << it << endl;
			abs_err_t = 0.0; rel_err_t = 0.0; smape_err_t = 0.0;

			#pragma omp parallel for collapse(3)
			for (long long int i = 1; i < N_3 + 1; i++)
				for (long long int j = 1; j < N_3 + 1; j++)
					for (long long int k = 1; k < N_3 + 1; k++)
					{
						Update_cell_central_with_error_calculation<ftype>(cube, cube2, dt_x, dt_y, dt_z, dt, dx, ab.first, it, i, j, k, abs_err_t, rel_err_t, smape_err_t);
					}

			comp_Ey[it] = cube2[1][1][1].Ey;
			comp_analit[it] = sin(-it * dt);

			for (long long int i = 1; i < N_3 + 1; i++)
				for (long long int j = 1; j < N_3 + 1; j++)
					for (long long int k = 1; k < N_3 + 1; k++)
					{
						cube[i][j][k].Ex = cube2[i][j][k].Ex;
						cube[i][j][k].Ey = cube2[i][j][k].Ey;
						cube[i][j][k].Ez = cube2[i][j][k].Ez;
						cube[i][j][k].Bx = cube2[i][j][k].Bx;
						cube[i][j][k].By = cube2[i][j][k].By;
						cube[i][j][k].Bz = cube2[i][j][k].Bz;
					}

			// Добавляем значение в вектор

			vec_abs_err[it] = (abs_err_t);

			vec_rel_err[it] = (rel_err_t * (double)100);

			vec_smape_err[it] = (smape_err_t / (double)(6 * N_3 * N_3 * N_3));

			// Изменяем значения на границе
			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{
					cube[N_3 + 1][i][j].Ex = cube[1][i][j].Ex;
					cube[N_3 + 1][i][j].Ey = cube[1][i][j].Ey;
					cube[N_3 + 1][i][j].Ez = cube[1][i][j].Ez;
					cube[N_3 + 1][i][j].Bx = cube[1][i][j].Bx;
					cube[N_3 + 1][i][j].By = cube[1][i][j].By;
					cube[N_3 + 1][i][j].Bz = cube[1][i][j].Bz;

					cube[0][i][j].Ex = cube[N_3][i][j].Ex;
					cube[0][i][j].Ey = cube[N_3][i][j].Ey;
					cube[0][i][j].Ez = cube[N_3][i][j].Ez;
					cube[0][i][j].Bx = cube[N_3][i][j].Bx;
					cube[0][i][j].By = cube[N_3][i][j].By;
					cube[0][i][j].Bz = cube[N_3][i][j].Bz;
				}
			}

			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{
					cube[i][N_3 + 1][j].Ex = cube[i][1][j].Ex;
					cube[i][N_3 + 1][j].Ey = cube[i][1][j].Ey;
					cube[i][N_3 + 1][j].Ez = cube[i][1][j].Ez;
					cube[i][N_3 + 1][j].Bx = cube[i][1][j].Bx;
					cube[i][N_3 + 1][j].By = cube[i][1][j].By;
					cube[i][N_3 + 1][j].Bz = cube[i][1][j].Bz;

					cube[i][0][j].Ex = cube[i][N_3][j].Ex;
					cube[i][0][j].Ey = cube[i][N_3][j].Ey;
					cube[i][0][j].Ez = cube[i][N_3][j].Ez;
					cube[i][0][j].Bx = cube[i][N_3][j].Bx;
					cube[i][0][j].By = cube[i][N_3][j].By;
					cube[i][0][j].Bz = cube[i][N_3][j].Bz;
				}
			}

			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{

					cube[i][j][N_3 + 1].Ex = cube[i][j][1].Ex;
					cube[i][j][N_3 + 1].Ey = cube[i][j][1].Ey;
					cube[i][j][N_3 + 1].Ez = cube[i][j][1].Ez;
					cube[i][j][N_3 + 1].Bx = cube[i][j][1].Bx;
					cube[i][j][N_3 + 1].By = cube[i][j][1].By;
					cube[i][j][N_3 + 1].Bz = cube[i][j][1].Bz;

					cube[i][j][0].Ex = cube[i][j][N_3].Ex;
					cube[i][j][0].Ey = cube[i][j][N_3].Ey;
					cube[i][j][0].Ez = cube[i][j][N_3].Ez;
					cube[i][j][0].Bx = cube[i][j][N_3].Bx;
					cube[i][j][0].By = cube[i][j][N_3].By;
					cube[i][j][0].Bz = cube[i][j][N_3].Bz;
				}
			}

			numb2 << ";" << dt*(double)it << ";";
			for (int j = 0; j < N_3 + 1; j++)
			{
				numb2 << cube[j][0][0].Ey << ";";
			}
			numb2 << endl;

		}
	}
	
	//вывод векторов в файл

	ofstream numb1(name_file);
	
	if (sizeof(cube[0][0][0].Ey) == 2)		numb1 << "Тип данных" << "half" << endl;
	if (sizeof(cube[0][0][0].Ey) == 4)		numb1 << "Тип данных" << "float" << endl;
	if (sizeof(cube[0][0][0].Ey) == 8)		numb1 << "Тип данных" << "double" << endl;
	numb1 << "Central operator" << endl;
	numb1 << "Тип суммирования " << type_sum << endl;
	numb1 << "T =" << T << endl;
	numb1 << "N = " << N_3 << endl;
	numb1 << "dt = " << dt << endl;
	numb1 << "Nt = " << (int)(T / dt) << endl << endl;

	for (int s = 0; s < Nt; s++)
	{
		numb1 << dt * (double)s << ";" << vec_abs_err[s] << ";" << vec_rel_err[s] << ";" << vec_smape_err[s] << endl;
	}
	numb1 << endl << endl;

	numb1.close();

	//ofstream numb1(name_file);

	//if (sizeof(cube[0][0][0].Ey) == 2)		numb1 << "Тип данных" << "half" << endl;
	//if (sizeof(cube[0][0][0].Ey) == 4)		numb1 << "Тип данных" << "float" << endl;
	//if (sizeof(cube[0][0][0].Ey) == 8)		numb1 << "Тип данных" << "double" << endl;
	//numb1 << "Central operator" << endl;
	//numb1 << "Тип суммирования " << type_sum << endl;
	//numb1 << "N = " << N_3 << endl;
	//numb1 << "T = " << T << endl;
	//numb1 << "Nt " << Nt << endl;
	//numb1 << "dt = " << dt << endl << endl;

	//for (int s = 0; s < Nt; s++)
	//{
	//	numb1 << dt*(double)s << ";" << comp_Ey[s]<< ";" << comp_analit[s] << endl;
	//}
	//numb1 << endl << endl;

	//numb1.close();




	numb2.close();


	vector<vector<vector<Component<ftype>>>> analitical_cube(N_3 + 1, vector<vector<Component<ftype>>>(N_3 + 1, vector<Component<ftype>>(N_3 + 1)));

	double t = (double)Nt * dt;
	for (long long int i = 0; i < N_3 + 1; i++)
	{
		for (long long int j = 0; j < N_3 + 1; j++)
		{
			for (long long int k = 0; k < N_3 + 1; k++)
			{
				ftype x = (ftype)(ab.first + (ftype)i * dx);

				analitical_cube[i][j][k].Bz = sin(x - t);
				analitical_cube[i][j][k].Ey = sin(x - t);
			}
		}
	}

	 //СЧИТАЕМ ОШИБКУ

	double abs_err = 0.0, rel_err = 0.0, rel_err_new = 0.0;
	double smape_err = 0.0;
	double temp, temp1 = 0.0;
	int i_max = 0, j_max = 0, k_max = 0;
	int i_max_rel = 0, j_max_rel = 0, k_max_rel = 0;
	double critical_val;

	if (sizeof(cube[0][0][0].Ey) == 2)
		critical_val = 0.1;
	else critical_val = 0.0001;

	for (long long int i = 0; i < N_3; i++)
	{
		for (long long int j = 0; j < N_3; j++)
		{
			for (long long int k = 0; k < N_3; k++)
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
		}
	}
	 //НОВАЯ ВЕРСИЯ ПОДСЧЕТА ОТНОСИТЕЛЬНОЙ ОШИБКИ
	rel_err_new = max(rel_err_new, abs((double)analitical_cube[i_max][j_max][k_max].Ey - (double)cube[i_max][j_max][k_max].Ey / (double)analitical_cube[i_max][j_max][k_max].Ey));
	rel_err_new = max(rel_err_new, abs((double)analitical_cube[i_max][j_max][k_max].Bz - (double)cube[i_max][j_max][k_max].Bz / (double)analitical_cube[i_max][j_max][k_max].Bz));


	cout << endl << "i_max = " << i_max << "  j_max = " << j_max << " k_max = " << k_max << endl;
	cout << "Absolute error " << abs_err << endl;
	cout << "Relative error " << rel_err * (double)100 << endl;
	cout << "Relative error new version " << rel_err_new * (double)100 << endl;
	cout << "Smape error " << smape_err / (double)(6 * N_3 * N_3 * N_3) << endl;
	return 0.0;
}

template <class ftype>
void Min_time_multy_array(long long int _N, ftype T, ftype dt)
{
	ftype time = (ftype)1000000;

	for (long long int i = 0; i < 1; i++)
	{
		ftype temp = (ftype)(Maxwells_central_solution_for_measuring_time<ftype>(_N, T, dt, "Kahan"));
		if (temp < time)
			time = temp;
	}
	cout << time << " ms" << endl;
}

template <class ftype>
double Maxwells_central_solution_for_depended(long long int _N, ftype T, ftype _dt, string type_sum, double& abs_err_return, double& rel_err_return, double& smape_err_return)
{

	setlocale(LC_ALL, "Russian");

	long long int N_3 = _N;
	long long int Nt = T / _dt;

	vector<vector<vector<Component<ftype>>>> cube(N_3 + 2, vector<vector<Component<ftype>>>(N_3 + 2, vector<Component<ftype>>(N_3 + 2)));
	vector<vector<vector<Component<ftype>>>> cube2(N_3 + 2, vector<vector<Component<ftype>>>(N_3 + 2, vector<Component<ftype>>(N_3 + 2)));
	vector<vector<vector<Component<ftype>>>> compensator(N_3 + 2, vector<vector<Component<ftype>>>(N_3 + 2, vector<Component<ftype>>(N_3 + 2)));
	vector<vector<vector<Component<ftype>>>> compensator2(N_3 + 2, vector<vector<Component<ftype>>>(N_3 + 2, vector<Component<ftype>>(N_3 + 2)));

	pair<ftype, ftype> ab(0.0, 2.0 * M_PI), cd(0.0, 2.0 * M_PI), fg(0.0, 2.0 * M_PI);
	ftype dt = _dt, dx = (ab.second - ab.first) / (ftype)N_3, dy = (cd.second - cd.first) / (ftype)N_3,
		dz = (fg.second - fg.first) / (ftype)N_3;

	ftype dt_x = dt / (dx);
	ftype dt_y = dt / (dy);
	ftype dt_z = dt / (dz);

	Initializing_the_cube_central<ftype>(cube, N_3, dx, ab, cd, fg);

	if (type_sum == "Kahan")
	{
		for (long long int it = 0; it < Nt; it++)
		{
			if (it % 100 == 0)
				cout << it << endl;
			#pragma omp parallel for collapse(3)
			for (long long int i = 1; i < N_3 + 1; i++)
				for (long long int j = 1; j < N_3 + 1; j++)
					for (long long int k = 1; k < N_3 + 1; k++)
					{
						Update_cell_central_with_KahanSum<ftype>(cube, cube2, compensator, compensator2, dt_x, dt_y, dt_z, dt, dx, ab.first, it, i, j, k);
					}

			for (long long int i = 1; i < N_3 + 1; i++)
				for (long long int j = 1; j < N_3 + 1; j++)
					for (long long int k = 1; k < N_3 + 1; k++)
					{
						cube[i][j][k].Ex = cube2[i][j][k].Ex;
						cube[i][j][k].Ey = cube2[i][j][k].Ey;
						cube[i][j][k].Ez = cube2[i][j][k].Ez;
						cube[i][j][k].Bx = cube2[i][j][k].Bx;
						cube[i][j][k].By = cube2[i][j][k].By;
						cube[i][j][k].Bz = cube2[i][j][k].Bz;
					}

			for (long long int i = 0; i < N_3 + 2; i++)
				for (long long int j = 0; j < N_3 + 2; j++)
					for (long long int k = 0; k < N_3 + 2; k++)
					{
						compensator[i][j][k].Ex = compensator2[i][j][k].Ex;
						compensator[i][j][k].Ey = compensator2[i][j][k].Ey;
						compensator[i][j][k].Ez = compensator2[i][j][k].Ez;
						compensator[i][j][k].Bx = compensator2[i][j][k].Bx;
						compensator[i][j][k].By = compensator2[i][j][k].By;
						compensator[i][j][k].Bz = compensator2[i][j][k].Bz;
					}
			// Изменяем значения на границе
			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{
					cube[N_3 + 1][i][j].Ex = cube[1][i][j].Ex;
					cube[N_3 + 1][i][j].Ey = cube[1][i][j].Ey;
					cube[N_3 + 1][i][j].Ez = cube[1][i][j].Ez;
					cube[N_3 + 1][i][j].Bx = cube[1][i][j].Bx;
					cube[N_3 + 1][i][j].By = cube[1][i][j].By;
					cube[N_3 + 1][i][j].Bz = cube[1][i][j].Bz;

					cube[0][i][j].Ex = cube[N_3][i][j].Ex;
					cube[0][i][j].Ey = cube[N_3][i][j].Ey;
					cube[0][i][j].Ez = cube[N_3][i][j].Ez;
					cube[0][i][j].Bx = cube[N_3][i][j].Bx;
					cube[0][i][j].By = cube[N_3][i][j].By;
					cube[0][i][j].Bz = cube[N_3][i][j].Bz;
				}
			}

			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{
					cube[i][N_3 + 1][j].Ex = cube[i][1][j].Ex;
					cube[i][N_3 + 1][j].Ey = cube[i][1][j].Ey;
					cube[i][N_3 + 1][j].Ez = cube[i][1][j].Ez;
					cube[i][N_3 + 1][j].Bx = cube[i][1][j].Bx;
					cube[i][N_3 + 1][j].By = cube[i][1][j].By;
					cube[i][N_3 + 1][j].Bz = cube[i][1][j].Bz;

					cube[i][0][j].Ex = cube[i][N_3][j].Ex;
					cube[i][0][j].Ey = cube[i][N_3][j].Ey;
					cube[i][0][j].Ez = cube[i][N_3][j].Ez;
					cube[i][0][j].Bx = cube[i][N_3][j].Bx;
					cube[i][0][j].By = cube[i][N_3][j].By;
					cube[i][0][j].Bz = cube[i][N_3][j].Bz;
				}
			}

			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{

					cube[i][j][N_3 + 1].Ex = cube[i][j][1].Ex;
					cube[i][j][N_3 + 1].Ey = cube[i][j][1].Ey;
					cube[i][j][N_3 + 1].Ez = cube[i][j][1].Ez;
					cube[i][j][N_3 + 1].Bx = cube[i][j][1].Bx;
					cube[i][j][N_3 + 1].By = cube[i][j][1].By;
					cube[i][j][N_3 + 1].Bz = cube[i][j][1].Bz;

					cube[i][j][0].Ex = cube[i][j][N_3].Ex;
					cube[i][j][0].Ey = cube[i][j][N_3].Ey;
					cube[i][j][0].Ez = cube[i][j][N_3].Ez;
					cube[i][j][0].Bx = cube[i][j][N_3].Bx;
					cube[i][j][0].By = cube[i][j][N_3].By;
					cube[i][j][0].Bz = cube[i][j][N_3].Bz;
				}
			}
		}
	}
	else
	{
		for (long long int it = 0; it < Nt; it++)
		{
			if (it % 100 == 0)
				cout << it << endl;
			#pragma omp parallel for collapse(3)
			for (long long int i = 1; i < N_3 + 1; i++)
				for (long long int j = 1; j < N_3 + 1; j++)
					for (long long int k = 1; k < N_3 + 1; k++)
					{
						Update_cell_central<ftype>(cube, cube2, dt_x, dt_y, dt_z, dt, dx, ab.first, it, i, j, k);
					}

			for (long long int i = 1; i < N_3 + 1; i++)
				for (long long int j = 1; j < N_3 + 1; j++)
					for (long long int k = 1; k < N_3 + 1; k++)
					{
						cube[i][j][k].Ex = cube2[i][j][k].Ex;
						cube[i][j][k].Ey = cube2[i][j][k].Ey;
						cube[i][j][k].Ez = cube2[i][j][k].Ez;
						cube[i][j][k].Bx = cube2[i][j][k].Bx;
						cube[i][j][k].By = cube2[i][j][k].By;
						cube[i][j][k].Bz = cube2[i][j][k].Bz;
					}

			// Изменяем значения на границе
			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{
					cube[N_3 + 1][i][j].Ex = cube[1][i][j].Ex;
					cube[N_3 + 1][i][j].Ey = cube[1][i][j].Ey;
					cube[N_3 + 1][i][j].Ez = cube[1][i][j].Ez;
					cube[N_3 + 1][i][j].Bx = cube[1][i][j].Bx;
					cube[N_3 + 1][i][j].By = cube[1][i][j].By;
					cube[N_3 + 1][i][j].Bz = cube[1][i][j].Bz;

					cube[0][i][j].Ex = cube[N_3][i][j].Ex;
					cube[0][i][j].Ey = cube[N_3][i][j].Ey;
					cube[0][i][j].Ez = cube[N_3][i][j].Ez;
					cube[0][i][j].Bx = cube[N_3][i][j].Bx;
					cube[0][i][j].By = cube[N_3][i][j].By;
					cube[0][i][j].Bz = cube[N_3][i][j].Bz;
				}
			}

			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{
					cube[i][N_3 + 1][j].Ex = cube[i][1][j].Ex;
					cube[i][N_3 + 1][j].Ey = cube[i][1][j].Ey;
					cube[i][N_3 + 1][j].Ez = cube[i][1][j].Ez;
					cube[i][N_3 + 1][j].Bx = cube[i][1][j].Bx;
					cube[i][N_3 + 1][j].By = cube[i][1][j].By;
					cube[i][N_3 + 1][j].Bz = cube[i][1][j].Bz;

					cube[i][0][j].Ex = cube[i][N_3][j].Ex;
					cube[i][0][j].Ey = cube[i][N_3][j].Ey;
					cube[i][0][j].Ez = cube[i][N_3][j].Ez;
					cube[i][0][j].Bx = cube[i][N_3][j].Bx;
					cube[i][0][j].By = cube[i][N_3][j].By;
					cube[i][0][j].Bz = cube[i][N_3][j].Bz;
				}
			}

			for (long long int i = 0; i < N_3 + 2; i++)
			{
				for (long long int j = 0; j < N_3 + 2; j++)
				{

					cube[i][j][N_3 + 1].Ex = cube[i][j][1].Ex;
					cube[i][j][N_3 + 1].Ey = cube[i][j][1].Ey;
					cube[i][j][N_3 + 1].Ez = cube[i][j][1].Ez;
					cube[i][j][N_3 + 1].Bx = cube[i][j][1].Bx;
					cube[i][j][N_3 + 1].By = cube[i][j][1].By;
					cube[i][j][N_3 + 1].Bz = cube[i][j][1].Bz;

					cube[i][j][0].Ex = cube[i][j][N_3].Ex;
					cube[i][j][0].Ey = cube[i][j][N_3].Ey;
					cube[i][j][0].Ez = cube[i][j][N_3].Ez;
					cube[i][j][0].Bx = cube[i][j][N_3].Bx;
					cube[i][j][0].By = cube[i][j][N_3].By;
					cube[i][j][0].Bz = cube[i][j][N_3].Bz;
				}
			}
		}
	}

	vector<vector<vector<Component<ftype>>>> analitical_cube(N_3 + 1, vector<vector<Component<ftype>>>(N_3 + 1, vector<Component<ftype>>(N_3 + 1)));

	double t = (double)Nt * dt;
	for (long long int i = 0; i < N_3 + 1; i++)
	{
		for (long long int j = 0; j < N_3 + 1; j++)
		{
			for (long long int k = 0; k < N_3 + 1; k++)
			{
				ftype x = (ftype)(ab.first + (ftype)i * dx);

				analitical_cube[i][j][k].Bz = sin(x - t);
				analitical_cube[i][j][k].Ey = sin(x - t);
			}
		}
	}

	// СЧИТАЕМ ОШИБКУ

	double abs_err = 0.0, rel_err = 0.0, rel_err_new = 0.0;
	double smape_err = 0.0;
	double temp, temp1 = 0.0;
	int i_max = 0, j_max = 0, k_max = 0;
	int i_max_rel = 0, j_max_rel = 0, k_max_rel = 0;
	double critical_val;

	if (sizeof(cube[0][0][0].Ey) == 2)
		critical_val = 0.1;
	else critical_val = 0.0001;

	for (long long int i = 0; i < N_3; i++)
	{
		for (long long int j = 0; j < N_3; j++)
		{
			for (long long int k = 0; k < N_3; k++)
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
		}
	}
	abs_err_return = abs_err;
	rel_err_return = rel_err * (double)100;
	smape_err_return = smape_err / (double)(6 * N_3 * N_3 * N_3);

	return 0.0;
}

template <class ftype>
void Dependence_error_from_step(int N, ftype T, ftype dt_start, ftype dt_finish, int count_step, string type_sum)
{
	int N_new;

	ftype dt;
	ftype step = (ftype) ((dt_finish - dt_start) /count_step);

	pair<ftype, ftype> ab(0.0, 2.0 * M_PI), cd(0.0, 2.0 * M_PI), fg(0.0, 2.0 * M_PI);
	ftype dx = (ab.second - ab.first) / (ftype)N, dy = (cd.second - cd.first) / (ftype)N, dz = (fg.second - fg.first) / (ftype)N;

	vector<double> vec_abs_err(count_step), vec_rel_err(count_step), vec_smape_err(count_step);

	ftype constanta = (ftype)(dx * dx / dt_start); // dx * dx = constanta * dt  =>  N = (ab.second - ab.first) / dx

	double abs_err_return, rel_err_return, smape_err_return;
	for (int i = 0; i < count_step; i++)
	{
		cout << "it = " << i << endl;
		dt = dt_start + step * i;
		N_new = (ab.second - ab.first) / sqrt(constanta * dt);
		Check_Curant(N_new, T, dt);
		Maxwells_central_solution_for_depended<ftype>(N_new, T, dt, type_sum, abs_err_return, rel_err_return, smape_err_return);
		vec_abs_err[i] = abs_err_return; vec_rel_err[i] = rel_err_return; vec_smape_err[i] = smape_err_return;
	}

	// вывод векторов в файл
	ofstream numb1("dependence error from dt.csv");

	if (sizeof(dt) == 2)		numb1 << "Тип данных" << "half" << endl;
	if (sizeof(dt) == 4)		numb1 << "Тип данных" << "float" << endl;
	if (sizeof(dt) == 8)		numb1 << "Тип данных" << "double" << endl;
	numb1 << "Central operator" << endl;
	numb1 << "Тип суммирования " << type_sum << endl;
	numb1 << "T =" << T << endl;
	numb1 << "N = " << N << endl;
	numb1 << "dt_start = " << dt_start << endl;
	numb1 << "dt_finish = " << dt_finish << endl;
	numb1 << "dx*dx/dt = " << constanta << endl << endl;
	numb1 << "count_step = " << count_step << endl << endl;

	for (int s = 0; s < count_step; s++)
	{
		numb1 << dt_start + step * (double)s << ";" << vec_abs_err[s] << ";" << vec_rel_err[s] << ";" << vec_smape_err[s] << endl;
	}
	numb1 << endl << endl;

	numb1.close();
}

template <class ftype>
void Dependence_error_from_step_reverse(int N, ftype T, ftype dt_start, ftype dt_finish, int count_step, string type_sum)
{
	int N_new;

	ftype dt;
	ftype step = (ftype)((dt_start - dt_finish) / count_step);

	pair<ftype, ftype> ab(0.0, 2.0 * M_PI), cd(0.0, 2.0 * M_PI), fg(0.0, 2.0 * M_PI);
	ftype dx = (ab.second - ab.first) / (ftype)N, dy = (cd.second - cd.first) / (ftype)N, dz = (fg.second - fg.first) / (ftype)N;

	vector<double> vec_abs_err(count_step), vec_rel_err(count_step), vec_smape_err(count_step);

	ftype constanta = (ftype)(dx * dx / dt_start); // dx * dx = constanta * dt  =>  N = (ab.second - ab.first) / dx
	cout << "REVERSE" << endl;
	double abs_err_return, rel_err_return, smape_err_return;
	for (int i = 0; i < count_step; i++)
	{
		cout << "it = " << i << endl;
		dt = dt_finish + step * i;
		N_new = (ab.second - ab.first) / sqrt(constanta * dt);
		Check_Curant(N_new, T, dt);
		Maxwells_central_solution_for_depended<ftype>(N_new, T, dt, type_sum, abs_err_return, rel_err_return, smape_err_return);
		vec_abs_err[i] = abs_err_return; vec_rel_err[i] = rel_err_return; vec_smape_err[i] = smape_err_return;
	}
	cout << dt<<endl;
	//вывод векторов в файл
	ofstream numb1("dependence error from dt reverse.csv");

	if (sizeof(dt) == 2)		numb1 << "Тип данных" << "half" << endl;
	if (sizeof(dt) == 4)		numb1 << "Тип данных" << "float" << endl;
	if (sizeof(dt) == 8)		numb1 << "Тип данных" << "double" << endl;
	numb1 << "Central operator" << endl;
	numb1 << "Тип суммирования " << type_sum << endl;
	numb1 << "T =" << T << endl;
	numb1 << "N = " << N << endl;
	numb1 << "dt_start = " << dt_start << endl;
	numb1 << "dt_finish = " << dt_finish << endl;
	numb1 << "dx*dx/dt = " << constanta << endl;
	numb1 << "count_step = " << count_step << endl << endl;

	for (int s = 0; s < count_step; s++)
	{
		numb1 << dt_finish + step * (double)s << ";" << vec_abs_err[s] << ";" << vec_rel_err[s] << ";" << vec_smape_err[s] << endl;
	}
	numb1 << endl << endl;

	numb1.close();
}
