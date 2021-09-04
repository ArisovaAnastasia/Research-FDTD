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
void Initializing_central_two_dimen(vector<vector<Component<ftype>>>& cube, long long int _N, const ftype _dx,
	pair <ftype, ftype> ab, pair <ftype, ftype> cd)
{
	ftype x;
	for (long long int i = 0; i < _N; i++)
	{
		for (long long int j = 0; j < _N; j++)
		{
			x = (ftype)(ab.first + _dx * (ftype)i);

			cube[i][j].Ex = (ftype)0.0;
			cube[i][j].Ey = sin(x);
			cube[i][j].Ez = (ftype)0.0;
			cube[i][j].Bx = (ftype)0.0;
			cube[i][j].By = (ftype)0.0;
			cube[i][j].Bz = sin(x);
		}
	}
	//��������� �������� ������� �������� 
	for (long long int i = 0; i < _N + 2; i++) {
			cube[_N][i].Ex = cube[0][i].Ex;
			cube[_N][i].Ey = cube[0][i].Ey;
			cube[_N][i].Ez = cube[0][i].Ez;
			cube[_N][i].Bx = cube[0][i].Bx;
			cube[_N][i].By = cube[0][i].By;
			cube[_N][i].Bz = cube[0][i].Bz;

			cube[_N + 1][i].Ex = cube[1][i].Ex;
			cube[_N + 1][i].Ey = cube[1][i].Ey;
			cube[_N + 1][i].Ez = cube[1][i].Ez;
			cube[_N + 1][i].Bx = cube[1][i].Bx;
			cube[_N + 1][i].By = cube[1][i].By;
			cube[_N + 1][i].Bz = cube[1][i].Bz;
	}
	for (long long int i = 0; i < _N + 2; i++) {
			cube[i][_N].Ex = cube[i][0].Ex;
			cube[i][_N].Ey = cube[i][0].Ey;
			cube[i][_N].Ez = cube[i][0].Ez;
			cube[i][_N].Bx = cube[i][0].Bx;
			cube[i][_N].By = cube[i][0].By;
			cube[i][_N].Bz = cube[i][0].Bz;

			cube[i][_N + 1].Ex = cube[i][1].Ex;
			cube[i][_N + 1].Ey = cube[i][1].Ey;
			cube[i][_N + 1].Ez = cube[i][1].Ez;
			cube[i][_N + 1].Bx = cube[i][1].Bx;
			cube[i][_N + 1].By = cube[i][1].By;
			cube[i][_N + 1].Bz = cube[i][1].Bz;
	}
}

template <class ftype>
void Update_cell_central_two_dimen(vector<vector<Component<ftype>>>& cube, vector<vector<Component<ftype>>>& cube2,
	ftype _dt_x, ftype _dt_y, ftype dt, ftype dx, ftype ab_first, long long int it,
	long long int i, long long int j)
{
	ftype tEx, tEy, tEz, tBx, tBy, tBz;
	ftype  J = (ftype)0;

	tBx = (cube[i][j - 1].Ez - cube[i][j + 1].Ez) * (_dt_y) / (ftype)2.0 + cube[i][j].Bx;
	tBy = (cube[i + 1][j].Ez - cube[i - 1][j].Ez) * (_dt_x) / (ftype)2.0 + cube[i][j].By;
	tEz = (cube[i + 1][j].By - cube[i - 1][j].By) * (_dt_x) / (ftype)2.0
		- (cube[i][j + 1].Bx - cube[i][j - 1].Bx) * (_dt_y) / (ftype)2.0 + cube[i][j].Ez - J * dt;
	
	cube2[i][j].Bx = tBx;
	cube2[i][j].By = tBy;
	cube2[i][j].Ez = tEz;

	tEx = (cube[i][j + 1].Bz - cube[i][j - 1].Bz) * (_dt_y) / (ftype)2.0 + cube[i][j].Ex - J * dt;
	tEy = - (cube[i + 1][j].Bz - cube[i - 1][j].Bz) * (_dt_x) / (ftype)2.0 + cube[i][j].Ey - J * dt;
	tBz = (cube[i - 1][j].Ey - cube[i + 1][j].Ey) * (_dt_x) / (ftype)2.0
		+ (cube[i][j + 1].Ex - cube[i][j - 1].Ex) * (_dt_y) / (ftype)2.0 + cube[i][j].Bz;
	

	cube2[i][j].Ex = tEx;
	cube2[i][j].Ey = tEy;
	cube2[i][j].Bz = tBz;
}

template <class ftype>
double Maxwells_central_solution_two_dimen(long long int _N, ftype T, ftype _dt, string type_sum)
{
	setlocale(LC_ALL, "Russian");

	long long int N_3 = _N; 	long long int Nt = T / _dt;

	vector<vector<Component<ftype>>> cube(N_3 + 2, vector<Component<ftype>>(N_3 + 2));
	vector<vector<Component<ftype>>> cube2(N_3 + 2, vector<Component<ftype>>(N_3 + 2));	
	/*vector<vector<vector<Component<ftype>>>> compensator(N_3 + 2, vector<vector<Component<ftype>>>(N_3 + 2, vector<Component<ftype>>(N_3 + 2)));
	vector<vector<vector<Component<ftype>>>> compensator2(N_3 + 2, vector<vector<Component<ftype>>>(N_3 + 2, vector<Component<ftype>>(N_3 + 2)));*/

	pair<ftype, ftype> ab(0.0, 2.0 * M_PI), cd(0.0, 2.0 * M_PI);
	ftype dt = _dt, dx = (ab.second - ab.first) / (ftype)N_3, dy = (cd.second - cd.first) / (ftype)N_3;
	ftype dt_x = dt / (dx), dt_y = dt / (dy);

	Initializing_central_two_dimen<ftype>(cube, N_3, dx, ab, cd);

	if (type_sum == "Kahan")
	{

	}
	else
	{
		for (long long int it = 0; it < Nt; it++)
		{
			if (it % 100 == 0)
				cout << it << endl;
#pragma omp parallel for collapse(2)
			for (long long int i = 1; i < N_3 + 1; i++)
				for (long long int j = 1; j < N_3 + 1; j++)
						Update_cell_central_two_dimen<ftype>(cube, cube2, dt_x, dt_y, dt, dx, ab.first, it, i, j);
					

			for (long long int i = 1; i < N_3 + 1; i++)
				for (long long int j = 1; j < N_3 + 1; j++)
				{
					cube[i][j].Ex = cube2[i][j].Ex;
					cube[i][j].Ey = cube2[i][j].Ey;
					cube[i][j].Ez = cube2[i][j].Ez;
					cube[i][j].Bx = cube2[i][j].Bx;
					cube[i][j].By = cube2[i][j].By;
					cube[i][j].Bz = cube2[i][j].Bz;
				}

			//��������� �������� ������� �������� 
			for (long long int i = 0; i < _N + 2; i++) {
				cube[0][i].Ex = cube[_N][i].Ex;
				cube[0][i].Ey = cube[_N][i].Ey;
				cube[0][i].Ez = cube[_N][i].Ez;
				cube[0][i].Bx = cube[_N][i].Bx;
				cube[0][i].By = cube[_N][i].By;
				cube[0][i].Bz = cube[_N][i].Bz;

				cube[_N + 1][i].Ex = cube[1][i].Ex;
				cube[_N + 1][i].Ey = cube[1][i].Ey;
				cube[_N + 1][i].Ez = cube[1][i].Ez;
				cube[_N + 1][i].Bx = cube[1][i].Bx;
				cube[_N + 1][i].By = cube[1][i].By;
				cube[_N + 1][i].Bz = cube[1][i].Bz;
			}
			for (long long int i = 0; i < _N + 2; i++) {
				cube[i][0].Ex = cube[i][_N].Ex;
				cube[i][0].Ey = cube[i][_N].Ey;
				cube[i][0].Ez = cube[i][_N].Ez;
				cube[i][0].Bx = cube[i][_N].Bx;
				cube[i][0].By = cube[i][_N].By;
				cube[i][0].Bz = cube[i][_N].Bz;

				cube[i][_N + 1].Ex = cube[i][1].Ex;
				cube[i][_N + 1].Ey = cube[i][1].Ey;
				cube[i][_N + 1].Ez = cube[i][1].Ez;
				cube[i][_N + 1].Bx = cube[i][1].Bx;
				cube[i][_N + 1].By = cube[i][1].By;
				cube[i][_N + 1].Bz = cube[i][1].Bz;
			}
		}
	}

	vector<vector<Component<ftype>>> analitical_cube(N_3 + 1, vector<Component<ftype>>(N_3 + 1));

	double t = (double)Nt * dt;
	for (long long int i = 0; i < N_3 + 1; i++)
	{
		for (long long int j = 0; j < N_3 + 1; j++)
		{
			ftype x = (ftype)(ab.first + (ftype)i * dx);

			analitical_cube[i][j].Bz = sin(x - t);
			analitical_cube[i][j].Ey = sin(x - t);
		}
	}

	// ������� ������

	double abs_err = 0.0, rel_err = 0.0, rel_err_new = 0.0;
	double smape_err = 0.0;
	double temp, temp1 = 0.0;
	int i_max = 0, j_max = 0;
	int i_max_rel = 0, j_max_rel = 0;
	double critical_val;

	if (sizeof(cube[0][0].Ey) == 2)
		critical_val = 0.1;
	else critical_val = 0.0001;

	for (long long int i = 0; i < N_3; i++)
	{
		for (long long int j = 0; j < N_3; j++)
		{
			temp1 = abs((double)analitical_cube[i][j].Bz - (double)cube[i][j].Bz);
			if (temp1 > abs_err)
			{
				abs_err = temp1;
				i_max = i; j_max = j;
			}
			temp1 = abs((double)analitical_cube[i][j].Ey - (double)cube[i][j].Ey);
			if (temp1 > abs_err)
			{
				abs_err = temp1;
				i_max = i; j_max = j;
			}

			if ((double)analitical_cube[i][j].Bz > critical_val)
				rel_err = max(rel_err, abs((double)analitical_cube[i][j].Ey - (double)cube[i][j].Ey) / (double)analitical_cube[i][j].Ey);
			if ((double)analitical_cube[i][j].Bz > critical_val)
				rel_err = max(rel_err, abs((double)analitical_cube[i][j].Bz - (double)cube[i][j].Bz) / (double)analitical_cube[i][j].Bz);

			temp = (double)2 * abs((double)analitical_cube[i][j].Bz - (double)cube[i][j].Bz) /
				(abs((double)analitical_cube[i][j].Bz) + abs((double)cube[i][j].Bz));
			smape_err += temp;
			temp = (double)2 * abs((double)analitical_cube[i][j].Ey - (double)cube[i][j].Ey) /
				(abs((double)analitical_cube[i][j].Ey) + abs((double)cube[i][j].Ey));
			smape_err += temp;
		}
	}
	// ����� ������ �������� ������������� ������
	rel_err_new = max(rel_err_new, abs((double)analitical_cube[i_max][j_max].Ey - (double)cube[i_max][j_max].Ey / (double)analitical_cube[i_max][j_max].Ey));
	rel_err_new = max(rel_err_new, abs((double)analitical_cube[i_max][j_max].Bz - (double)cube[i_max][j_max].Bz / (double)analitical_cube[i_max][j_max].Bz));


	cout << endl << "i_max = " << i_max << "  j_max = " << j_max << endl;
	cout << "Absolute error " << abs_err << endl;
	cout << "Relative error " << rel_err * (double)100 << endl;
	cout << "Relative error new version " << rel_err_new * (double)100 << endl;
	cout << "Smape error " << smape_err / (double)(6 * N_3 * N_3 * N_3) << endl;
	return 0.0;
}

template <class ftype>
double Maxwells_central_solution_two_dimen_graph(long long int _N, ftype T, ftype _dt, string type_sum)
{
	setlocale(LC_ALL, "Russian");

	long long int N_3 = _N; 	long long int Nt = T / _dt;

	vector<vector<Component<ftype>>> cube(N_3 + 2, vector<Component<ftype>>(N_3 + 2));
	vector<vector<Component<ftype>>> cube2(N_3 + 2, vector<Component<ftype>>(N_3 + 2));

	pair<ftype, ftype> ab(0.0, 2.0 * M_PI), cd(0.0, 2.0 * M_PI);
	ftype dt = _dt, dx = (ab.second - ab.first) / (ftype)N_3, dy = (cd.second - cd.first) / (ftype)N_3;
	ftype dt_x = dt / (dx), dt_y = dt / (dy);

	ofstream numb1("graph_solution_2.csv");
	numb1 << ";" << "y" << endl;
	numb1 << "x" << ";" << ";";
	for (int i = 0; i < N_3 + 1; i++)
	{
		numb1 << dx * (ftype)i << ";";
	}
	numb1 << endl;

	Initializing_central_two_dimen<ftype>(cube, N_3, dx, ab, cd);

	if (type_sum == "Kahan")
	{

	}
	else
	{
		for (long long int it = 0; it < Nt; it++)
		{
			if (it % 100 == 0)
				cout << it << endl;
#pragma omp parallel for collapse(2)
			for (long long int i = 1; i < N_3 + 1; i++)
				for (long long int j = 1; j < N_3 + 1; j++)
					Update_cell_central_two_dimen<ftype>(cube, cube2, dt_x, dt_y, dt, dx, ab.first, it, i, j);


			for (long long int i = 1; i < N_3 + 1; i++)
				for (long long int j = 1; j < N_3 + 1; j++)
				{
					cube[i][j].Ex = cube2[i][j].Ex;
					cube[i][j].Ey = cube2[i][j].Ey;
					cube[i][j].Ez = cube2[i][j].Ez;
					cube[i][j].Bx = cube2[i][j].Bx;
					cube[i][j].By = cube2[i][j].By;
					cube[i][j].Bz = cube2[i][j].Bz;
				}

			//��������� �������� ������� �������� 
			for (long long int i = 0; i < _N + 2; i++) {
				cube[0][i].Ex = cube[_N][i].Ex;
				cube[0][i].Ey = cube[_N][i].Ey;
				cube[0][i].Ez = cube[_N][i].Ez;
				cube[0][i].Bx = cube[_N][i].Bx;
				cube[0][i].By = cube[_N][i].By;
				cube[0][i].Bz = cube[_N][i].Bz;

				cube[_N + 1][i].Ex = cube[1][i].Ex;
				cube[_N + 1][i].Ey = cube[1][i].Ey;
				cube[_N + 1][i].Ez = cube[1][i].Ez;
				cube[_N + 1][i].Bx = cube[1][i].Bx;
				cube[_N + 1][i].By = cube[1][i].By;
				cube[_N + 1][i].Bz = cube[1][i].Bz;
			}
			for (long long int i = 0; i < _N + 2; i++) {
				cube[i][0].Ex = cube[i][_N].Ex;
				cube[i][0].Ey = cube[i][_N].Ey;
				cube[i][0].Ez = cube[i][_N].Ez;
				cube[i][0].Bx = cube[i][_N].Bx;
				cube[i][0].By = cube[i][_N].By;
				cube[i][0].Bz = cube[i][_N].Bz;

				cube[i][_N + 1].Ex = cube[i][1].Ex;
				cube[i][_N + 1].Ey = cube[i][1].Ey;
				cube[i][_N + 1].Ez = cube[i][1].Ez;
				cube[i][_N + 1].Bx = cube[i][1].Bx;
				cube[i][_N + 1].By = cube[i][1].By;
				cube[i][_N + 1].Bz = cube[i][1].Bz;
			}
			
			ftype t = dt * (ftype)it;
			numb1 << ";" << t << ";";
			for (int j = 0; j < N_3 + 1; j++)
			{
				numb1 << cube[dx * (ftype)j][0].Ey << ";";
			}
			numb1 << endl;			
		}
	}
	numb1.close();

	ofstream numb2("graph_solution.csv");
	numb2 << ";" << "y" << endl;
	numb2 << "x" << ";"<< ";";
	for (int i = 0; i < N_3 + 1; i++)
	{
		numb2 << dx * (ftype)i<< ";";
	}
	numb2 << endl;
	for (int i = 0; i < N_3 + 1; i++)
	{
		ftype y = dy * (ftype)i;
		numb2 << ";" << y<< ";";
		for (int j = 0; j < N_3 + 1; j++)
		{
			numb2 << cube[j][i].Ey<< ";";
		}
		numb2 << endl;
	}

	numb2.close();




	vector<vector<Component<ftype>>> analitical_cube(N_3 + 1, vector<Component<ftype>>(N_3 + 1));

	double t = (double)Nt * dt;
	for (long long int i = 0; i < N_3 + 1; i++)
	{
		for (long long int j = 0; j < N_3 + 1; j++)
		{
			ftype x = (ftype)(ab.first + (ftype)i * dx);

			analitical_cube[i][j].Bz = sin(x - t);
			analitical_cube[i][j].Ey = sin(x - t);
		}
	}
	// ������� ������

	double abs_err = 0.0, rel_err = 0.0, rel_err_new = 0.0;
	double smape_err = 0.0;
	double temp, temp1 = 0.0;
	int i_max = 0, j_max = 0;
	int i_max_rel = 0, j_max_rel = 0;
	double critical_val;

	if (sizeof(cube[0][0].Ey) == 2)
		critical_val = 0.1;
	else critical_val = 0.0001;

	for (long long int i = 0; i < N_3; i++)
	{
		for (long long int j = 0; j < N_3; j++)
		{
			temp1 = abs((double)analitical_cube[i][j].Bz - (double)cube[i][j].Bz);
			if (temp1 > abs_err)
			{
				abs_err = temp1;
				i_max = i; j_max = j;
			}
			temp1 = abs((double)analitical_cube[i][j].Ey - (double)cube[i][j].Ey);
			if (temp1 > abs_err)
			{
				abs_err = temp1;
				i_max = i; j_max = j;
			}

			if ((double)analitical_cube[i][j].Bz > critical_val)
				rel_err = max(rel_err, abs((double)analitical_cube[i][j].Ey - (double)cube[i][j].Ey) / (double)analitical_cube[i][j].Ey);
			if ((double)analitical_cube[i][j].Bz > critical_val)
				rel_err = max(rel_err, abs((double)analitical_cube[i][j].Bz - (double)cube[i][j].Bz) / (double)analitical_cube[i][j].Bz);

			temp = (double)2 * abs((double)analitical_cube[i][j].Bz - (double)cube[i][j].Bz) /
				(abs((double)analitical_cube[i][j].Bz) + abs((double)cube[i][j].Bz));
			smape_err += temp;
			temp = (double)2 * abs((double)analitical_cube[i][j].Ey - (double)cube[i][j].Ey) /
				(abs((double)analitical_cube[i][j].Ey) + abs((double)cube[i][j].Ey));
			smape_err += temp;
		}
	}
	// ����� ������ �������� ������������� ������
	rel_err_new = max(rel_err_new, abs((double)analitical_cube[i_max][j_max].Ey - (double)cube[i_max][j_max].Ey / (double)analitical_cube[i_max][j_max].Ey));
	rel_err_new = max(rel_err_new, abs((double)analitical_cube[i_max][j_max].Bz - (double)cube[i_max][j_max].Bz / (double)analitical_cube[i_max][j_max].Bz));


	cout << endl << "i_max = " << i_max << "  j_max = " << j_max << endl;
	cout << "Absolute error " << abs_err << endl;
	cout << "Relative error " << rel_err * (double)100 << endl;
	cout << "Relative error new version " << rel_err_new * (double)100 << endl;
	cout << "Smape error " << smape_err / (double)(6 * N_3 * N_3 * N_3) << endl;
	return 0.0;
}
