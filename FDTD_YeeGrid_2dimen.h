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

template <class ftype>
void Initializing_FDTD_two_dimen(vector<vector<Component<ftype>>>& cube, int Nx, int Ny, const ftype dx, const ftype dt,
	pair <ftype, ftype> ab)
{
	ftype x1, x2;
	for (long long int i = 0; i < Nx; i++)
	{
		for (long long int j = 0; j < Ny; j++)
		{
			x1 = (ftype)(ab.first + dx * (ftype)i);
			x2 = (ftype)(ab.first + dx * (ftype)i + dx / 2.0);


			cube[i][j].Ex = (ftype)0.0;
			cube[i][j].Ey = sin(x1);
			cube[i][j].Ez = (ftype)0.0;
			cube[i][j].Bx = (ftype)0.0;
			cube[i][j].By = (ftype)0.0;
			cube[i][j].Bz = sin(x2 - dt / (ftype)2.0);
		}
	}
	//Заполняем отдельно границу куба
	for (long long int i = 0; i < Ny + 2; i++) {
		cube[Nx][i].Ex = cube[0][i].Ex;
		cube[Nx][i].Ey = cube[0][i].Ey;
		cube[Nx][i].Ez = cube[0][i].Ez;
		cube[Nx][i].Bx = cube[0][i].Bx;
		cube[Nx][i].By = cube[0][i].By;
		cube[Nx][i].Bz = cube[0][i].Bz;

		cube[Nx + 1][i].Ex = cube[1][i].Ex;
		cube[Nx + 1][i].Ey = cube[1][i].Ey;
		cube[Nx + 1][i].Ez = cube[1][i].Ez;
		cube[Nx + 1][i].Bx = cube[1][i].Bx;
		cube[Nx + 1][i].By = cube[1][i].By;
		cube[Nx + 1][i].Bz = cube[1][i].Bz;
	}
	for (long long int i = 0; i < Nx + 2; i++) {
		cube[i][Ny].Ex = cube[i][0].Ex;
		cube[i][Ny].Ey = cube[i][0].Ey;
		cube[i][Ny].Ez = cube[i][0].Ez;
		cube[i][Ny].Bx = cube[i][0].Bx;
		cube[i][Ny].By = cube[i][0].By;
		cube[i][Ny].Bz = cube[i][0].Bz;

		cube[i][Ny + 1].Ex = cube[i][1].Ex;
		cube[i][Ny + 1].Ey = cube[i][1].Ey;
		cube[i][Ny + 1].Ez = cube[i][1].Ez;
		cube[i][Ny + 1].Bx = cube[i][1].Bx;
		cube[i][Ny + 1].By = cube[i][1].By;
		cube[i][Ny + 1].Bz = cube[i][1].Bz;
	}
}

template <class ftype>
void Initializing_FDTD_two_dimen_Gauss(vector<vector<Component<ftype>>>& cube, int Nx, int Ny, const ftype dx, const ftype dt,
	pair <ftype, ftype> ab)
{
	ftype x1, x2, y0, y1, y2; 
	y0 = (ab.second - ab.first) / 2.0;
	ftype delta = M_PI / 3.0;
	for (long long int i = 0; i < Nx; i++)
	{
		for (long long int j = 0; j < Ny; j++)
		{
			x1 = (ftype)(ab.first + dx * (ftype)i + (ftype)0.5 * dx);
			x2 = (ftype)(ab.first + dx * (ftype)i);

			y1 = (ftype)(ab.first + dx * (ftype)j + (ftype)0.5 * dx);
			y2 = (ftype)(ab.first + dx * (ftype)j);



			
		}
	}
	//Заполняем отдельно границу куба
	for (long long int i = 0; i < Ny + 2; i++) {
		cube[Nx][i].Ex = cube[0][i].Ex;
		cube[Nx][i].Ey = cube[0][i].Ey;
		cube[Nx][i].Ez = cube[0][i].Ez;
		cube[Nx][i].Bx = cube[0][i].Bx;
		cube[Nx][i].By = cube[0][i].By;
		cube[Nx][i].Bz = cube[0][i].Bz;

		cube[Nx + 1][i].Ex = cube[1][i].Ex;
		cube[Nx + 1][i].Ey = cube[1][i].Ey;
		cube[Nx + 1][i].Ez = cube[1][i].Ez;
		cube[Nx + 1][i].Bx = cube[1][i].Bx;
		cube[Nx + 1][i].By = cube[1][i].By;
		cube[Nx + 1][i].Bz = cube[1][i].Bz;
	}
	for (long long int i = 0; i < Nx + 2; i++) {
		cube[i][Ny].Ex = cube[i][0].Ex;
		cube[i][Ny].Ey = cube[i][0].Ey;
		cube[i][Ny].Ez = cube[i][0].Ez;
		cube[i][Ny].Bx = cube[i][0].Bx;
		cube[i][Ny].By = cube[i][0].By;
		cube[i][Ny].Bz = cube[i][0].Bz;

		cube[i][Ny + 1].Ex = cube[i][1].Ex;
		cube[i][Ny + 1].Ey = cube[i][1].Ey;
		cube[i][Ny + 1].Ez = cube[i][1].Ez;
		cube[i][Ny + 1].Bx = cube[i][1].Bx;
		cube[i][Ny + 1].By = cube[i][1].By;
		cube[i][Ny + 1].Bz = cube[i][1].Bz;
	}
}

template <class ftype>
void Update_electric_field_two_dimen(vector<vector<Component<ftype>>>& cube,
	ftype _dt_x, ftype _dt_y,
	long long int i, long long int j)
{
	ftype tEx, tEy, tEz;

	tEx = (cube[i][j + 1].Bz - cube[i][j].Bz) * (_dt_y) + cube[i][j].Ex;
	tEy = -(cube[i + 1][j].Bz - cube[i][j].Bz) * (_dt_x) + cube[i][j].Ey;
	tEz = (cube[i + 1][j].By - cube[i][j].By) * (_dt_x)
		-(cube[i][j + 1].Bx - cube[i][j].Bx) * (_dt_y) + cube[i][j].Ez;

	cube[i][j].Ex = tEx;
	cube[i][j].Ey = tEy;
	cube[i][j].Ez = tEz;
}

template <class ftype>
void Update_magnetic_field_two_dimen(vector<vector<Component<ftype>>>& cube,
	ftype _dt_x, ftype _dt_y,
	long long int i, long long int j)
{
	ftype tBx, tBy, tBz;

	tBx = (cube[i][j - 1].Ez - cube[i][j].Ez) * (_dt_y) + cube[i][j].Bx;
	tBy = (cube[i][j].Ez - cube[i - 1][j].Ez) * (_dt_x) + cube[i][j].By;
	tBz = (cube[i - 1][j].Ey - cube[i][j].Ey) * (_dt_x)
		+(cube[i][j].Ex - cube[i][j - 1].Ex) * (_dt_y) + cube[i][j].Bz;

	cube[i][j].Bx = tBx;
	cube[i][j].By = tBy;
	cube[i][j].Bz = tBz;
}

template <class ftype>
double FDTD_two_dimen(long long int _N, ftype T, ftype _dt, string type_sum)
{
	setlocale(LC_ALL, "Russian");
	long long int N_3 = _N; 	long long int Nt = T / _dt;

	vector<vector<Component<ftype>>> cube(N_3 + 2, vector<Component<ftype>>(N_3 + 2));
	/*vector<vector<vector<Component<ftype>>>> compensator(N_3 + 2, vector<vector<Component<ftype>>>(N_3 + 2, vector<Component<ftype>>(N_3 + 2)));
	vector<vector<vector<Component<ftype>>>> compensator2(N_3 + 2, vector<vector<Component<ftype>>>(N_3 + 2, vector<Component<ftype>>(N_3 + 2)));*/
	pair<ftype, ftype> ab(0.0, 2.0 * M_PI), cd(0.0, 2.0 * M_PI);
	ftype dt = _dt, dx = (ab.second - ab.first) / (ftype)N_3, dy = (cd.second - cd.first) / (ftype)N_3;
	ftype dt_x = dt / (dx), dt_y = dt / (dy);

	vector<double> vec_energy(Nt + 1);

	Initializing_FDTD_two_dimen_Gauss<ftype>(cube, N_3, N_3, dx, dt, ab);

	if (type_sum == "Kahan")
	{

	}
	else
	{
		vec_energy[0] = CalculateEnergyPML(cube, N_3, N_3, 0);
		for (long long int it = 0; it < Nt; it++)
		{
			if (it % 100 == 0) cout << it << endl;

#pragma omp parallel for collapse(2)
			for (long long int i = 1; i < N_3 + 1; i++)
				for (long long int j = 1; j < N_3 + 1; j++)
					Update_electric_field_two_dimen<ftype>(cube, dt_x, dt_y, i, j);

#pragma omp parallel for collapse(2)
			for (long long int i = 1; i < N_3 + 1; i++)
				for (long long int j = 1; j < N_3 + 1; j++)
					Update_magnetic_field_two_dimen<ftype>(cube, dt_x, dt_y, i, j);

			// Border conditions
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

			vec_energy[it + 1] = CalculateEnergyPML(cube, N_3, N_3, 0);
		}
	}

	// ������� ������

		//вывод вектора энергии в файл
	ofstream numb_energy("Energy_graph_without_PML.csv");
	for (int s = 0; s < Nt; s++)
	{
		numb_energy << dt * (double)(s + 1) << ";" << vec_energy[s] << endl;
	}
	numb_energy << endl << endl;
	numb_energy.close();


	ofstream numb2("graph_solution2.csv");
	numb2 << ";" << "y" << endl;
	numb2 << "x" << ";" << ";";
	for (int i = 1; i < N_3 + 1; i++)
	{
		numb2 << dx * (ftype)i << ";";
	}
	numb2 << endl;
	for (int i = 1; i < N_3 + 1; i++)
	{
		ftype y = dy * (ftype)i;
		numb2 << ";" << y << ";";
		for (int j = 1; j < N_3 + 1; j++)
		{
			numb2 << cube[j][i].Ey << ";";
		}
		numb2 << endl;
	}

	numb2.close();

	return 0.0;
}