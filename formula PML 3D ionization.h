#pragma once
#define _USE_MATH_DEFINES

#include "type data.h"
#include "cmath"

#define w 2.35e15 //c^-1
#define t_0 3.33333e-15 //c
#define w_a 4.13e16 // c^-1
#define Ng 0.016

using namespace std;

template <typename type_data>
class data3d;

template <class ftype>
ftype moduleE(data3d<Component<ftype>>& cube, int i, int j, int k) {

	return sqrt(cube(i, j, k).Ex * cube(i, j, k).Ex + cube(i, j, k).Ey * cube(i, j, k).Ey + cube(i, j, k).Ez * cube(i, j, k).Ez);
}

template <class ftype>
ftype probability_ionization(data3d<Component<ftype>>& cube, int i, int j, int k) {

	//probability of hydrogen ionization
	if (moduleE(cube, i, j, k) < 2.0e-5) {
		return (ftype)0.0;
	}
	ftype coeff = (ftype)1.;
	ftype t_moduleE = coeff / moduleE(cube, i, j, k);
	ftype temp;

	const double w_a_t_0 = w_a * t_0;
	temp = (ftype)(4. * w_a_t_0)* (t_moduleE * (ftype)exp(-(ftype)(2. / 3.) * t_moduleE));
	
	return temp;
}

template <class ftype>
void Update_electric_field_three_dimen_ionization(data3d<Component<ftype>>& cube, data3d<ftype>& Ne,
	ftype dt, ftype dt_x, ftype dt_y, ftype dt_z, int i, int j, int k)
{
	ftype tEx, tEy, tEz;

	tEx = (cube(i, j + 1, k).Bz - cube(i, j, k).Bz) * (dt_y)
		-(cube(i, j, k + 1).By - cube(i, j, k).By) * (dt_z)+cube(i, j, k).Ex
		- (ftype)(4. * M_PI) * dt * cube(i, j, k).Jx;

	tEy = (cube(i, j, k + 1).Bx - cube(i, j, k).Bx) * (dt_z)
		-(cube(i + 1, j, k).Bz - cube(i, j, k).Bz) * (dt_x)+cube(i, j, k).Ey
		- (ftype)(4. * M_PI) * dt * cube(i, j, k).Jy;

	tEz = (cube(i + 1, j, k).By - cube(i, j, k).By) * (dt_x)
		-(cube(i, j + 1, k).Bx - cube(i, j, k).Bx) * (dt_y)+cube(i, j, k).Ez
		- (ftype)(4. * M_PI) * dt * cube(i, j, k).Jz;

	cube(i, j, k).Ex = tEx;
	cube(i, j, k).Ey = tEy;
	cube(i, j, k).Ez = tEz;
}

template <class ftype>
void Update_electric_currents_three_dimen_ionization(data3d<Component<ftype>>& cube, data3d<ftype>& Ne,
	ftype dt, ftype dt_x, ftype dt_y, ftype dt_z, int i, int j, int k)
{
	ftype tJx, tJy, tJz, tNe, tProb_ion;
	const double w_t_0 = w * t_0;
	tJx = cube(i, j, k).Jx + (ftype)(w_t_0 * w_t_0 / (4. * M_PI)) * Ne(i, j, k) * dt * cube(i, j, k).Ex;
	tJy = cube(i, j, k).Jy + (ftype)(w_t_0 * w_t_0 / (4. * M_PI)) * Ne(i, j, k) * dt * cube(i, j, k).Ey;
	tJz = cube(i, j, k).Jz + (ftype)(w_t_0 * w_t_0 / (4. * M_PI)) * Ne(i, j, k) * dt * cube(i, j, k).Ez;


	tProb_ion = probability_ionization(cube, i, j, k);

	//tNe = Ne(i, j, k) * ((ftype)1. - dt * tProb_ion) + dt * tProb_ion * (ftype)Ng;

	ftype k1 = tProb_ion * (ftype)(Ng - Ne(i, j, k));
	ftype k2 = tProb_ion * (ftype)(Ng - Ne(i, j, k) - dt / 2. * k1);
	ftype k3 = tProb_ion * (ftype)(Ng - Ne(i, j, k) - dt / 2. * k2);
	ftype k4 = tProb_ion * (ftype)(Ng - Ne(i, j, k) - dt * k3);

	tNe = Ne(i, j, k) + dt / 6. * (k1 + 2. * k2 + 2. * k3 + k4);

	
	if (moduleE(cube, i, j, k) > 1.0e-4) {
		//cout << tProb_ion << " "<< (ftype)exp(-(ftype)(2. / 3.) / moduleE(cube, i, j, k)) <<endl;
	}
	cube(i, j, k).Jx = tJx;
	cube(i, j, k).Jy = tJy;
	cube(i, j, k).Jz = tJz;

	Ne(i, j, k) = tNe;
}


