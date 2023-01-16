#include <CL/sycl.hpp>
#include "type data.h"
#include <cmath> 

#define IsPrint false 
#define IsPrint2 false 

template <class ftypePML>
ftypePML distanceH(int N, int delta, int i)
{
	ftypePML dist = 0.0;

	if (delta == 0)
		return 0.0;

	if (i < delta + 1)
		dist = (ftypePML)(delta + 1 - i);

	if (i > N + delta)
		dist = (ftypePML)(i - delta - N) - (ftypePML)0.5;

	return dist / (ftypePML)delta;
}

template <class ftypePML>
ftypePML distanceE(int N, int delta, int i)
{ //���� ���������� ������������ PML-����
	ftypePML dist = 0.0;

	if (delta == 0)
		return 0.0;

	if (i < delta + 1)
		dist = (ftypePML)(delta + 1 - i) - (ftypePML)0.5;

	if (i > N + delta)
		dist = (ftypePML)(i - delta - N);

	return dist / (ftypePML)delta;
}

template <class ftype_init, class ftypePML>
void Initializing_cube_sycl_half_version(Fields<ftype_init>& cube, SplitFields<ftypePML>& cube_split,
	int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	for (int i = 0; i < Nx + 2 * delta_x + 2; i++)
		for (int j = 0; j < Ny + 2 * delta_y + 2; j++)
			for (int k = 0; k < Nz + 2 * delta_z + 2; k++) {
				cube_split.Exy(i, j, k) = (ftypePML)0.0;
				cube_split.Exz(i, j, k) = (ftypePML)0.0;
				cube_split.Eyx(i, j, k) = (ftypePML)0.0;
				cube_split.Eyz(i, j, k) = (ftypePML)0.0;
				cube_split.Ezx(i, j, k) = (ftypePML)0.0;
				cube_split.Ezy(i, j, k) = (ftypePML)0.0;

				cube_split.Bxy(i, j, k) = (ftypePML)0.0;
				cube_split.Bxz(i, j, k) = (ftypePML)0.0;
				cube_split.Byx(i, j, k) = (ftypePML)0.0;
				cube_split.Byz(i, j, k) = (ftypePML)0.0;
				cube_split.Bzx(i, j, k) = (ftypePML)0.0;
				cube_split.Bzy(i, j, k) = (ftypePML)0.0;

				cube.Ex(i, j, k) = cl::sycl::detail::float2Half(0.0);
				cube.Ey(i, j, k) = cl::sycl::detail::float2Half(0.0);
				cube.Ez(i, j, k) = cl::sycl::detail::float2Half(0.0);
				cube.Bx(i, j, k) = cl::sycl::detail::float2Half(0.0);
				cube.By(i, j, k) = cl::sycl::detail::float2Half(0.0);
				cube.Bz(i, j, k) = cl::sycl::detail::float2Half(0.0);
			}
}

//template <class ftype_init>
//void Initializing_Sigma_sycl_half_version(Sigma<ftype_init>& Sigma,
//	int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z, double n,
//	double sigma_x, double sigma_y)
//{
//	double var_Sigma_max_x = sigma_x;
//	double var_Sigma_max_y = sigma_y;
//	double var_Sigma_max_z = sigma_y;
//
//
//
//	for (int i = 1; i < Nx + 2 * delta_x + 1; i++) {
//		for (int j = 1; j < Ny + 2 * delta_y + 1; j++) {
//			for (int k = 1; k < Nz + 2 * delta_z + 1; k++) {
//				Sigma.sigmaEx(i, j, k) = float_to_half(var_Sigma_max_x * powf(distanceE<double>(Nx, delta_x, i), n));
//				Sigma.sigmaBx(i, j, k) = float_to_half(var_Sigma_max_x * powf(distanceH<double>(Nx, delta_x, i), n));
//
//				Sigma.sigmaEy(i, j, k) = float_to_half(var_Sigma_max_y * powf(distanceE<double>(Ny, delta_y, j), n));
//				Sigma.sigmaBy(i, j, k) = float_to_half(var_Sigma_max_y * powf(distanceH<double>(Ny, delta_y, j), n));
//
//				Sigma.sigmaEz(i, j, k) = float_to_half(var_Sigma_max_z * powf(distanceE<double>(Nz, delta_z, k), n));
//				Sigma.sigmaBz(i, j, k) = float_to_half(var_Sigma_max_z * powf(distanceH<double>(Nz, delta_z, k), n));
//			}
//		}
//	}
//}

//template <class ftype_init>
//void Initializing_Sigma_Coeff_sycl_half_version(Coefficient<ftype_init>& Coeff, Sigma<ftype_init>& Sigma, int Nx, int Ny, int Nz,
//	int delta_x, int delta_y, int delta_z, double n,
//	double sigma_x, double sigma_y, double dt) {
//
//	float var_Sigma_max_x = sigma_x;
//	float var_Sigma_max_y = sigma_y;
//	float var_Sigma_max_z = sigma_y;
//
//	float sigmaEx, sigmaEy, sigmaEz, sigmaBx, sigmaBy, sigmaBz;
//
//	float _Exy1, _Exz1, _Ezx1, _Ezy1, _Eyx1, _Eyz1;
//	float _Bxy1, _Bxz1, _Bzx1, _Bzy1, _Byx1, _Byz1;
//
//	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
//		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
//			for (int k = 1; k < Nz + 2 * delta_z + 1; k++) {
//				sigmaEx = var_Sigma_max_x * powf(distanceE<float>(Nx, delta_x, i), n);
//				sigmaBx = var_Sigma_max_x * powf(distanceH<float>(Nx, delta_x, i), n);
//
//				sigmaEy = var_Sigma_max_y * powf(distanceE<float>(Ny, delta_y, j), n);
//				sigmaBy = var_Sigma_max_y * powf(distanceH<float>(Ny, delta_y, j), n);
//
//				sigmaEz = var_Sigma_max_z * powf(distanceE<float>(Nz, delta_z, k), n);
//				sigmaBz = var_Sigma_max_z * powf(distanceH<float>(Nz, delta_z, k), n);
//
//				Sigma.sigmaEx(i, j, k) = float_to_half(sigmaEx);
//				Sigma.sigmaBx(i, j, k) = float_to_half(sigmaBx);
//
//				Sigma.sigmaEy(i, j, k) = float_to_half(sigmaEy);
//				Sigma.sigmaBy(i, j, k) = float_to_half(sigmaBy);
//
//				Sigma.sigmaEz(i, j, k) = float_to_half(sigmaEz);
//				Sigma.sigmaBz(i, j, k) = float_to_half(sigmaBz);
//
//				if ((i >= delta_x + 1) && (i < Nx + delta_x + 1) && (j >= delta_y + 1) && (j < Ny + delta_y + 1) && (k >= delta_z + 1) && (k < Nz + delta_z + 1))
//				{
//
//				}
//				else {
//					_Exy1 = exp(-dt * sigmaEy);
//					_Exz1 = exp(-dt * sigmaEz);
//					_Eyx1 = exp(-dt * sigmaEx);
//					_Eyz1 = exp(-dt * sigmaEz);
//					_Ezx1 = exp(-dt * sigmaEx);
//					_Ezy1 = exp(-dt * sigmaEy);
//
//					_Bxy1 = exp(-dt * sigmaBy);
//					_Bxz1 = exp(-dt * sigmaBz);
//					_Byx1 = exp(-dt * sigmaBx);
//					_Byz1 = exp(-dt * sigmaBz);
//					_Bzx1 = exp(-dt * sigmaBx);
//					_Bzy1 = exp(-dt * sigmaBy);
//
//					Coeff.Exy1(i, j, k) = float_to_half(_Exy1);
//					Coeff.Exz1(i, j, k) = float_to_half(_Exz1);
//					Coeff.Eyx1(i, j, k) = float_to_half(_Eyx1);
//					Coeff.Eyz1(i, j, k) = float_to_half(_Eyz1);
//					Coeff.Ezx1(i, j, k) = float_to_half(_Ezx1);
//					Coeff.Ezy1(i, j, k) = float_to_half(_Ezy1);
//
//					Coeff.Bxy1(i, j, k) = float_to_half(_Bxy1);
//					Coeff.Bxz1(i, j, k) = float_to_half(_Bxz1);
//					Coeff.Byx1(i, j, k) = float_to_half(_Byx1);
//					Coeff.Byz1(i, j, k) = float_to_half(_Byz1);
//					Coeff.Bzx1(i, j, k) = float_to_half(_Bzx1);
//					Coeff.Bzy1(i, j, k) = float_to_half(_Bzy1);
//
//					if (abs(sigmaEx) > FLT_EPSILON) {
//						Coeff.Eyx2(i, j, k) = float_to_half(1.0 / sigmaEx - _Eyx1 / sigmaEx);
//						Coeff.Ezx2(i, j, k) = float_to_half(1.0 / sigmaEx - _Ezx1 / sigmaEx);
//					}
//					else {
//						Coeff.Eyx2(i, j, k) = float_to_half(dt);
//						Coeff.Ezx2(i, j, k) = float_to_half(dt);
//					}
//					if (abs(sigmaEy) > FLT_EPSILON) {
//						Coeff.Exy2(i, j, k) = float_to_half(1.0 / sigmaEy - _Exy1 / sigmaEy);
//						Coeff.Ezy2(i, j, k) = float_to_half(1.0 / sigmaEy - _Ezy1 / sigmaEy);
//					}
//					else {
//						Coeff.Exy2(i, j, k) = float_to_half(dt);
//						Coeff.Ezy2(i, j, k) = float_to_half(dt);
//					}
//					if (abs(sigmaEz) > FLT_EPSILON)
//					{
//						Coeff.Exz2(i, j, k) = float_to_half(1.0 / sigmaEz - _Exz1 / sigmaEz);
//						Coeff.Eyz2(i, j, k) = float_to_half(1.0 / sigmaEz - _Eyz1 / sigmaEz);
//					}
//					else {
//						Coeff.Exz2(i, j, k) = float_to_half(dt);
//						Coeff.Eyz2(i, j, k) = float_to_half(dt);
//					}
//					if (abs(sigmaBx) > FLT_EPSILON)
//					{
//						Coeff.Byx2(i, j, k) = float_to_half(1.0 / sigmaBx - _Byx1 / sigmaBx);
//						Coeff.Bzx2(i, j, k) = float_to_half(1.0 / sigmaBx - _Bzx1 / sigmaBx);
//					}
//					else {
//						Coeff.Byx2(i, j, k) = float_to_half(dt);
//						Coeff.Bzx2(i, j, k) = float_to_half(dt);
//					}
//					if (abs(sigmaBy) > FLT_EPSILON)
//					{
//						Coeff.Bxy2(i, j, k) = float_to_half(1.0 / sigmaBy - _Bxy1 / sigmaBy);
//						Coeff.Bzy2(i, j, k) = float_to_half(1.0 / sigmaBy - _Bzy1 / sigmaBy);
//					}
//					else {
//						Coeff.Bxy2(i, j, k) = float_to_half(dt);
//						Coeff.Bzy2(i, j, k) = float_to_half(dt);
//					}
//					if (abs(sigmaBz) > FLT_EPSILON)
//					{
//						Coeff.Bxz2(i, j, k) = float_to_half(1.0 / sigmaBz - _Bxz1 / sigmaBz);
//						Coeff.Byz2(i, j, k) = float_to_half(1.0 / sigmaBz - _Byz1 / sigmaBz);
//					}
//					else {
//						Coeff.Bxz2(i, j, k) = float_to_half(dt);
//						Coeff.Byz2(i, j, k) = float_to_half(dt);
//					}
//				}
//			}
//}

template <class ftype_init, class ftype_sigma>
void Initializing_Sigma_Coeff_sycl_half_version_2_0(Coefficient<ftype_init>& Coeff, Sigma<ftype_sigma>& Sigma, int Nx, int Ny, int Nz,
	int delta_x, int delta_y, int delta_z, double n,
	ftype_sigma sigma_x, ftype_sigma sigma_y, double dt) {

	ftype_sigma var_Sigma_max_x = sigma_x;
	ftype_sigma var_Sigma_max_y = sigma_y;
	ftype_sigma var_Sigma_max_z = sigma_y;

	ofstream numb_sgm("sigma_new.txt"), numb_coef("coeff_new.csv");


	ftype_sigma sigmaEx, sigmaEy, sigmaEz, sigmaBx, sigmaBy, sigmaBz;

	ftype_init Exy1, Exz1, Ezx1, Ezy1, Eyx1, Eyz1;
	ftype_init Bxy1, Bxz1, Bzx1, Bzy1, Byx1, Byz1;

	if (sizeof(ftype_init) == 2) {

		for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
			for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
				for (int k = 1; k < Nz + 2 * delta_z + 1; k++) {
					sigmaEx = var_Sigma_max_x * powf(distanceE<ftype_init>(Nx, delta_x, i), n);
					sigmaBx = var_Sigma_max_x * powf(distanceH<ftype_init>(Nx, delta_x, i), n);

					sigmaEy = var_Sigma_max_y * powf(distanceE<ftype_init>(Ny, delta_y, j), n);
					sigmaBy = var_Sigma_max_y * powf(distanceH<ftype_init>(Ny, delta_y, j), n);

					sigmaEz = var_Sigma_max_z * powf(distanceE<ftype_init>(Nz, delta_z, k), n);
					sigmaBz = var_Sigma_max_z * powf(distanceH<ftype_init>(Nz, delta_z, k), n);

					Sigma.sigmaEx(i, j, k) = sigmaEx;
					Sigma.sigmaBx(i, j, k) = sigmaBx;

					Sigma.sigmaEy(i, j, k) = sigmaEy;
					Sigma.sigmaBy(i, j, k) = sigmaBy;

					Sigma.sigmaEz(i, j, k) = sigmaEz;
					Sigma.sigmaBz(i, j, k) = sigmaBz;

					//numb_sgm << sigmaEx <<" "<< sigmaEy << " " << sigmaEz << " " << sigmaBx << " " << sigmaBy << " " << sigmaBz<<" ";

					if ((i >= delta_x + 1) && (i < Nx + delta_x + 1) && (j >= delta_y + 1) && (j < Ny + delta_y + 1) && (k >= delta_z + 1) && (k < Nz + delta_z + 1))
					{

					}
					else {
						Exy1 = exp(-dt * sigmaEy);
						Exz1 = exp(-dt * sigmaEz);
						Eyx1 = exp(-dt * sigmaEx);
						Eyz1 = exp(-dt * sigmaEz);
						Ezx1 = exp(-dt * sigmaEx);
						Ezy1 = exp(-dt * sigmaEy);

						Bxy1 = exp(-dt * sigmaBy);
						Bxz1 = exp(-dt * sigmaBz);
						Byx1 = exp(-dt * sigmaBx);
						Byz1 = exp(-dt * sigmaBz);
						Bzx1 = exp(-dt * sigmaBx);
						Bzy1 = exp(-dt * sigmaBy);

						Coeff.Exy1(i, j, k) = cl::sycl::detail::float2Half(Exy1);
						Coeff.Exz1(i, j, k) = cl::sycl::detail::float2Half(Exz1);
						Coeff.Eyx1(i, j, k) = cl::sycl::detail::float2Half(Eyx1);
						Coeff.Eyz1(i, j, k) = cl::sycl::detail::float2Half(Eyz1);
						Coeff.Ezx1(i, j, k) = cl::sycl::detail::float2Half(Ezx1);
						Coeff.Ezy1(i, j, k) = cl::sycl::detail::float2Half(Ezy1);

						Coeff.Bxy1(i, j, k) = cl::sycl::detail::float2Half(Bxy1);
						Coeff.Bxz1(i, j, k) = cl::sycl::detail::float2Half(Bxz1);
						Coeff.Byx1(i, j, k) = cl::sycl::detail::float2Half(Byx1);
						Coeff.Byz1(i, j, k) = cl::sycl::detail::float2Half(Byz1);
						Coeff.Bzx1(i, j, k) = cl::sycl::detail::float2Half(Bzx1);
						Coeff.Bzy1(i, j, k) = cl::sycl::detail::float2Half(Bzy1);

						if (sigmaEx != (ftype_sigma)0.0) {
							Coeff.Eyx2(i, j, k) = cl::sycl::detail::float2Half(1.0 / sigmaEx - Eyx1 / sigmaEx);
							Coeff.Ezx2(i, j, k) = cl::sycl::detail::float2Half(1.0 / sigmaEx - Ezx1 / sigmaEx);
						}
						else {
							Coeff.Eyx2(i, j, k) = cl::sycl::detail::float2Half(dt);
							Coeff.Ezx2(i, j, k) = cl::sycl::detail::float2Half(dt);
						}
						if (sigmaEy != (ftype_sigma)0.0) {
							Coeff.Exy2(i, j, k) = cl::sycl::detail::float2Half(1.0 / sigmaEy - Exy1 / sigmaEy);
							Coeff.Ezy2(i, j, k) = cl::sycl::detail::float2Half(1.0 / sigmaEy - Ezy1 / sigmaEy);
						}
						else {
							Coeff.Exy2(i, j, k) = cl::sycl::detail::float2Half(dt);
							Coeff.Ezy2(i, j, k) = cl::sycl::detail::float2Half(dt);
						}
						if (sigmaEz != (ftype_sigma)0.0)
						{
							Coeff.Exz2(i, j, k) = cl::sycl::detail::float2Half(1.0 / sigmaEz - Exz1 / sigmaEz);
							Coeff.Eyz2(i, j, k) = cl::sycl::detail::float2Half(1.0 / sigmaEz - Eyz1 / sigmaEz);
						}
						else {
							Coeff.Exz2(i, j, k) = cl::sycl::detail::float2Half(dt);
							Coeff.Eyz2(i, j, k) = cl::sycl::detail::float2Half(dt);
						}
						if (sigmaBx != (ftype_sigma)0.0)
						{
							Coeff.Byx2(i, j, k) = cl::sycl::detail::float2Half(1.0 / sigmaBx - Byx1 / sigmaBx);
							Coeff.Bzx2(i, j, k) = cl::sycl::detail::float2Half(1.0 / sigmaBx - Bzx1 / sigmaBx);
						}
						else {
							Coeff.Byx2(i, j, k) = cl::sycl::detail::float2Half(dt);
							Coeff.Bzx2(i, j, k) = cl::sycl::detail::float2Half(dt);
						}
						if (sigmaBy != (ftype_sigma)0.0)
						{
							Coeff.Bxy2(i, j, k) = cl::sycl::detail::float2Half(1.0 / sigmaBy - Bxy1 / sigmaBy);
							Coeff.Bzy2(i, j, k) = cl::sycl::detail::float2Half(1.0 / sigmaBy - Bzy1 / sigmaBy);
						}
						else {
							Coeff.Bxy2(i, j, k) = cl::sycl::detail::float2Half(dt);
							Coeff.Bzy2(i, j, k) = cl::sycl::detail::float2Half(dt);
						}
						if (sigmaBz != (ftype_sigma)0.0)
						{
							Coeff.Bxz2(i, j, k) = cl::sycl::detail::float2Half(1.0 / sigmaBz - Bxz1 / sigmaBz);
							Coeff.Byz2(i, j, k) = cl::sycl::detail::float2Half(1.0 / sigmaBz - Byz1 / sigmaBz);
						}
						else {
							Coeff.Bxz2(i, j, k) = cl::sycl::detail::float2Half(dt);
							Coeff.Byz2(i, j, k) = cl::sycl::detail::float2Half(dt);
						}
					}
				}
	}
	else {
		for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
			for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
				for (int k = 1; k < Nz + 2 * delta_z + 1; k++) {
					sigmaEx = var_Sigma_max_x * powf(distanceE<ftype_init>(Nx, delta_x, i), n);
					sigmaBx = var_Sigma_max_x * powf(distanceH<ftype_init>(Nx, delta_x, i), n);

					sigmaEy = var_Sigma_max_y * powf(distanceE<ftype_init>(Ny, delta_y, j), n);
					sigmaBy = var_Sigma_max_y * powf(distanceH<ftype_init>(Ny, delta_y, j), n);

					sigmaEz = var_Sigma_max_z * powf(distanceE<ftype_init>(Nz, delta_z, k), n);
					sigmaBz = var_Sigma_max_z * powf(distanceH<ftype_init>(Nz, delta_z, k), n);

					Sigma.sigmaEx(i, j, k) = sigmaEx;
					Sigma.sigmaBx(i, j, k) = sigmaBx;

					Sigma.sigmaEy(i, j, k) = sigmaEy;
					Sigma.sigmaBy(i, j, k) = sigmaBy;

					Sigma.sigmaEz(i, j, k) = sigmaEz;
					Sigma.sigmaBz(i, j, k) = sigmaBz;

					numb_sgm << sigmaEx << " " << sigmaEy << " " << sigmaEz << " " << sigmaBx << " " << sigmaBy << " " << sigmaBz << " ";

					if ((i >= delta_x + 1) && (i < Nx + delta_x + 1) && (j >= delta_y + 1) && (j < Ny + delta_y + 1) && (k >= delta_z + 1) && (k < Nz + delta_z + 1))
					{

					}
					else {
						Exy1 = exp(-dt * sigmaEy);
						Exz1 = exp(-dt * sigmaEz);
						Eyx1 = exp(-dt * sigmaEx);
						Eyz1 = exp(-dt * sigmaEz);
						Ezx1 = exp(-dt * sigmaEx);
						Ezy1 = exp(-dt * sigmaEy);

						Bxy1 = exp(-dt * sigmaBy);
						Bxz1 = exp(-dt * sigmaBz);
						Byx1 = exp(-dt * sigmaBx);
						Byz1 = exp(-dt * sigmaBz);
						Bzx1 = exp(-dt * sigmaBx);
						Bzy1 = exp(-dt * sigmaBy);

						numb_coef << Exy1 << "; " << Exz1 << " ;" << Ezx1 << "; " << Ezy1 << "; " << Eyx1 << " ;" << Eyz1 << "; "
							<< Bxy1 << " ;" << Bxz1 << " ;" << Bzx1 << " ;" << Bzy1 << " ;" << Byx1 << "; " << Byz1 << " ;";
						Coeff.Exy1(i, j, k) = (ftype_init)Exy1;
						Coeff.Exz1(i, j, k) = (ftype_init)Exz1;
						Coeff.Eyx1(i, j, k) = (ftype_init)Eyx1;
						Coeff.Eyz1(i, j, k) = (ftype_init)Eyz1;
						Coeff.Ezx1(i, j, k) = (ftype_init)Ezx1;
						Coeff.Ezy1(i, j, k) = (ftype_init)Ezy1;

						Coeff.Bxy1(i, j, k) = (ftype_init)Bxy1;
						Coeff.Bxz1(i, j, k) = (ftype_init)Bxz1;
						Coeff.Byx1(i, j, k) = (ftype_init)Byx1;
						Coeff.Byz1(i, j, k) = (ftype_init)Byz1;
						Coeff.Bzx1(i, j, k) = (ftype_init)Bzx1;
						Coeff.Bzy1(i, j, k) = (ftype_init)Bzy1;

						if (Sigma.sigmaEx(i, j, k) != (ftype_sigma)0.0) {
							Coeff.Eyx2(i, j, k) = 1.0 / Sigma.sigmaEx(i, j, k) - Eyx1 / Sigma.sigmaEx(i, j, k);
							Coeff.Ezx2(i, j, k) = 1.0 / Sigma.sigmaEx(i, j, k) - Ezx1 / Sigma.sigmaEx(i, j, k);
						}
						else {
							Coeff.Eyx2(i, j, k) = dt;
							Coeff.Ezx2(i, j, k) = dt;
						}
						if (Sigma.sigmaEy(i, j, k) != (ftype_sigma)0.0) {
							Coeff.Exy2(i, j, k) = 1.0 / Sigma.sigmaEy(i, j, k) - Exy1 / Sigma.sigmaEy(i, j, k);
							Coeff.Ezy2(i, j, k) = 1.0 / Sigma.sigmaEy(i, j, k) - Ezy1 / Sigma.sigmaEy(i, j, k);
						}
						else {
							Coeff.Exy2(i, j, k) = dt;
							Coeff.Ezy2(i, j, k) = dt;
						}
						if (Sigma.sigmaEz(i, j, k) != (ftype_sigma)0.0)
						{
							Coeff.Exz2(i, j, k) = 1.0 / Sigma.sigmaEz(i, j, k) - Exz1 / Sigma.sigmaEz(i, j, k);
							Coeff.Eyz2(i, j, k) = 1.0 / Sigma.sigmaEz(i, j, k) - Eyz1 / Sigma.sigmaEz(i, j, k);
						}
						else {
							Coeff.Exz2(i, j, k) = dt;
							Coeff.Eyz2(i, j, k) = dt;
						}
						if (Sigma.sigmaBx(i, j, k) != (ftype_sigma)0.0)
						{
							Coeff.Byx2(i, j, k) = 1.0 / Sigma.sigmaBx(i, j, k) - Byx1 / Sigma.sigmaBx(i, j, k);
							Coeff.Bzx2(i, j, k) = 1.0 / Sigma.sigmaBx(i, j, k) - Bzx1 / Sigma.sigmaBx(i, j, k);
						}
						else {
							Coeff.Byx2(i, j, k) = dt;
							Coeff.Bzx2(i, j, k) = dt;
						}
						if (Sigma.sigmaBy(i, j, k) != (ftype_sigma)0.0)
						{
							Coeff.Bxy2(i, j, k) = 1.0 / Sigma.sigmaBy(i, j, k) - Bxy1 / Sigma.sigmaBy(i, j, k);
							Coeff.Bzy2(i, j, k) = 1.0 / Sigma.sigmaBy(i, j, k) - Bzy1 / Sigma.sigmaBy(i, j, k);
						}
						else {
							Coeff.Bxy2(i, j, k) = dt;
							Coeff.Bzy2(i, j, k) = dt;
						}
						if (Sigma.sigmaBz(i, j, k) != (ftype_sigma)0.0)
						{
							Coeff.Bxz2(i, j, k) = 1.0 / Sigma.sigmaBz(i, j, k) - Bxz1 / Sigma.sigmaBz(i, j, k);
							Coeff.Byz2(i, j, k) = 1.0 / Sigma.sigmaBz(i, j, k) - Byz1 / Sigma.sigmaBz(i, j, k);
						}
						else {
							Coeff.Bxz2(i, j, k) = dt;
							Coeff.Byz2(i, j, k) = dt;
						}
						numb_coef << Coeff.Exy2(i, j, k) << " ;" << Coeff.Exz2(i, j, k) << " ;" << Coeff.Ezx2(i, j, k) << "; " << Coeff.Ezy2(i, j, k) << " ;" << Coeff.Eyx2(i, j, k) << "; " << Coeff.Eyz2(i, j, k) << " ;"
							<< Coeff.Bxy2(i, j, k) << "; " << Coeff.Bxz2(i, j, k) << " ;" << Coeff.Bzx2(i, j, k) << " ;" << Coeff.Bzy2(i, j, k) << " ;" << Coeff.Byx2(i, j, k) << " ;" << Coeff.Byz2(i, j, k) << "; ";

					}
				}
	}

	numb_sgm.close();
	numb_coef.close();
}


template <class ftype_init>
void Add_Currents_sycl_half_version(ftype_init& Ex, ftype_init& Ey, ftype_init& Ez, ftype_init& Bx, ftype_init& By, ftype_init& Bz,
	int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z, ftype_init dx, ftype_init dy, ftype_init dz, const ftype_init dt, int it,
	std::pair <ftype_init, ftype_init> ab, std::pair <ftype_init, ftype_init> cd, std::pair <ftype_init, ftype_init> fg)
{
	ftype_init x1, x2, y1, y2, z1, z2, t1, t2;

	ftype_init ax_transverse = ab.second * (ftype)3. / (ftype)4.; // поперечное
	ftype_init ay_transverse = cd.second * (ftype)1. / (ftype)6.;
	ftype_init az_transverse = fg.second * (ftype)1. / (ftype)6.;

	//ftype ax_transverse = ab.second * (ftype)1. / (ftype)2.; // поперечное
	//ftype ay_transverse = cd.second * (ftype)1. / (ftype)2.;
	//ftype az_transverse = fg.second * (ftype)1. / (ftype)2.;

	ftype_init tp_x = ab.second / (ftype_init)6.0;
	ftype_init tp_y = cd.second / (ftype_init)6.0;
	ftype_init tp_z = fg.second / (ftype_init)6.0;

	ftype_init lambda_ = (ftype_init)2. * (ftype_init)M_PI / (ftype_init)4.; // wvelength               //?
	ftype_init k_ = (ftype_init)2. * (ftype_init)M_PI / lambda_; // wavenumber = 2 pi/ wavelength
	ftype_init omega_ = (ftype_init)2. * (ftype_init)M_PI / lambda_; // 2 pi c /wavelength 
	ftype_init A = (ftype_init)1.; //amplitude

	t1 = dt * (ftype_init)it;
	t2 = t1 + (ftype_init)0.5 * dt;

	ftype_init x0 = (ab.second - ab.first) / (ftype_init)2.;
	ftype_init y0 = (cd.second - cd.first) / (ftype_init)2.;
	ftype_init z0 = (fg.second - fg.first) / (ftype_init)2.;
	ftype_init t0_x = (ftype_init)3. * tp_x;
	ftype_init t0_y = (ftype_init)3. * tp_y;
	ftype_init t0_z = (ftype_init)3. * tp_z;

	int index_start_x = delta_x + 1, index_start_y = delta_y + 1, index_start_z = delta_z + 1;

	int offset = 5;
	x1 = (ftype_init)(ab.first + dx * (ftype_init)offset + (ftype_init)0.5 * dx);
	x2 = (ftype_init)(ab.first + dx * (ftype_init)offset);

	for (int j = 0; j < Ny; j++)
		for (int k = 0; k < Nz; k++) {

			y1 = (ftype_init)(cd.first + dy * (ftype_init)j + (ftype_init)0.5 * dy);
			y2 = (ftype_init)(cd.first + dy * (ftype_init)j);

			z1 = (ftype_init)(fg.first + dz * (ftype_init)k + (ftype_init)0.5 * dz);
			z2 = (ftype_init)(fg.first + dz * (ftype_init)k);

			// Ez and By has the same y coordinate
			Ey(offset + index_start_x, j + index_start_y, k + index_start_z) += A * exp(-(t1 - t0_x) * (t1 - t0_x) / (tp_x * tp_x))
				* exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse)) *
				sin(omega_ * t1 - k_ * x1); // sin(phase), phase = omega * t - k * x
			Bz(offset + index_start_x, j + index_start_y, k + index_start_z) += A * exp(-(t2 - t0_x) * (t2 - t0_x) / (tp_x * tp_x))
				* exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse)) *
				sin(omega_ * t2 - k_ * x2); // sin(phase), phase = omega * t - k * x
		}

}

template <class ftype_init>
void Calculate_Currents_sycl_half_version(char comp, sycl::queue& q, Fields<ftype_init>& cube, int it, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z,
	ftype_init dx, ftype_init dy, ftype_init dz, ftype_init dt,
	std::pair <ftype_init, ftype_init> ab, std::pair <ftype_init, ftype_init> cd, std::pair <ftype_init, ftype_init> fg)
{
	ftype_init tp_x;
	ftype_init tp_y;
	ftype_init tp_z;

	ftype_init lambda_;
	ftype_init k_ ;
	ftype_init omega_ ;
	ftype_init A ; //amplitude
	ftype_init x0;
	ftype_init y0 ;
	ftype_init z0;
	ftype_init t0_x;
	ftype_init t0_y;
	ftype_init t0_z ;

	q.submit([&](sycl::handler& h) {
		sycl::stream out(4096, 4096, h);
		
		tp_x = ab.second / ftype_init(6.0f);
		tp_y = cd.second / ftype_init(6.0f);
		tp_z = fg.second / ftype_init(6.0f);

		lambda_ = ftype_init(2.0f) * ftype_init((float)M_PI)
			/ ftype_init(4.0f);
		k_ = ftype_init(2.0f) * ftype_init((float)M_PI) / lambda_;
		omega_ = ftype_init(2.0f) * ftype_init((float)M_PI) / lambda_;
		A = ftype_init(10.0f); //amplitude
		x0 = (ab.second - ab.first) / ftype_init(2.0f);
		y0 = (cd.second - cd.first) / ftype_init(2.0f);
		z0 = (fg.second - fg.first) / ftype_init(2.0f);
		t0_x = ftype_init(3.0f) * tp_x;
		t0_y = ftype_init(3.0f) * tp_y;
		t0_z = ftype_init(3.0f) * tp_z;
	}).wait_and_throw();

	int index_start_x = delta_x + 1;
	int index_start_y = delta_y + 1;
	int index_start_z = delta_z + 1;

	int offset = 5;
	float foffset = 5.0;

	switch (comp) {
		case 'X': {
			/*
			q.submit([&](sycl::handler& h) {
				sycl::stream out(4096, 4096, h);
				h.parallel_for(range<1>(2), [=](auto index) {

				half a = half(float(0.1 + 0.2));
				half f = half(float(0.2 + 0.3));
				half c = a + f;

				out << half(a) << " " << a << sycl::endl;
				out << half(f) << " " << f << sycl::endl;
				out << half(c) << " result  " << c << sycl::endl;
					});
				}).wait_and_throw();
				*/
			q.submit([&](sycl::handler& h) {
				sycl::stream out(4096, 4096, h);

				h.parallel_for(range<2>(Ny, Nz), [=](auto index) {

					ftype_init ax_transverse = ab.second * ftype_init(3.0f / 4.0f);
					ftype_init ay_transverse = cd.second * ftype_init(1.0f / 6.0f);
					ftype_init az_transverse = fg.second * ftype_init(1.0f / 6.0f);

					ftype_init t1 = dt * ftype_init((float)it);
					ftype_init t2 = t1 + ftype_init(0.5f) * dt;
					ftype_init x1 = ab.first + dx * ftype_init((float)foffset)
						+ ftype_init(0.5f) * dx;
					ftype_init x2 = ab.first + dx * ftype_init((float)foffset);

					int j = index[0];
					int k = index[1];

					auto kernel_cube = cube;

					ftype_init y1 = (cd.first + dy * ftype_init((float)j) + ftype_init(0.5f) * dy);
					ftype_init y2 = (cd.first + dy * ftype_init((float)j));

					ftype_init z1 = (fg.first + dz * ftype_init((float)k) + ftype_init(0.5f) * dz);
					ftype_init z2 = (fg.first + dz * ftype_init((float)k));

					ftype_init e, b;

					// Ez and By has the same y coordinate
					 
					 e = A * exp(-(t1 - t0_x) * (t1 - t0_x) / (tp_x * tp_x))
						* exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse)) * exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse))
						* sin(omega_ * t1 - k_ * x1);
					 b = A * exp(-(t2 - t0_x) * (t2 - t0_x) / (tp_x * tp_x))
						* exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse)) * exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse)) *
						sin(omega_ * t2 - k_ * x2);

					kernel_cube.Ey(offset + index_start_x, j + index_start_y, k + index_start_z) += e;
					kernel_cube.Bz(offset + index_start_x, j + index_start_y, k + index_start_z) += b;

					//out << float(e) << sycl::endl;
					});
				}).wait_and_throw();
				
			break;
		}

		case 'x': {
			ftype_init ax_transverse = ab.second * (ftype_init)3. / (ftype_init)4.;
			ftype_init ay_transverse = cd.second * (ftype_init)1. / (ftype_init)6.;
			ftype_init az_transverse = fg.second * (ftype_init)1. / (ftype_init)6.;

			ftype_init t1 = dt * (ftype_init)it;
			ftype_init t2 = t1 + (ftype_init)0.5 * dt;
			ftype_init x1 = (ftype_init)(ab.first + dx * (ftype_init)offset + (ftype_init)0.5 * dx);
			ftype_init x2 = (ftype_init)(ab.first + dx * (ftype_init)offset);

			for (int j = 0; j < Ny; j++)
				for (int k = 0; k < Nz; k++) {

					ftype_init y1 = (ftype_init)(cd.first + dy * (ftype_init)j + (ftype_init)0.5 * dy);
					ftype_init y2 = (ftype_init)(cd.first + dy * (ftype_init)j);

					ftype_init z1 = (ftype_init)(fg.first + dz * (ftype_init)k + (ftype_init)0.5 * dz);
					ftype_init z2 = (ftype_init)(fg.first + dz * (ftype_init)k);

					// Ez and By has the same y coordinate
					cube.Ey(Nx + delta_x - offset, j + index_start_y, k + index_start_z) -= A * exp(-(t1 - t0_x) * (t1 - t0_x) / (tp_x * tp_x))
						* exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse)) * exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse))
						* sin(omega_ * t1 - k_ * x1); // sin(phase), phase = omega * t - k * x
					cube.Bz(Nx + delta_x - offset, j + index_start_y, k + index_start_z) += A * exp(-(t2 - t0_x) * (t2 - t0_x) / (tp_x * tp_x))
						* exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse)) * exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse)) *
						sin(omega_ * t2 - k_ * x2); // sin(phase), phase = omega * t - k * x
				}
			break;
		}

		case 'Y': {

			ftype_init ax_transverse = ab.second * (ftype_init)1. / (ftype_init)6.;
			ftype_init ay_transverse = cd.second * (ftype_init)3. / (ftype_init)4.;
			ftype_init az_transverse = fg.second * (ftype_init)1. / (ftype_init)6.;

			ftype_init t1 = dt * (ftype_init)it;
			ftype_init t2 = t1 + (ftype_init)0.5 * dt;
			ftype_init y1 = (ftype_init)(cd.first + dy * (ftype_init)offset + (ftype_init)0.5 * dy);
			ftype_init y2 = (ftype_init)(cd.first + dy * (ftype_init)offset);

			for (int i = 0; i < Nx; i++)
				for (int k = 0; k < Nz; k++) {

					ftype_init x1 = (ftype_init)(ab.first + dx * (ftype_init)i + (ftype_init)0.5 * dx);
					ftype_init x2 = (ftype_init)(ab.first + dx * (ftype_init)i);

					ftype_init z1 = (ftype_init)(fg.first + dz * (ftype_init)k + (ftype_init)0.5 * dz);
					ftype_init z2 = (ftype_init)(fg.first + dz * (ftype_init)k);

					// Ez and By has the same y coordinate
					cube.Ez(i + index_start_x, offset + index_start_y, k + index_start_z) += A * exp(-(t1 - t0_y) * (t1 - t0_y) / (tp_y * tp_y))
						* exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse)) * exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse)) *
						sin(omega_ * t1 - k_ * y1); // sin(phase), phase = omega * t - k * x
					cube.Bx(i + index_start_x, offset + index_start_y, k + index_start_z) += A * exp(-(t2 - t0_y) * (t2 - t0_y) / (tp_y * tp_y))
						* exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse)) * exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse)) *
						sin(omega_ * t2 - k_ * y2); // sin(phase), phase = omega * t - k * x
				}
			break;
		}

		case 'y': {

			ftype_init ax_transverse = ab.second * (ftype_init)1. / (ftype_init)6.;
			ftype_init ay_transverse = cd.second * (ftype_init)3. / (ftype_init)4.;
			ftype_init az_transverse = fg.second * (ftype_init)1. / (ftype_init)6.;

			ftype_init t1 = dt * (ftype_init)it;
			ftype_init t2 = t1 + (ftype_init)0.5 * dt;
			ftype_init y1 = (ftype_init)(cd.first + dy * (ftype_init)offset + (ftype_init)0.5 * dy);
			ftype_init y2 = (ftype_init)(cd.first + dy * (ftype_init)offset);

			for (int i = 0; i < Nx; i++)
				for (int k = 0; k < Nz; k++) {

					ftype_init x1 = (ftype_init)(ab.first + dx * (ftype_init)i + (ftype_init)0.5 * dx);
					ftype_init x2 = (ftype_init)(ab.first + dx * (ftype_init)i);

					ftype_init z1 = (ftype_init)(fg.first + dz * (ftype_init)k + (ftype_init)0.5 * dz);
					ftype_init z2 = (ftype_init)(fg.first + dz * (ftype_init)k);

					// Ez and By has the same y coordinate
					cube.Ez(i + index_start_x, Ny + delta_y - offset, k + index_start_z) -= A * exp(-(t1 - t0_y) * (t1 - t0_y) / (tp_y * tp_y))
						* exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse)) * exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse)) *
						sin(omega_ * t1 - k_ * y1); // sin(phase), phase = omega * t - k * x
					cube.Bx(i + index_start_x, Ny + delta_y - offset, k + index_start_z) += A * exp(-(t2 - t0_y) * (t2 - t0_y) / (tp_y * tp_y))
						* exp(-(z2 - z0) * (z2 - z0) / (az_transverse * az_transverse)) * exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse)) *
						sin(omega_ * t2 - k_ * y2); // sin(phase), phase = omega * t - k * x
				}
			break;
		}

		case 'Z': {
			ftype_init ax_transverse = ab.second * (ftype_init)1. / (ftype_init)6.;
			ftype_init ay_transverse = cd.second * (ftype_init)1. / (ftype_init)6.;
			ftype_init az_transverse = fg.second * (ftype_init)3. / (ftype_init)4.;

			ftype_init t1 = dt * (ftype_init)it;
			ftype_init t2 = t1 + (ftype_init)0.5 * dt;
			ftype_init z1 = (ftype_init)(fg.first + dz * (ftype_init)offset + (ftype_init)0.5 * dz);
			ftype_init z2 = (ftype_init)(fg.first + dz * (ftype_init)offset);

			for (int i = 0; i < Nx; i++)
				for (int j = 0; j < Ny; j++) {

					ftype_init x1 = (ftype_init)(ab.first + dx * (ftype_init)i + (ftype_init)0.5 * dx);
					ftype_init x2 = (ftype_init)(ab.first + dx * (ftype_init)i);

					ftype_init y1 = (ftype_init)(cd.first + dy * (ftype_init)j + (ftype_init)0.5 * dy);
					ftype_init y2 = (ftype_init)(cd.first + dy * (ftype_init)j);

					// Ez and By has the same y coordinate
					cube.Ex(i + index_start_x, j + index_start_y, offset + index_start_z) += A * exp(-(t1 - t0_z) * (t1 - t0_z) / (tp_z * tp_z))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse)) * exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse)) *
						sin(omega_ * t1 - k_ * z1); // sin(phase), phase = omega * t - k * x
					cube.By(i + index_start_x, j + index_start_y, offset + index_start_z) += A * exp(-(t2 - t0_z) * (t2 - t0_z) / (tp_z * tp_z))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse))* exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse)) *
						sin(omega_ * t2 - k_ * z2); // sin(phase), phase = omega * t - k * x
				}
			break;
		}

		case 'z': {
			ftype_init ax_transverse = ab.second * (ftype_init)1. / (ftype_init)6.;
			ftype_init ay_transverse = cd.second * (ftype_init)1. / (ftype_init)6.;
			ftype_init az_transverse = fg.second * (ftype_init)3. / (ftype_init)4.;

			ftype_init t1 = dt * (ftype_init)it;
			ftype_init t2 = t1 + (ftype_init)0.5 * dt;
			ftype_init z1 = (ftype_init)(fg.first + dz * (ftype_init)offset + (ftype_init)0.5 * dz);
			ftype_init z2 = (ftype_init)(fg.first + dz * (ftype_init)offset);

			for (int i = 0; i < Nx; i++)
				for (int j = 0; j < Ny; j++) {

					ftype_init x1 = (ftype_init)(ab.first + dx * (ftype_init)i + (ftype_init)0.5 * dx);
					ftype_init x2 = (ftype_init)(ab.first + dx * (ftype_init)i);

					ftype_init y1 = (ftype_init)(cd.first + dy * (ftype_init)j + (ftype_init)0.5 * dy);
					ftype_init y2 = (ftype_init)(cd.first + dy * (ftype_init)j);

					// Ez and By has the same y coordinate
					cube.Ex(i + index_start_x, j + index_start_y, Nz + delta_z - offset) -= A * exp(-(t1 - t0_z) * (t1 - t0_z) / (tp_z * tp_z))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse)) * exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse)) *
						sin(omega_ * t1 - k_ * z1); // sin(phase), phase = omega * t - k * x
					cube.By(i + index_start_x, j + index_start_y, Nz + delta_z - offset) += A * exp(-(t2 - t0_z) * (t2 - t0_z) / (tp_z * tp_z))
						* exp(-(x2 - x0) * (x2 - x0) / (ax_transverse * ax_transverse)) * exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse)) *
						sin(omega_ * t2 - k_ * z2); // sin(phase), phase = omega * t - k * x
				}
			break;
		}
	}
}

// DOMAIN
template <class ftype_init>
void Update_electric_domain_half_version(Fields<ftype_init> cube,
	ftype_init dt_x, ftype_init dt_y, ftype_init dt_z, int i, int j, int k)
{
	ftype_init tEx, tEy, tEz;

	tEx = (cube.Bz(i, j + 1, k) - cube.Bz(i, j, k)) * (dt_y)
		-(cube.By(i, j, k + 1) - cube.By(i, j, k)) * (dt_z)+cube.Ex(i, j, k);

	tEy = (cube.Bx(i, j, k + 1) - cube.Bx(i, j, k)) * (dt_z)
		-(cube.Bz(i + 1, j, k) - cube.Bz(i, j, k)) * (dt_x)+cube.Ey(i, j, k);

	tEz = (cube.By(i + 1, j, k) - cube.By(i, j, k)) * (dt_x)
		-(cube.Bx(i, j + 1, k) - cube.Bx(i, j, k)) * (dt_y)+cube.Ez(i, j, k);

	cube.Ex(i, j, k) = tEx;
	cube.Ey(i, j, k) = tEy;
	cube.Ez(i, j, k) = tEz;
}

template <class ftype_init>
void Update_magnetic_domain_half_version(Fields<ftype_init> cube,
	ftype_init dt_x, ftype_init dt_y, ftype_init dt_z, int i, int j, int k)
{
	ftype_init tBx, tBy, tBz;

	tBx = -(cube.Ez(i, j, k) - cube.Ez(i, j - 1, k)) * (dt_y)
		+(cube.Ey(i, j, k) - cube.Ey(i, j, k - 1)) * (dt_z)+cube.Bx(i, j, k);

	tBy = -(cube.Ex(i, j, k) - cube.Ex(i, j, k - 1)) * (dt_z)
		+(cube.Ez(i, j, k) - cube.Ez(i - 1, j, k)) * (dt_x)+cube.By(i, j, k);

	tBz = -(cube.Ey(i, j, k) - cube.Ey(i - 1, j, k)) * (dt_x)
		+(cube.Ex(i, j, k) - cube.Ex(i, j - 1, k)) * (dt_y)+cube.Bz(i, j, k);

	cube.Bx(i, j, k) = tBx;
	cube.By(i, j, k) = tBy;
	cube.Bz(i, j, k) = tBz;
}


// PML
template <class ftype, class ftypePML>
void Run_PML_electric(SplitFields<ftypePML> cube_split, Coefficient<ftypePML> Coeff,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z) {
	// corner on z=0
	if (IsPrint) std::cout << "corner on z=0" << std::endl;
	for (int i = 1; i < delta_x + 1; i++)
		for (int j = 1; j < delta_y + 1; j++)
			for (int k = 1; k < delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);

	for (int i = Nx + delta_x + 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < delta_y + 1; j++)
			for (int k = 1; k < delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);

	for (int i = 1; i < delta_x + 1; i++)
		for (int j = Ny + delta_y + 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);

	for (int i = Nx + delta_x + 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = Ny + delta_y + 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	// corner on z=N
	if (IsPrint) std::cout << "corner on z=n" << std::endl;
	for (int i = 1; i < delta_x + 1; i++)
		for (int j = 1; j < delta_y + 1; j++)
			for (int k = Nz + delta_z + 1; k < Nz + 2 * delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);

	for (int i = Nx + delta_x + 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < delta_y + 1; j++)
			for (int k = Nz + delta_z + 1; k < Nz + 2 * delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);

	for (int i = 1; i < delta_x + 1; i++)
		for (int j = Ny + delta_y + 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = Nz + delta_z + 1; k < Nz + 2 * delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);

	for (int i = Nx + delta_x + 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = Ny + delta_y + 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = Nz + delta_z + 1; k < Nz + 2 * delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	//XY  z=0
	if (IsPrint) std::cout << "XY  z=0" << std::endl;
	for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
		for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
			for (int k = 1; k < delta_z; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	//XZ  y=0
	if (IsPrint) std::cout << "XZ  y=0" << std::endl;
	for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
		for (int j = 1; j < delta_y; j++)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	//YZ  x=0
	if (IsPrint) std::cout << "YZ  x=0" << std::endl;
	for (int i = 1; i < delta_x; i++)
		for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	//XY  z=N
	if (IsPrint) std::cout << "XY  z=N" << std::endl;
	for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
		for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
			for (int k = Nz + delta_z + 1; k < Nz + 2 * delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	//XZ  y=N
	if (IsPrint) std::cout << "XZ  y=N" << std::endl;
	for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
		for (int j = Ny + delta_y + 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	//YZ  x=N
	if (IsPrint) std::cout << "YZ  x=N" << std::endl;
	for (int i = Nx + delta_x + 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);


	// колонны вдоль ребер куба

	if (IsPrint) std::cout << "X y=0 z=0" << std::endl;
	//X y=0 z=0
	for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
		for (int j = 1; j < delta_y + 1; j++)
			for (int k = 1; k < delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	if (IsPrint) std::cout << "X y=n z=0" << std::endl;
	//X y=n z=0
	for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
		for (int j = Ny + delta_y + 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	if (IsPrint) std::cout << "X y=n z=n" << std::endl;
	//X y=0 z=n
	for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
		for (int j = 1; j < delta_y + 1; j++)
			for (int k = Nz + delta_z + 1; k < Nz + 2 * delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	if (IsPrint) std::cout << "X y=n z=n" << std::endl;
	//X y=n z=n
	for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
		for (int j = Ny + delta_y + 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = Nz + delta_z + 1; k < Nz + 2 * delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	if (IsPrint) std::cout << "Y x=0 z=0" << std::endl;
	//Y x=0 z=0
	for (int i = 1; i < delta_x + 1; i++)
		for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
			for (int k = 1; k < delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	if (IsPrint) std::cout << "Y x=n z=0" << std::endl;
	//Y x=n z=0
	for (int i = Nx + delta_x + 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
			for (int k = 1; k < delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	if (IsPrint) std::cout << "Y x=0 z=n" << std::endl;
	//Y x=0 z=n
	for (int i = 1; i < delta_x + 1; i++)
		for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
			for (int k = Nz + delta_z + 1; k < Nz + 2 * delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	if (IsPrint) std::cout << "Y x=n z=n" << std::endl;
	//Y x=n z=n
	for (int i = Nx + delta_x + 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
			for (int k = Nz + delta_z + 1; k < Nz + 2 * delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	if (IsPrint) std::cout << "Z x=0 y=0" << std::endl;
	//Z x=0 y=0
	for (int i = 1; i < delta_x + 1; i++)
		for (int j = 1; j < delta_y + 1; j++)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	if (IsPrint) std::cout << "Z x=n y=0" << std::endl;
	//Z x=n y=0
	for (int i = Nx + delta_x + 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < delta_y + 1; j++)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	if (IsPrint) std::cout << "Z x=0 y=n" << std::endl;
	//Z x=0 y=n
	for (int i = 1; i < delta_x + 1; i++)
		for (int j = Ny + delta_y + 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	if (IsPrint) std::cout << "Z x=n y=n" << std::endl;
	//Z x=n y=n
	for (int i = Nx + delta_x + 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = Ny + delta_y + 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
				Update_electric_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
}

template <class ftype, class ftypePML>
void Run_PML_magnetic(SplitFields<ftypePML> cube_split, Coefficient<ftypePML> Coeff,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z) {

	// corner on z=0
	for (int i = 1; i < delta_x + 1; i++)
		for (int j = 1; j < delta_y + 1; j++)
			for (int k = 1; k < delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);

	for (int i = Nx + delta_x + 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < delta_y + 1; j++)
			for (int k = 1; k < delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);

	for (int i = 1; i < delta_x + 1; i++)
		for (int j = Ny + delta_y + 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);

	for (int i = Nx + delta_x + 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = Ny + delta_y + 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	// corner on z=N
	for (int i = 1; i < delta_x + 1; i++)
		for (int j = 1; j < delta_y + 1; j++)
			for (int k = Nz + delta_z + 1; k < Nz + 2 * delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);

	for (int i = Nx + delta_x + 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < delta_y + 1; j++)
			for (int k = Nz + delta_z + 1; k < Nz + 2 * delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);

	for (int i = 1; i < delta_x + 1; i++)
		for (int j = Ny + delta_y + 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = Nz + delta_z + 1; k < Nz + 2 * delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);

	for (int i = Nx + delta_x + 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = Ny + delta_y + 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = Nz + delta_z + 1; k < Nz + 2 * delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);

	//XY  z=0
	for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
		for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
			for (int k = 1; k < delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	//XZ  y=0
	for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
		for (int j = 1; j < delta_y + 1; j++)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	//YZ  x=0
	for (int i = 1; i < delta_x + 1; i++)
		for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	//XY  z=N
	for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
		for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
			for (int k = Nz + delta_z + 2; k < Nz + 2 * delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	//XZ  y=N
	for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
		for (int j = Ny + delta_y + 2; j < Ny + 2 * delta_y + 1; j++)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	//YZ  x=N
	for (int i = Nx + delta_x + 2; i < Nx + 2 * delta_x + 1; i++)
		for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);

	// колонны вдоль ребер куба
	//X y=0 z=0
	for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
		for (int j = 1; j < delta_y + 1; j++)
			for (int k = 1; k < delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	//X y=n z=0
	for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
		for (int j = Ny + delta_y + 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	//X y=0 z=n
	for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
		for (int j = 1; j < delta_y + 1; j++)
			for (int k = Nz + delta_z + 1; k < Nz + 2 * delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	//X y=n z=n
	for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
		for (int j = Ny + delta_y + 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = Nz + delta_z + 1; k < Nz + 2 * delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);

	//Y x=0 z=0
	for (int i = 1; i < delta_x + 1; i++)
		for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
			for (int k = 1; k < delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	//Y x=n z=0
	for (int i = Nx + delta_x + 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
			for (int k = 1; k < delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	//Y x=0 z=n
	for (int i = 1; i < delta_x + 1; i++)
		for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
			for (int k = Nz + delta_z + 1; k < Nz + 2 * delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	//Y x=n z=n
	for (int i = Nx + delta_x + 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
			for (int k = Nz + delta_z + 1; k < Nz + 2 * delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);

	//Z x=0 y=0
	for (int i = 1; i < delta_x + 1; i++)
		for (int j = 1; j < delta_y + 1; j++)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	//Z x=n y=0
	for (int i = Nx + delta_x + 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < delta_y + 1; j++)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	//Z x=0 y=n
	for (int i = 1; i < delta_x + 1; i++)
		for (int j = Ny + delta_y + 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
	//Z x=n y=n
	for (int i = Nx + delta_x + 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = Ny + delta_y + 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
				Update_magnetic_PML<ftype, ftypePML>(cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
}

template <class ftype, class ftypePML>
inline void Update_electric_PML(SplitFields<ftypePML> cube_split, Coefficient<ftypePML> Coeff,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k)
{
	ftypePML tExy, tExz, tEyx, tEyz, tEzx, tEzy;

	//if (IsPrint)
	//	std::cout << "Update_electric_PML: " << " i=" << i << " j=" << j << " k=" << k << std::endl;

	//c1*(A + b*c2)
	tExy = cube_split.Exy(i, j, k) * Coeff.Exy1(i, j, k) + ((cube_split.Bzx(i, j + 1, k) + cube_split.Bzy(i, j + 1, k)) - (cube_split.Bzx(i, j, k) + cube_split.Bzy(i, j, k))) * Coeff.Exy2(i, j, k) * (_1dy);

	tExz = cube_split.Exz(i, j, k) * Coeff.Exz1(i, j, k) - ((cube_split.Byx(i, j, k + 1) + cube_split.Byz(i, j, k + 1)) - (cube_split.Byx(i, j, k) + cube_split.Byz(i, j, k))) * Coeff.Exz2(i, j, k) * (_1dz);

	tEyx = cube_split.Eyx(i, j, k) * Coeff.Eyx1(i, j, k) - ((cube_split.Bzx(i + 1, j, k) + cube_split.Bzy(i + 1, j, k)) - (cube_split.Bzx(i, j, k) + cube_split.Bzy(i, j, k))) * Coeff.Eyx2(i, j, k) * (_1dx);

	tEyz = cube_split.Eyz(i, j, k) * Coeff.Eyz1(i, j, k) + ((cube_split.Bxy(i, j, k + 1) + cube_split.Bxz(i, j, k + 1)) - (cube_split.Bxy(i, j, k) + cube_split.Bxz(i, j, k))) * Coeff.Eyz2(i, j, k) * (_1dz);

	tEzx = cube_split.Ezx(i, j, k) * Coeff.Ezx1(i, j, k) + ((cube_split.Byx(i + 1, j, k) + cube_split.Byz(i + 1, j, k)) - (cube_split.Byx(i, j, k) + cube_split.Byz(i, j, k))) * Coeff.Ezx2(i, j, k) * (_1dx);

	tEzy = cube_split.Ezy(i, j, k) * Coeff.Ezy1(i, j, k) - ((cube_split.Bxy(i, j + 1, k) + cube_split.Bxz(i, j + 1, k)) - (cube_split.Bxy(i, j, k) + cube_split.Bxz(i, j, k))) * Coeff.Ezy2(i, j, k) * (_1dy);


	if (std::abs(tExy) <= 1e-150) {
		tExy = 0.0;
	}
	if (std::abs(tExz) <= 1e-150) {
		tExz = 0.0;
	}
	if (std::abs(tEyx) <= 1e-150) {
		tEyx = 0.0;
	}
	if (std::abs(tEyz) <= 1e-150) {
		tEyz = 0.0;
	}
	if (std::abs(tEzx) <= 1e-150) {
		tEzx = 0.0;
	}
	if (std::abs(tEzy) <= 1e-150) {
		tEzy = 0.0;
	}

	cube_split.Exy(i, j, k) = tExy;
	cube_split.Exz(i, j, k) = tExz;
	cube_split.Eyx(i, j, k) = tEyx;
	cube_split.Eyz(i, j, k) = tEyz;
	cube_split.Ezx(i, j, k) = tEzx;
	cube_split.Ezy(i, j, k) = tEzy;

	/*

	if (sizeof(ftype) == 2 && sizeof(ftypePML) == 4) {
		ftypePML temp = tExy + tExz;
		cube.Ex(i, j, k) = cl::sycl::detail::float2Half(tExy + tExz);
		if (temp != 0) {
			out << temp << "  " << cube.Ex(i, j, k) << cl::sycl::endl;
		}
		temp = tEyx + tEyz;
		cube.Ey(i, j, k) = cl::sycl::detail::float2Half(tEyx + tEyz);
		//if (temp != 0) {
		//	out << temp << "  " << cube.Ey(i, j, k) << cl::sycl::endl;
		//}

		temp = tEzx + tEzy;
		cube.Ez(i, j, k) = cl::sycl::detail::float2Half(tEzx + tEzy);
		//if (temp != 0) {
		//	out << temp << "  " << cube.Ez(i, j, k) << cl::sycl::endl;
		//}
	}
	else {
		cube.Ex(i, j, k) = (tExy + tExz);
		if (cube.Ex(i, j, k) > 0) {
			//out << tExy + tExz << "  " << cube.Ex(i, j, k) << cl::sycl::endl;
		}

		cube.Ey(i, j, k) = (tEyx + tEyz);
		if (cube.Ey(i, j, k) > 0) {
			//out << tEyx + tEyz << "  " << cube.Ey(i, j, k) << cl::sycl::endl;
		}
		cube.Ez(i, j, k) = (tEzx + tEzy);
		if (cube.Ez(i, j, k) > 0) {
			//out << tEzx + tEzy << "  " << cube.Ez(i, j, k) << cl::sycl::endl;
		}
	}

	*/
}

template <class ftype, class ftypePML>
inline void Update_magnetic_PML(SplitFields<ftypePML> cube_split, Coefficient<ftypePML> Coeff,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k)
{
	ftypePML tBxy, tBxz, tByx, tByz, tBzx, tBzy;

	//if (IsPrint)
	//	std::cout << "Update_magnetic_PML: " << " i=" << i << " j=" << j << " k=" << k << std::endl;

	tBxy = cube_split.Bxy(i, j, k) * Coeff.Bxy1(i, j, k) - ((cube_split.Ezx(i, j, k) + cube_split.Ezy(i, j, k)) - (cube_split.Ezx(i, j - 1, k) + cube_split.Ezy(i, j - 1, k))) * Coeff.Bxy2(i, j, k) * (_1dy);

	tBxz = cube_split.Bxz(i, j, k) * Coeff.Bxz1(i, j, k) + ((cube_split.Eyx(i, j, k) + cube_split.Eyz(i, j, k)) - (cube_split.Eyx(i, j, k - 1) + cube_split.Eyz(i, j, k - 1))) * Coeff.Bxz2(i, j, k) * (_1dz);

	tByx = cube_split.Byx(i, j, k) * Coeff.Byx1(i, j, k) + ((cube_split.Ezx(i, j, k) + cube_split.Ezy(i, j, k)) - (cube_split.Ezx(i - 1, j, k) + cube_split.Ezy(i - 1, j, k))) * Coeff.Byx2(i, j, k) * (_1dx);

	tByz = cube_split.Byz(i, j, k) * Coeff.Byz1(i, j, k) - ((cube_split.Exy(i, j, k) + cube_split.Exz(i, j, k)) - (cube_split.Exy(i, j, k - 1) + cube_split.Exz(i, j, k - 1))) * Coeff.Byz2(i, j, k) * (_1dz);

	tBzx = cube_split.Bzx(i, j, k) * Coeff.Bzx1(i, j, k) - ((cube_split.Eyx(i, j, k) + cube_split.Eyz(i, j, k)) - (cube_split.Eyx(i - 1, j, k) + cube_split.Eyz(i - 1, j, k))) * Coeff.Bzx2(i, j, k) * (_1dx);

	tBzy = cube_split.Bzy(i, j, k) * Coeff.Bzy1(i, j, k) + ((cube_split.Exy(i, j, k) + cube_split.Exz(i, j, k)) - (cube_split.Exy(i, j - 1, k) + cube_split.Exz(i, j - 1, k))) * Coeff.Bzy2(i, j, k) * (_1dy);
	if (std::abs(tBxy) <= 1e-150) {
		tBxy = 0.0;
	}
	if (std::abs(tBxz) <= 1e-150) {
		tBxz = 0.0;
	}
	if (std::abs(tByx) <= 1e-150) {
		tByx = 0.0;
	}
	if (std::abs(tByz) <= 1e-150) {
		tByz = 0.0;
	}
	if (std::abs(tBzx) <= 1e-150) {
		tBzx = 0.0;
	}
	if (std::abs(tBzy) <= 1e-150) {
		tBzy = 0.0;
	}

	cube_split.Bxy(i, j, k) = tBxy;
	cube_split.Bxz(i, j, k) = tBxz;
	cube_split.Byx(i, j, k) = tByx;
	cube_split.Byz(i, j, k) = tByz;
	cube_split.Bzx(i, j, k) = tBzx;
	cube_split.Bzy(i, j, k) = tBzy;
}

// PML BOUND
template <class ftype, class ftypePML>
inline void Update_electric_PML_BOUND_half_version(Fields<ftype> cube,
	SplitFields<ftypePML> cube_split, Coefficient<ftypePML> Coeff,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k, sycl::stream out)
{
	ftypePML tExy, tExz, tEyx, tEyz, tEzx, tEzy;

	//c1*(A + b*c2)
	tExy = cube_split.Exy(i, j, k) * Coeff.Exy1(i, j, k) + (float(cube.Bz(i, j + 1, k)) - float(cube.Bz(i, j, k))) * Coeff.Exy2(i, j, k) * (_1dy);

	tExz = cube_split.Exz(i, j, k) * Coeff.Exz1(i, j, k) - (float(cube.By(i, j, k + 1)) - float(cube.By(i, j, k))) * Coeff.Exz2(i, j, k) * (_1dz);

	tEyx = cube_split.Eyx(i, j, k) * Coeff.Eyx1(i, j, k) - (float(cube.Bz(i + 1, j, k)) - float(cube.Bz(i, j, k))) * Coeff.Eyx2(i, j, k) * (_1dx);

	tEyz = cube_split.Eyz(i, j, k) * Coeff.Eyz1(i, j, k) + (float(cube.Bx(i, j, k + 1)) - float(cube.Bx(i, j, k))) * Coeff.Eyz2(i, j, k) * (_1dz);

	tEzx = cube_split.Ezx(i, j, k) * Coeff.Ezx1(i, j, k) + (float(cube.By(i + 1, j, k)) - float(cube.By(i, j, k))) * Coeff.Ezx2(i, j, k) * (_1dx);

	tEzy = cube_split.Ezy(i, j, k) * Coeff.Ezy1(i, j, k) - (float(cube.Bx(i, j + 1, k)) - float(cube.Bx(i, j, k))) * Coeff.Ezy2(i, j, k) * (_1dy);


	//if (std::abs(tExy) > 0.0) {
	//	out << tExy;
	//}
	//if (std::abs(tExz) > 0.0) {
	//	out << tExz;
	//}
	//if (std::abs(tEyx) > 0.0) {
	//	out << tEyx;
	//}
	//if (std::abs(tEyz) > 0.0) {
	//	out << tEyz;
	//}
	//if (std::abs(tEzx) > 0.0) {
	//	out << tEzx;
	//}
	//if (std::abs(tEzy) > 0.0) {
	//	out << tEzy;
	//}

	cube_split.Exy(i, j, k) = tExy;
	cube_split.Exz(i, j, k) = tExz;
	cube_split.Eyx(i, j, k) = tEyx;
	cube_split.Eyz(i, j, k) = tEyz;
	cube_split.Ezx(i, j, k) = tEzx;
	cube_split.Ezy(i, j, k) = tEzy;
}

template <class ftype, class ftypePML>
inline void Update_magnetic_PML_BOUND_half_version(Fields<ftype> cube,
	SplitFields<ftypePML> cube_split, Coefficient<ftypePML> Coeff,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k, sycl::stream out)
{
	ftypePML tBxy, tBxz, tByx, tByz, tBzx, tBzy;

	tBxy = cube_split.Bxy(i, j, k) * Coeff.Bxy1(i, j, k) - (float(cube.Ez(i, j, k)) - float(cube.Ez(i, j - 1, k))) * Coeff.Bxy2(i, j, k) * (_1dy);

	tBxz = cube_split.Bxz(i, j, k) * Coeff.Bxz1(i, j, k) + (float(cube.Ey(i, j, k)) - float(cube.Ey(i, j, k - 1))) * Coeff.Bxz2(i, j, k) * (_1dz);

	tByx = cube_split.Byx(i, j, k) * Coeff.Byx1(i, j, k) + (float(cube.Ez(i, j, k)) - float(cube.Ez(i - 1, j, k))) * Coeff.Byx2(i, j, k) * (_1dx);

	tByz = cube_split.Byz(i, j, k) * Coeff.Byz1(i, j, k) - (float(cube.Ex(i, j, k)) - float(cube.Ex(i, j, k - 1))) * Coeff.Byz2(i, j, k) * (_1dz);

	tBzx = cube_split.Bzx(i, j, k) * Coeff.Bzx1(i, j, k) - (float(cube.Ey(i, j, k)) - float(cube.Ey(i - 1, j, k))) * Coeff.Bzx2(i, j, k) * (_1dx);

	tBzy = cube_split.Bzy(i, j, k) * Coeff.Bzy1(i, j, k) + (float(cube.Ex(i, j, k)) - float(cube.Ex(i, j - 1, k))) * Coeff.Bzy2(i, j, k) * (_1dy);

	//if (std::abs(tBxy) <= 1e-150) {
	//	tBxy = 0.0;
	//}
	//if (std::abs(tBxz) <= 1e-150) {
	//	tBxz = 0.0;
	//}
	//if (std::abs(tByx) <= 1e-150) {
	//	tByx = 0.0;
	//}
	//if (std::abs(tByz) <= 1e-150) {
	//	tByz = 0.0;
	//}
	//if (std::abs(tBzx) <= 1e-150) {
	//	tBzx = 0.0;
	//}
	//if (std::abs(tBzy) <= 1e-150) {
	//	tBzy = 0.0;
	//}

	cube_split.Bxy(i, j, k) = tBxy;
	cube_split.Bxz(i, j, k) = tBxz;
	cube_split.Byx(i, j, k) = tByx;
	cube_split.Byz(i, j, k) = tByz;
	cube_split.Bzx(i, j, k) = tBzx;
	cube_split.Bzy(i, j, k) = tBzy;
}

template <class ftype, class ftypePML>
inline void Update_magnetic_PML_BOUND_float_version(Fields<ftype> cube,
	SplitFields<ftypePML> cube_split, Coefficient<ftypePML> Coeff,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k, sycl::stream out)
{
	ftypePML tBxy, tBxz, tByx, tByz, tBzx, tBzy;

	tBxy = cube_split.Bxy(i, j, k) * Coeff.Bxy1(i, j, k) - ((ftypePML)cube.Ez(i, j, k) - (ftypePML)cube.Ez(i, j - 1, k)) * Coeff.Bxy2(i, j, k) * (_1dy);

	tBxz = cube_split.Bxz(i, j, k) * Coeff.Bxz1(i, j, k) + ((ftypePML)cube.Ey(i, j, k) - (ftypePML)cube.Ey(i, j, k - 1)) * Coeff.Bxz2(i, j, k) * (_1dz);

	tByx = cube_split.Byx(i, j, k) * Coeff.Byx1(i, j, k) + ((ftypePML)cube.Ez(i, j, k) - (ftypePML)cube.Ez(i - 1, j, k)) * Coeff.Byx2(i, j, k) * (_1dx);

	tByz = cube_split.Byz(i, j, k) * Coeff.Byz1(i, j, k) - ((ftypePML)cube.Ex(i, j, k) - (ftypePML)cube.Ex(i, j, k - 1)) * Coeff.Byz2(i, j, k) * (_1dz);

	tBzx = cube_split.Bzx(i, j, k) * Coeff.Bzx1(i, j, k) - ((ftypePML)cube.Ey(i, j, k) - (ftypePML)cube.Ey(i - 1, j, k)) * Coeff.Bzx2(i, j, k) * (_1dx);

	tBzy = cube_split.Bzy(i, j, k) * Coeff.Bzy1(i, j, k) + ((ftypePML)cube.Ex(i, j, k) - (ftypePML)cube.Ex(i, j - 1, k)) * Coeff.Bzy2(i, j, k) * (_1dy);

	//if (std::abs(tBxy) > 0.0) {
	//	out << tBxy;
	//}
	//if (std::abs(tBxz) > 0.0) {
	//	out << tBxz;
	//}
	//if (std::abs(tByx) > 0.0) {
	//	out << tByx;
	//}
	//if (std::abs(tByz) > 0.0) {
	//	out << tByz;
	//}
	//if (std::abs(tBzx) > 0.0) {
	//	out << tBzx;
	//}
	//if (std::abs(tBzy) > 0.0) {
	//	out << tBzy;
	//}

	cube_split.Bxy(i, j, k) = tBxy;
	cube_split.Bxz(i, j, k) = tBxz;
	cube_split.Byx(i, j, k) = tByx;
	cube_split.Byz(i, j, k) = tByz;
	cube_split.Bzx(i, j, k) = tBzx;
	cube_split.Bzy(i, j, k) = tBzy;
}


template <class ftype, class ftypePML>
inline void Update_electric_PML_BOUND_float_version(Fields<ftype> cube,
	SplitFields<ftypePML> cube_split, Coefficient<ftypePML> Coeff,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k, sycl::stream out)
{
	ftypePML tExy, tExz, tEyx, tEyz, tEzx, tEzy;

	//c1*(A + b*c2)
	tExy = cube_split.Exy(i, j, k) * Coeff.Exy1(i, j, k) + ((ftypePML)cube.Bz(i, j + 1, k) - (ftypePML)cube.Bz(i, j, k)) * Coeff.Exy2(i, j, k) * (_1dy);

	tExz = cube_split.Exz(i, j, k) * Coeff.Exz1(i, j, k) - ((ftypePML)cube.By(i, j, k + 1) - (ftypePML)cube.By(i, j, k)) * Coeff.Exz2(i, j, k) * (_1dz);

	tEyx = cube_split.Eyx(i, j, k) * Coeff.Eyx1(i, j, k) - ((ftypePML)cube.Bz(i + 1, j, k) - (ftypePML)cube.Bz(i, j, k)) * Coeff.Eyx2(i, j, k) * (_1dx);

	tEyz = cube_split.Eyz(i, j, k) * Coeff.Eyz1(i, j, k) + ((ftypePML)cube.Bx(i, j, k + 1) - (ftypePML)cube.Bx(i, j, k)) * Coeff.Eyz2(i, j, k) * (_1dz);

	tEzx = cube_split.Ezx(i, j, k) * Coeff.Ezx1(i, j, k) + ((ftypePML)cube.By(i + 1, j, k) - (ftypePML)cube.By(i, j, k)) * Coeff.Ezx2(i, j, k) * (_1dx);

	tEzy = cube_split.Ezy(i, j, k) * Coeff.Ezy1(i, j, k) - ((ftypePML)cube.Bx(i, j + 1, k) - (ftypePML)cube.Bx(i, j, k)) * Coeff.Ezy2(i, j, k) * (_1dy);

	//if (std::abs(tExy) > 0.0) {
	//	out << tExy;
	//}
	//if (std::abs(tExz) > 0.0) {
	//	out << tExz;
	//}
	//if (std::abs(tEyx) > 0.0) {
	//	out << tEyx;
	//}
	//if (std::abs(tEyz) > 0.0) {
	//	out << tEyz;
	//}
	//if (std::abs(tEzx) > 0.0) {
	//	out << tEzx;
	//}
	//if (std::abs(tEzy) > 0.0) {
	//	out << tEzy;
	//}

	cube_split.Exy(i, j, k) = tExy;
	cube_split.Exz(i, j, k) = tExz;
	cube_split.Eyx(i, j, k) = tEyx;
	cube_split.Eyz(i, j, k) = tEyz;
	cube_split.Ezx(i, j, k) = tEzx;
	cube_split.Ezy(i, j, k) = tEzy;
}

//REVERSE MATCHING

template <class ftype, class ftypePML>
void CubeSplit2Cube_electric_float_version(Fields<ftype> cube,
	SplitFields<ftypePML> cube_split, int i, int j, int k) {

	cube.Ex(i, j, k) = (ftype)(cube_split.Exy(i, j, k) + cube_split.Exz(i, j, k));
	cube.Ey(i, j, k) = (ftype)(cube_split.Eyx(i, j, k) + cube_split.Eyz(i, j, k));
	cube.Ez(i, j, k) = (ftype)(cube_split.Ezx(i, j, k) + cube_split.Ezy(i, j, k));

}

template <class ftype, class ftypePML>
void CubeSplit2Cube_electric_half_version(Fields<ftype> cube,
	SplitFields<ftypePML> cube_split, int i, int j, int k) {

	cube.Ex(i, j, k) = ftype((cube_split.Exy(i, j, k) + cube_split.Exz(i, j, k)));
	cube.Ey(i, j, k) = ftype((cube_split.Eyx(i, j, k) + cube_split.Eyz(i, j, k)));
	cube.Ez(i, j, k) = ftype((cube_split.Ezx(i, j, k) + cube_split.Ezy(i, j, k)));

}

template <class ftype, class ftypePML>
inline void CubeSplit2Cube_magnetic_float_version(Fields<ftype> cube,
	SplitFields<ftypePML> cube_split, int i, int j, int k) {

	cube.Bx(i, j, k) = (ftype)(cube_split.Bxy(i, j, k) + cube_split.Bxz(i, j, k));
	cube.By(i, j, k) = (ftype)(cube_split.Byx(i, j, k) + cube_split.Byz(i, j, k));
	cube.Bz(i, j, k) = (ftype)(cube_split.Bzx(i, j, k) + cube_split.Bzy(i, j, k));
}

template <class ftype, class ftypePML>
inline void CubeSplit2Cube_magnetic_half_version(Fields<ftype> cube,
	SplitFields<ftypePML> cube_split, int i, int j, int k) {

	cube.Bx(i, j, k) = ftype((cube_split.Bxy(i, j, k) + cube_split.Bxz(i, j, k)));
	cube.By(i, j, k) = ftype((cube_split.Byx(i, j, k) + cube_split.Byz(i, j, k)));
	cube.Bz(i, j, k) = ftype((cube_split.Bzx(i, j, k) + cube_split.Bzy(i, j, k)));
}





//template <class ftype, class ftypePML>
//inline void Update_electric_field_PML_sycl_half_version_bound(Fields<ftype> cube,
//	SplitFields<ftypePML> cube_split, Coefficient<ftypePML> Coeff,
//	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k)
//{
//	ftypePML tExy, tExz, tEyx, tEyz, tEzx, tEzy;
//
//	//c1*(A + b*c2)
//	tExy = cube_split.Exy(i, j, k) * Coeff.Exy1(i, j, k) + (float(cube.Bz(i, j + 1, k)) - float(cube.Bz(i, j, k))) * Coeff.Exy2(i, j, k) * (_1dy);
//
//	tExz = cube_split.Exz(i, j, k) * Coeff.Exz1(i, j, k) - (float(cube.By(i, j, k + 1)) - float(cube.By(i, j, k))) * Coeff.Exz2(i, j, k) * (_1dz);
//
//	tEyx = cube_split.Eyx(i, j, k) * Coeff.Eyx1(i, j, k) - (float(cube.Bz(i + 1, j, k)) - float(cube.Bz(i, j, k))) * Coeff.Eyx2(i, j, k) * (_1dx);
//
//	tEyz = cube_split.Eyz(i, j, k) * Coeff.Eyz1(i, j, k) + (float(cube.Bx(i, j, k + 1)) - float(cube.Bx(i, j, k))) * Coeff.Eyz2(i, j, k) * (_1dz);
//
//	tEzx = cube_split.Ezx(i, j, k) * Coeff.Ezx1(i, j, k) + (float(cube.By(i + 1, j, k)) - float(cube.By(i, j, k))) * Coeff.Ezx2(i, j, k) * (_1dx);
//
//	tEzy = cube_split.Ezy(i, j, k) * Coeff.Ezy1(i, j, k) - (float(cube.Bx(i, j + 1, k)) - float(cube.Bx(i, j, k))) * Coeff.Ezy2(i, j, k) * (_1dy);
//
//
//	if (std::abs(tExy) <= 1e-150) {
//		tExy = 0.0;
//	}
//	if (std::abs(tExz) <= 1e-150) {
//		tExz = 0.0;
//	}
//	if (std::abs(tEyx) <= 1e-150) {
//		tEyx = 0.0;
//	}
//	if (std::abs(tEyz) <= 1e-150) {
//		tEyz = 0.0;
//	}
//	if (std::abs(tEzx) <= 1e-150) {
//		tEzx = 0.0;
//	}
//	if (std::abs(tEzy) <= 1e-150) {
//		tEzy = 0.0;
//	}
//
//	cube_split.Exy(i, j, k) = tExy;
//	cube_split.Exz(i, j, k) = tExz;
//	cube_split.Eyx(i, j, k) = tEyx;
//	cube_split.Eyz(i, j, k) = tEyz;
//	cube_split.Ezx(i, j, k) = tEzx;
//	cube_split.Ezy(i, j, k) = tEzy;
//}
//





//template <class ftype, class ftypePML>
//void Loop_Update_PMLBound_electric(int st1, int fn1, int st2, int fn2, int st3, int fn3, Fields<ftype> cube,
//	SplitFields<ftypePML> cube_split, Coefficient<ftypePML> Coeff,
//	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz) {
//	for (int i = st1; i < fn1; i += 1)
//		for (int j = st2; j < fn2; j += 1)
//			for (int k = st3; k < fn3; k += 1) {
//				Update_electric_field_PML_sycl_half_version_bound<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
//			}
//}

template <class ftype, class ftypePML>
void Update_PMLBound_electric(Fields<ftype> cube,
	SplitFields<ftypePML> cube_split, Coefficient<ftypePML> Coeff, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z, ftypePML _1dx, ftypePML _1dy, ftypePML _1dz) {
	/*
	//XY
	for (int i = delta_x + 1; i < Nx + delta_x + 1; i += 1)
		for (int j = delta_y + 1; j < Ny + delta_y + 1; j += 1)
			for (int k = delta_z; k < delta_z + 1; k += 1) {
				Update_electric_field_PML_sycl_half_version_bound<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);

			}
	//XZ
	for (int i = delta_x + 1; i < Nx + delta_x + 1; i += 1)
		for (int j = delta_y; j < delta_y + 1; j += 1)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k += 1) {
				Update_electric_field_PML_sycl_half_version_bound<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
			}
	//YZ
	for (int i = delta_x; i < delta_x + 1; i += 1)
		for (int j = delta_y + 1; j < Ny + delta_y + 1; j += 1)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k += 1) {
				Update_electric_field_PML_sycl_half_version_bound<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
			}
			*/

	if (sizeof(ftype) == 2) {
		// YZ
		for (int i = delta_x; i != 0; i = 0)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++) {
					Update_electric_field_PML_sycl_bound_2_half_version<ftype, ftypePML>
						(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
				}

		// XY
		for (int k = delta_z; k != 0; k = 0)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int i = delta_x + 1; i < Nx + delta_x + 1; i++) {
					Update_electric_field_PML_sycl_bound_2_half_version<ftype, ftypePML>
						(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
				}
		// XZ
		for (int j = delta_y; j != 0; j = 0)
			for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++) {
					Update_electric_field_PML_sycl_bound_2_half_version<ftype, ftypePML>
						(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
				}
	}
	else {
		// YZ
		for (int i = delta_x; i != 0; i = 0)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++) {
					Update_electric_field_PML_sycl_bound_2_float_version<ftype, ftypePML>
						(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
				}

		// XY
		for (int k = delta_z; k != 0; k = 0)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int i = delta_x + 1; i < Nx + delta_x + 1; i++) {
					Update_electric_field_PML_sycl_bound_2_float_version<ftype, ftypePML>
						(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
				}
		// XZ
		for (int j = delta_y; j != 0; j = 0)
			for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++) {
					Update_electric_field_PML_sycl_bound_2_float_version<ftype, ftypePML>
						(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
				}
	}

	
}

/*
template <class ftype, class ftypePML>
inline void Update_magnetic_field_PML_sycl_half_version(Fields<ftype> cube,
	SplitFields<ftypePML> cube_split, Coefficient<ftypePML> Coeff,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k)
{
	ftypePML tBxy, tBxz, tByx, tByz, tBzx, tBzy;
	if (sizeof(ftype) == 2 && sizeof(ftypePML) == 4) {
		tBxy = cube_split.Bxy(i, j, k) * Coeff.Bxy1(i, j, k) - (float(cube.Ez(i, j, k)) - float(cube.Ez(i, j - 1, k))) * Coeff.Bxy2(i, j, k) * (_1dy);

		tBxz = cube_split.Bxz(i, j, k) * Coeff.Bxz1(i, j, k) + (float(cube.Ey(i, j, k)) - float(cube.Ey(i, j, k - 1))) * Coeff.Bxz2(i, j, k) * (_1dz);
		
		tByx = cube_split.Byx(i, j, k) * Coeff.Byx1(i, j, k) + (float(cube.Ez(i, j, k)) - float(cube.Ez(i - 1, j, k))) * Coeff.Byx2(i, j, k) * (_1dx);

		tByz = cube_split.Byz(i, j, k) * Coeff.Byz1(i, j, k) - (float(cube.Ex(i, j, k)) - float(cube.Ex(i, j, k - 1))) * Coeff.Byz2(i, j, k) * (_1dz);

		tBzx = cube_split.Bzx(i, j, k) * Coeff.Bzx1(i, j, k) - (float(cube.Ey(i, j, k)) - float(cube.Ey(i - 1, j, k))) * Coeff.Bzx2(i, j, k) * (_1dx);

		tBzy = cube_split.Bzy(i, j, k) * Coeff.Bzy1(i, j, k) + (float(cube.Ex(i, j, k)) - float(cube.Ex(i, j - 1, k))) * Coeff.Bzy2(i, j, k) * (_1dy);
	} else {
		tBxy = cube_split.Bxy(i, j, k) * Coeff.Bxy1(i, j, k) - (cube.Ez(i, j, k) - cube.Ez(i, j - 1, k)) * Coeff.Bxy2(i, j, k) * (_1dy);

		tBxz = cube_split.Bxz(i, j, k) * Coeff.Bxz1(i, j, k) + (cube.Ey(i, j, k) - cube.Ey(i, j, k - 1)) * Coeff.Bxz2(i, j, k) * (_1dz);

		tByx = cube_split.Byx(i, j, k) * Coeff.Byx1(i, j, k) + (cube.Ez(i, j, k) - cube.Ez(i - 1, j, k)) * Coeff.Byx2(i, j, k) * (_1dx);

		tByz = cube_split.Byz(i, j, k) * Coeff.Byz1(i, j, k) - (cube.Ex(i, j, k) - cube.Ex(i, j, k - 1)) * Coeff.Byz2(i, j, k) * (_1dz);

		tBzx = cube_split.Bzx(i, j, k) * Coeff.Bzx1(i, j, k) - (cube.Ey(i, j, k) - cube.Ey(i - 1, j, k)) * Coeff.Bzx2(i, j, k) * (_1dx);

		tBzy = cube_split.Bzy(i, j, k) * Coeff.Bzy1(i, j, k) + (cube.Ex(i, j, k) - cube.Ex(i, j - 1, k)) * Coeff.Bzy2(i, j, k) * (_1dy);
	}
	if (std::abs(tBxy) <= 1e-150) {
		tBxy = 0.0;
	}
	if (std::abs(tBxz) <= 1e-150) {
		tBxz = 0.0;
	}
	if (std::abs(tByx) <= 1e-150) {
		tByx = 0.0;
	}
	if (std::abs(tByz) <= 1e-150) {
		tByz = 0.0;
	}
	if (std::abs(tBzx) <= 1e-150) {
		tBzx = 0.0;
	}
	if (std::abs(tBzy) <= 1e-150) {
		tBzy = 0.0;
	}

	cube_split.Bxy(i, j, k) = tBxy;
	cube_split.Bxz(i, j, k) = tBxz;
	cube_split.Byx(i, j, k) = tByx;
	cube_split.Byz(i, j, k) = tByz;
	cube_split.Bzx(i, j, k) = tBzx;
	cube_split.Bzy(i, j, k) = tBzy;

	if (sizeof(ftype) == 2 && sizeof(ftypePML) == 4) {
		cube.Bx(i, j, k) = cl::sycl::detail::float2Half(tBxy + tBxz);
		cube.By(i, j, k) = cl::sycl::detail::float2Half(tByx + tByz);
		cube.Bz(i, j, k) = cl::sycl::detail::float2Half(tBzx + tBzy);
	}
	else {
		cube.Bx(i, j, k) = (tBxy + tBxz);
		cube.By(i, j, k) = (tByx + tByz);
		cube.Bz(i, j, k) = (tBzx + tBzy);
	}
}
*/


//template <class ftype, class ftypePML>
//inline void Update_magnetic_field_PML_sycl_half_version_bound(Fields<ftype> cube,
//	SplitFields<ftypePML> cube_split, Coefficient<ftypePML> Coeff,
//	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k)
//{
//	ftypePML tBxy, tBxz, tByx, tByz, tBzx, tBzy;
//
//	tBxy = cube_split.Bxy(i, j, k) * Coeff.Bxy1(i, j, k) - ((ftypePML)float(cube.Ez(i, j, k)) - (ftypePML)float(cube.Ez(i, j - 1, k))) * Coeff.Bxy2(i, j, k) * (_1dy);
//
//	tBxz = cube_split.Bxz(i, j, k) * Coeff.Bxz1(i, j, k) + ((ftypePML)float(cube.Ey(i, j, k)) - (ftypePML)float(cube.Ey(i, j, k - 1))) * Coeff.Bxz2(i, j, k) * (_1dz);
//	
//	tByx = cube_split.Byx(i, j, k) * Coeff.Byx1(i, j, k) + ((ftypePML)float(cube.Ez(i, j, k)) - (ftypePML)float(cube.Ez(i - 1, j, k))) * Coeff.Byx2(i, j, k) * (_1dx);
//
//	tByz = cube_split.Byz(i, j, k) * Coeff.Byz1(i, j, k) - ((ftypePML)float(cube.Ex(i, j, k)) - (ftypePML)float(cube.Ex(i, j, k - 1))) * Coeff.Byz2(i, j, k) * (_1dz);
//
//	tBzx = cube_split.Bzx(i, j, k) * Coeff.Bzx1(i, j, k) - ((ftypePML)float(cube.Ey(i, j, k)) - (ftypePML)float(cube.Ey(i - 1, j, k))) * Coeff.Bzx2(i, j, k) * (_1dx);
//
//	tBzy = cube_split.Bzy(i, j, k) * Coeff.Bzy1(i, j, k) + ((ftypePML)float(cube.Ex(i, j, k)) - (ftypePML)float(cube.Ex(i, j - 1, k))) * Coeff.Bzy2(i, j, k) * (_1dy);
//
//	if (std::abs(tBxy) <= 1e-150) {
//		tBxy = 0.0;
//	}
//	if (std::abs(tBxz) <= 1e-150) {
//		tBxz = 0.0;
//	}
//	if (std::abs(tByx) <= 1e-150) {
//		tByx = 0.0;
//	}
//	if (std::abs(tByz) <= 1e-150) {
//		tByz = 0.0;
//	}
//	if (std::abs(tBzx) <= 1e-150) {
//		tBzx = 0.0;
//	}
//	if (std::abs(tBzy) <= 1e-150) {
//		tBzy = 0.0;
//	}
//
//	cube_split.Bxy(i, j, k) = tBxy;
//	cube_split.Bxz(i, j, k) = tBxz;
//	cube_split.Byx(i, j, k) = tByx;
//	cube_split.Byz(i, j, k) = tByz;
//	cube_split.Bzx(i, j, k) = tBzx;
//	cube_split.Bzy(i, j, k) = tBzy;
//}



template <class ftype, class ftypePML>
void Update_PMLBound_magnetic(Fields<ftype> cube,
	SplitFields<ftypePML> cube_split, Coefficient<ftypePML> Coeff, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z, ftypePML _1dx, ftypePML _1dy, ftypePML _1dz) {
	/*
	//edge from n,0,0 by YZ
	for (int i = Nx + delta_x + 1; i < Nx + delta_x + 2; i += 1)
		for (int j = delta_y + 1; j < Ny + delta_y + 1; j += 1)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k += 1) {
				Update_magnetic_field_PML_sycl_half_version_bound<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
			}

	//edge from 0,m,0 by XZ
	for (int i = delta_x + 1; i < Nx + delta_x + 1; i += 1)
		for (int j = Ny + delta_y + 1; j < Ny + delta_y + 2; j += 1)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k += 1) {
				Update_magnetic_field_PML_sycl_half_version_bound<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
			}

	//edge from 0,0,k by XY
	for (int i = delta_x + 1; i < Nx + delta_x + 1; i += 1)
		for (int j = delta_y + 1; j < Ny + delta_y + 1; j += 1)
			for (int k = Nz + delta_z + 1; k < Nz + delta_z + 2; k += 1) {
				Update_magnetic_field_PML_sycl_half_version_bound<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
			}
	*/

	if (sizeof(ftype) == 2) {
		// YZ
		for (int i = Nx + delta_x + 1; i != Nx + 1; i = Nx + 1)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++) {
					Update_magnetic_field_PML_sycl_bound_2_half_version<ftype, ftypePML>
						(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
				}
		// XY
		for (int k = Nz + delta_z + 1; k != Nz + 1; k = Nz + 1)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int i = delta_x + 1; i < Nx + delta_x + 1; i++) {
					Update_magnetic_field_PML_sycl_bound_2_half_version<ftype, ftypePML>
						(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
				}
		// XZ
		for (int j = Ny + delta_y + 1; j != Ny + 1; j = Ny + 1)
			for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++) {
					Update_magnetic_field_PML_sycl_bound_2_half_version<ftype, ftypePML>
						(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
				}
	}
	else {
		// YZ
		for (int i = Nx + delta_x + 1; i != Nx + 1; i = Nx + 1)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++) {
					Update_magnetic_field_PML_sycl_bound_2_float_version<ftype, ftypePML>
						(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
				}
		// XY
		for (int k = Nz + delta_z + 1; k != Nz + 1; k = Nz + 1)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int i = delta_x + 1; i < Nx + delta_x + 1; i++) {
					Update_magnetic_field_PML_sycl_bound_2_float_version<ftype, ftypePML>
						(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
				}
		// XZ
		for (int j = Ny + delta_y + 1; j != Ny + 1; j = Ny + 1)
			for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++) {
					Update_magnetic_field_PML_sycl_bound_2_float_version<ftype, ftypePML>
						(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
				}
	}

	
}


template <class ftype, class ftypePML>
void Run_CubeSplit2Cube_electric(Fields<ftype> cube,
	SplitFields<ftypePML> cube_split, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z) {
	/*
	int x = Nx + 2 * delta_x, y = Ny + 2 * delta_y, z = Nz + 2 * delta_z;
	//XY
	for (int i = delta_x; i < Nx + delta_x + 1; i += 1)
		for (int j = delta_y; j < Ny + delta_y + 1; j += 1)
			for (int k = delta_z; k < delta_z + 1; k += 1) {
				CubeSplit2Cube_electric<ftype, ftypePML>(cube, cube_split, i, j, k);

			}
	//XZ
	for (int i = delta_x; i < Nx + delta_x + 1; i += 1)
		for (int j = delta_y; j < delta_y + 1; j += 1)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k += 1) {
				CubeSplit2Cube_electric<ftype, ftypePML>(cube, cube_split, i, j, k);
			}
	//YZ
	for (int i = delta_x; i < delta_x + 1; i += 1)
		for (int j = delta_y + 1; j < Ny + delta_y + 1; j += 1)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k += 1) {
				CubeSplit2Cube_electric<ftype, ftypePML>(cube, cube_split, i, j, k);
			}

	//edge from n,0,0 by YZ
	for (int i = Nx + delta_x + 1; i < Nx + delta_x + 2; i += 1)
		for (int j = delta_y; j < Ny + delta_y + 2; j += 1)
			for (int k = delta_z; k < Nz + delta_z + 2; k += 1) {
				CubeSplit2Cube_electric<ftype, ftypePML>(cube, cube_split, i, j, k);

			}
	//edge from 0,m,0 by XZ
	for (int i = delta_x; i < Nx + delta_x + 1; i += 1) // укоротили на один чтобы не было пересечений
		for (int j = Ny + delta_y + 1; j < Ny + delta_y + 2; j += 1)
			for (int k = delta_z; k < Nz + delta_z + 2; k += 1) {
				CubeSplit2Cube_electric<ftype, ftypePML>(cube, cube_split, i, j, k);
			}

	//edge from 0,0,k by XY
	for (int i = delta_x; i < Nx + delta_x + 1; i += 1)// укоротили на один чтобы не было пересечений
		for (int j = delta_y; j < Ny + delta_y + 1; j += 1)// укоротили на один чтобы не было пересечений
			for (int k = Nz + delta_z + 1; k < Nz + delta_z + 2; k += 1) {
				CubeSplit2Cube_electric<ftype, ftypePML>(cube, cube_split, i, j, k);
			}
			*/

	if (sizeof(ftype) == 2) {
		// YZ
		for (int i = delta_x; i != 0; i = 0)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
					CubeSplit2Cube_electric_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);

		// XY
		for (int k = delta_z; k != 0; k = 0)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
					CubeSplit2Cube_electric_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);
		// XZ
		for (int j = delta_y; j != 0; j = 0)
			for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
					CubeSplit2Cube_electric_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);

		// YZ x=n
		for (int i = Nx + delta_x + 1; i != Nx + 1; i = Nx + 1)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
					CubeSplit2Cube_electric_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);
		// XY z=n
		for (int k = Nz + delta_z + 1; k != Nz + 1; k = Nz + 1)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
					CubeSplit2Cube_electric_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);
		// XZ y=n
		for (int j = Ny + delta_y + 1; j != Ny + 1; j = Ny + 1)
			for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
					CubeSplit2Cube_electric_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);
	}
	else {
		// YZ
		for (int i = delta_x; i != 0; i = 0)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
					CubeSplit2Cube_electric_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);

		// XY
		for (int k = delta_z; k != 0; k = 0)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
					CubeSplit2Cube_electric_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);
		// XZ
		for (int j = delta_y; j != 0; j = 0)
			for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
					CubeSplit2Cube_electric_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);

		// YZ x=n
		for (int i = Nx + delta_x + 1; i != Nx + 1; i = Nx + 1)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
					CubeSplit2Cube_electric_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);
		// XY z=n
		for (int k = Nz + delta_z + 1; k != Nz + 1; k = Nz + 1)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
					CubeSplit2Cube_electric_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);
		// XZ y=n
		for (int j = Ny + delta_y + 1; j != Ny + 1; j = Ny + 1)
			for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
					CubeSplit2Cube_electric_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);
	}
	
	
}




template <class ftype, class ftypePML>
void Run_CubeSplit2Cube_magnetic(Fields<ftype> cube,
	SplitFields<ftypePML> cube_split, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z) {
	
	/*
	int x = Nx + 2 * delta_x, y = Ny + 2 * delta_y, z = Nz + 2 * delta_z;
	//XY
	for (int i = delta_x; i < Nx + delta_x + 1; i += 1)
		for (int j = delta_y; j < Ny + delta_y + 1; j += 1)
			for (int k = delta_z; k < delta_z + 1; k += 1) {
				CubeSplit2Cube_magnetic<ftype, ftypePML>(cube, cube_split, i, j, k);

			}
	//XZ
	for (int i = delta_x; i < Nx + delta_x + 1; i += 1)
		for (int j = delta_y; j < delta_y + 1; j += 1)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k += 1) {
				CubeSplit2Cube_magnetic<ftype, ftypePML>(cube, cube_split, i, j, k);
			}
	//YZ
	for (int i = delta_x; i < delta_x + 1; i += 1)
		for (int j = delta_y + 1; j < Ny + delta_y + 1; j += 1)
			for (int k = delta_z + 1; k < Nz + delta_z + 1; k += 1) {
				CubeSplit2Cube_magnetic<ftype, ftypePML>(cube, cube_split, i, j, k);
			}
	//edge from n,0,0 by YZ
	for (int i = Nx + delta_x + 1; i < Nx + delta_x + 2; i += 1)
		for (int j = delta_y; j < Ny + delta_y + 2; j += 1)
			for (int k = delta_z; k < Nz + delta_z + 2; k += 1) {
				CubeSplit2Cube_magnetic<ftype, ftypePML>(cube, cube_split, i, j, k);

			}
	//edge from 0,m,0 by XZ
	for (int i = delta_x; i < Nx + delta_x + 1; i += 1) // укоротили на один чтобы не было пересечений
		for (int j = Ny + delta_y + 1; j < Ny + delta_y + 2; j += 1)
			for (int k = delta_z; k < Nz + delta_z + 2; k += 1) {
				CubeSplit2Cube_magnetic<ftype, ftypePML>(cube, cube_split, i, j, k);
			}

	//edge from 0,0,k by XY
	for (int i = delta_x; i < Nx + delta_x + 1; i += 1)// укоротили на один чтобы не было пересечений
		for (int j = delta_y; j < Ny + delta_y + 1; j += 1)// укоротили на один чтобы не было пересечений
			for (int k = Nz + delta_z + 1; k < Nz + delta_z + 2; k += 1) {
				CubeSplit2Cube_magnetic<ftype, ftypePML>(cube, cube_split, i, j, k);
			}
*/

	if (sizeof(ftype) == 2) {
		// YZ
		for (int i = delta_x; i != 0; i = 0)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
					CubeSplit2Cube_magnetic_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);

		// XY
		for (int k = delta_z; k != 0; k = 0)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
					CubeSplit2Cube_magnetic_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);
		// XZ
		for (int j = delta_y; j != 0; j = 0)
			for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
					CubeSplit2Cube_magnetic_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);

		// YZ x=n
		for (int i = Nx + delta_x + 1; i != Nx + 1; i = Nx + 1)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
					CubeSplit2Cube_magnetic_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);
		// XY z=n
		for (int k = Nz + delta_z + 1; k != Nz + 1; k = Nz + 1)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
					CubeSplit2Cube_magnetic_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);
		// XZ y=n
		for (int j = Ny + delta_y + 1; j != Ny + 1; j = Ny + 1)
			for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
					CubeSplit2Cube_magnetic_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);
	}
	else {
		// YZ
		for (int i = delta_x; i != 0; i = 0)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
					CubeSplit2Cube_magnetic_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);

		// XY
		for (int k = delta_z; k != 0; k = 0)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
					CubeSplit2Cube_magnetic_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);
		// XZ
		for (int j = delta_y; j != 0; j = 0)
			for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
					CubeSplit2Cube_magnetic_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);

		// YZ x=n
		for (int i = Nx + delta_x + 1; i != Nx + 1; i = Nx + 1)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
					CubeSplit2Cube_magnetic_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);
		// XY z=n
		for (int k = Nz + delta_z + 1; k != Nz + 1; k = Nz + 1)
			for (int j = delta_y + 1; j < Ny + delta_y + 1; j++)
				for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
					CubeSplit2Cube_magnetic_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);
		// XZ y=n
		for (int j = Ny + delta_y + 1; j != Ny + 1; j = Ny + 1)
			for (int i = delta_x + 1; i < Nx + delta_x + 1; i++)
				for (int k = delta_z + 1; k < Nz + delta_z + 1; k++)
					CubeSplit2Cube_magnetic_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);
	}

	
}


template <class ftype, class ftypePML>
void Update_PeriodicBound_electric(Fields<ftype> cube, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z) {
	for (int i = delta_x + 1; i < Nx + delta_x + 1; i += 1) {

		cube.Ex(i, 1, 0) = cube.Ex(i, 1, 1);
		cube.Ex(i, 1, 2) = cube.Ex(i, 1, 1);
		cube.Ex(i, 0, 1) = cube.Ex(i, 1, 1);
		cube.Ex(i, 2, 1) = cube.Ex(i, 1, 1);

		cube.Ey(i, 1, 0) = cube.Ey(i, 1, 1);
		cube.Ey(i, 1, 2) = cube.Ey(i, 1, 1);
		cube.Ey(i, 0, 1) = cube.Ey(i, 1, 1);
		cube.Ey(i, 2, 1) = cube.Ey(i, 1, 1);

		cube.Ez(i, 1, 0) = cube.Ez(i, 1, 1);
		cube.Ez(i, 1, 2) = cube.Ez(i, 1, 1);
		cube.Ez(i, 0, 1) = cube.Ez(i, 1, 1);
		cube.Ez(i, 2, 1) = cube.Ez(i, 1, 1);
	}
}


template <class ftype>
void Graph_Solution_sycl(Fields<ftype>& cube, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z,
	ftype dx, ftype dy, ftype dz, std::string file_name1, std::string file_name2)
{
	ofstream numbEy(file_name1), numbBz(file_name2);

	numbEy << Nx + 2 * delta_x << ";" << Ny + 2 * delta_y << std::endl;
	numbBz << Nx + 2 * delta_x << ";" << Ny + 2 * delta_y << std::endl;

	for (int i = 1; i < Ny + 2 * delta_y + 1; i++)
	{
		ftype y = dy * (ftype)i;

		for (int j = 1; j < Nx + 2 * delta_x + 1; j++)
		{
			numbEy << cube.Ey(j, i, (Nz / 2 + delta_z) ? (Nz / 2 + delta_z) : 1) << ";";
			numbBz << cube.Bz(j, i, (Nz / 2 + delta_z) ? (Nz / 2 + delta_z) : 1) << ";";
		}
		numbEy << std::endl;
		numbBz << std::endl;

	}
	numbEy << std::endl;
	numbEy << "Nx = " << Nx << ";  Ny = " << Ny << ";   Nz = " << Nz << std::endl;
	numbEy << "delta_x = " << delta_x << ";  delta_y = " << delta_y << ";   delta_z = " << delta_z << std::endl;

	numbBz << std::endl;
	numbBz << "Nx = " << Nx << ";  Ny = " << Ny << ";   Nz = " << Nz << std::endl;
	numbBz << "delta_x = " << delta_x << ";  delta_y = " << delta_y << ";   delta_z = " << delta_z << std::endl;

	numbEy.close();
	numbBz.close();
}

template <class ftype>
double CalculateEnergy_sycl_half_version(Fields<ftype>& cube, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z, int it)
{
	double energy = 0.0, temp;
		ofstream numbEy("energy.csv");
		
	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
			{
				temp = 0.0;
				temp += (double)cube.Ex(i, j, k) * (double)cube.Ex(i, j, k)
					+ (double)cube.Ey(i, j, k) * (double)cube.Ey(i, j, k)
					+ (double)cube.Ez(i, j, k) * (double)cube.Ez(i, j, k);
				temp += (double)cube.Bx(i, j, k) * (double)cube.Bx(i, j, k)
					+ (double)cube.By(i, j, k) * (double)cube.By(i, j, k)
					+ (double)cube.Bz(i, j, k) * (double)cube.Bz(i, j, k);
				if (i == 14 && j == 25 && k == 41) {
					std::cout << cube.Ex(i, j, k) << " " << cube.Ey(i, j, k) << " " << cube.Ez(i, j, k)
						<< " " << cube.Bx(i, j, k) << " " << cube.By(i, j, k) << " " << cube.Bz(i, j, k) << std::endl;
				}
				energy += temp;
				if (it == 477) {
					numbEy << i<<";"<<j<<";"<<k<<";"<<energy << std::endl;
					if (temp > 10)
						std::cout << i << " " << j << " " << k << " " << temp << std::endl;
				}
			}
		numbEy.close();
	return energy;
}

template <class ftype>
float CalculateEnergy_sycl_float_version(sycl::queue& q, Fields<ftype>& cube, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	float host_energy = 0.0;
	try {
		float* energy = malloc_shared<float>(1, q);
		*energy = 0.f;

		q.submit([&](sycl::handler& h) {

			h.parallel_for(range<3>(Nx, Ny, Nz),
				reduction(energy, 0.f, std::plus<float>()),
				[=](auto index, auto& energy) {

					auto kernel_cube = cube;

					int i = index[0] + delta_x + 1;
					int j = index[1] + delta_y + 1;
					int k = index[2] + delta_z + 1;

					float e, b;

					e = kernel_cube.Ex(i, j, k) * kernel_cube.Ex(i, j, k)
						+ kernel_cube.Ey(i, j, k) * kernel_cube.Ey(i, j, k)
						+ kernel_cube.Ez(i, j, k) * kernel_cube.Ez(i, j, k);
					b = kernel_cube.Bx(i, j, k) * kernel_cube.Bx(i, j, k)
						+ kernel_cube.By(i, j, k) * kernel_cube.By(i, j, k)
						+ kernel_cube.Bz(i, j, k) * kernel_cube.Bz(i, j, k);

					energy += e + b;
				});
			}).wait_and_throw();

			//for (int i = 1 + delta_x; i < Nx + delta_x + 1; i++)
			//	for (int j = 1 + delta_y; j < Ny + delta_y + 1; j++)
			//		for (int k = 1 + delta_z; k < Nz + delta_z + 1; k++)
			//		{
			//			temp = 0.0;
			//			temp += (double)cube.Ex(i, j, k) * (double)cube.Ex(i, j, k)
			//				+ (double)cube.Ey(i, j, k) * (double)cube.Ey(i, j, k)
			//				+ (double)cube.Ez(i, j, k) * (double)cube.Ez(i, j, k);
			//			temp += (double)cube.Bx(i, j, k) * (double)cube.Bx(i, j, k)
			//				+ (double)cube.By(i, j, k) * (double)cube.By(i, j, k)
			//				+ (double)cube.Bz(i, j, k) * (double)cube.Bz(i, j, k);
			//			energy += temp;
			//		}
		host_energy = *energy;
		free(energy, q);
	}
	catch (sycl::exception const& e) {
		std::cout << "ERROR\n" << e.get_cl_code() << std::endl << e.what();
		std::terminate();
	}

	return host_energy;
}

template <class ftype>
float CalculateEnergy_sycl_half_version(sycl::queue& q, Fields<ftype>& cube, int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z)
{
	float host_energy = 0.0;
	try {
		float* energy = malloc_shared<float>(1, q);
		*energy = 0.f;

		q.submit([&](sycl::handler& h) {
			sycl::stream out(4096, 4096, h);
			h.parallel_for(range<3>(Nx, Ny, Nz),
				reduction(energy, 0.f, std::plus<float>()),
				[=](auto index, auto& energy) {

					auto kernel_cube = cube;

					int i = index[0] + delta_x + 1;
					int j = index[1] + delta_y + 1;
					int k = index[2] + delta_z + 1;

					//float e, b;
					//e = kernel_cube.Ex(i, j, k) * kernel_cube.Ex(i, j, k)
					//	+ kernel_cube.Ey(i, j, k) * kernel_cube.Ey(i, j, k)
					//	+ kernel_cube.Ez(i, j, k) * kernel_cube.Ez(i, j, k);
					//b = kernel_cube.Bx(i, j, k) * kernel_cube.Bx(i, j, k)
					//	+ kernel_cube.By(i, j, k) * kernel_cube.By(i, j, k)
					//	+ kernel_cube.Bz(i, j, k) * kernel_cube.Bz(i, j, k);

					float e, b;
					e = float(kernel_cube.Ex(i, j, k)) * float(kernel_cube.Ex(i, j, k))
						+ float(kernel_cube.Ey(i, j, k)) * float(kernel_cube.Ey(i, j, k))
						+ float(kernel_cube.Ez(i, j, k)) * float(kernel_cube.Ez(i, j, k));
					b = float(kernel_cube.Bx(i, j, k)) * float(kernel_cube.Bx(i, j, k))
						+ float(kernel_cube.By(i, j, k)) * float(kernel_cube.By(i, j, k))
						+ float(kernel_cube.Bz(i, j, k)) * float(kernel_cube.Bz(i, j, k));

					energy += e + b;
					//if (e > 0. || b > 0.)
						//out << e << " " << b <<   sycl::endl;
				});
			}).wait_and_throw();

			host_energy = *energy;
			free(energy, q);
	}
	catch (sycl::exception const& e) {
		std::cout << "ERROR\n" << e.get_cl_code() << std::endl << e.what();
		std::terminate();
	}

	return host_energy;
}

template <class ftypePML, class type = float>
void Initializing_Sigma_float(Sigma<type>& Sigma,
	int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z, double n,
	ftypePML sigma_x, ftypePML sigma_y)
{
	ftypePML var_Sigma_max_x = sigma_x;
	ftypePML var_Sigma_max_y = sigma_y;
	ftypePML var_Sigma_max_z = sigma_y;

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++) {
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++) {
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++) {
				Sigma.sigmaEx(i, j, k) = var_Sigma_max_x * powf(distanceE<ftypePML>(Nx, delta_x, i), n);
				Sigma.sigmaBx(i, j, k) = var_Sigma_max_x * powf(distanceH<ftypePML>(Nx, delta_x, i), n);

				Sigma.sigmaEy(i, j, k) = var_Sigma_max_y * powf(distanceE<ftypePML>(Ny, delta_y, j), n);
				Sigma.sigmaBy(i, j, k) = var_Sigma_max_y * powf(distanceH<ftypePML>(Ny, delta_y, j), n);

				Sigma.sigmaEz(i, j, k) = var_Sigma_max_z * powf(distanceE<ftypePML>(Nz, delta_z, k), n);
				Sigma.sigmaBz(i, j, k) = var_Sigma_max_z * powf(distanceH<ftypePML>(Nz, delta_z, k), n);
			}
		}
	}
	//cout << sizeof(Sigma(10, 18, 11).sigmaH_x) << sizeof(double) << sizeof(float) << endl;
}

template <class ftypePML, class ftype = double>
void Initializing_Coeff_float(Coefficient<ftypePML>& Coeff, Sigma<ftype>& Sigma, int Nx, int Ny, int Nz,
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

					Exy1 = exp(-dt * Sigma.sigmaEy(i, j, k));
					Exz1 = exp(-dt * Sigma.sigmaEz(i, j, k));
					Eyx1 = exp(-dt * Sigma.sigmaEx(i, j, k));
					Eyz1 = exp(-dt * Sigma.sigmaEz(i, j, k));
					Ezx1 = exp(-dt * Sigma.sigmaEx(i, j, k));
					Ezy1 = exp(-dt * Sigma.sigmaEy(i, j, k));

					Bxy1 = exp(-dt * Sigma.sigmaBy(i, j, k));
					Bxz1 = exp(-dt * Sigma.sigmaBz(i, j, k));
					Byx1 = exp(-dt * Sigma.sigmaBx(i, j, k));
					Byz1 = exp(-dt * Sigma.sigmaBz(i, j, k));
					Bzx1 = exp(-dt * Sigma.sigmaBx(i, j, k));
					Bzy1 = exp(-dt * Sigma.sigmaBy(i, j, k));

					Coeff.Exy1(i, j, k) = (ftypePML)Exy1;
					Coeff.Exz1(i, j, k) = (ftypePML)Exz1;
					Coeff.Eyx1(i, j, k) = (ftypePML)Eyx1;
					Coeff.Eyz1(i, j, k) = (ftypePML)Eyz1;
					Coeff.Ezx1(i, j, k) = (ftypePML)Ezx1;
					Coeff.Ezy1(i, j, k) = (ftypePML)Ezy1;

					Coeff.Bxy1(i, j, k) = (ftypePML)Bxy1;
					Coeff.Bxz1(i, j, k) = (ftypePML)Bxz1;
					Coeff.Byx1(i, j, k) = (ftypePML)Byx1;
					Coeff.Byz1(i, j, k) = (ftypePML)Byz1;
					Coeff.Bzx1(i, j, k) = (ftypePML)Bzx1;
					Coeff.Bzy1(i, j, k) = (ftypePML)Bzy1;

					if (Sigma.sigmaEx(i, j, k) != (ftypePML)0.0) {
						Coeff.Eyx2(i, j, k) = 1.0 / Sigma.sigmaEx(i, j, k) - Eyx1 / Sigma.sigmaEx(i, j, k);
						Coeff.Ezx2(i, j, k) = 1.0 / Sigma.sigmaEx(i, j, k) - Ezx1 / Sigma.sigmaEx(i, j, k);
					}
					else {
						Coeff.Eyx2(i, j, k) = dt;
						Coeff.Ezx2(i, j, k) = dt;
					}
					if (Sigma.sigmaEy(i, j, k) != (ftypePML)0.0) {
						Coeff.Exy2(i, j, k) = 1.0 / Sigma.sigmaEy(i, j, k) - Exy1 / Sigma.sigmaEy(i, j, k);
						Coeff.Ezy2(i, j, k) = 1.0 / Sigma.sigmaEy(i, j, k) - Ezy1 / Sigma.sigmaEy(i, j, k);
					}
					else {
						Coeff.Exy2(i, j, k) = dt;
						Coeff.Ezy2(i, j, k) = dt;
					}
					if (Sigma.sigmaEz(i, j, k) != (ftypePML)0.0)
					{
						Coeff.Exz2(i, j, k) = 1.0 / Sigma.sigmaEz(i, j, k) - Exz1 / Sigma.sigmaEz(i, j, k);
						Coeff.Eyz2(i, j, k) = 1.0 / Sigma.sigmaEz(i, j, k) - Eyz1 / Sigma.sigmaEz(i, j, k);
					}
					else {
						Coeff.Exz2(i, j, k) = dt;
						Coeff.Eyz2(i, j, k) = dt;
					}
					if (Sigma.sigmaBx(i, j, k) != (ftypePML)0.0)
					{
						Coeff.Byx2(i, j, k) = 1.0 / Sigma.sigmaBx(i, j, k) - Byx1 / Sigma.sigmaBx(i, j, k);
						Coeff.Bzx2(i, j, k) = 1.0 / Sigma.sigmaBx(i, j, k) - Bzx1 / Sigma.sigmaBx(i, j, k);
					}
					else {
						Coeff.Byx2(i, j, k) = dt;
						Coeff.Bzx2(i, j, k) = dt;
					}
					if (Sigma.sigmaBy(i, j, k) != (ftypePML)0.0)
					{
						Coeff.Bxy2(i, j, k) = 1.0 / Sigma.sigmaBy(i, j, k) - Bxy1 / Sigma.sigmaBy(i, j, k);
						Coeff.Bzy2(i, j, k) = 1.0 / Sigma.sigmaBy(i, j, k) - Bzy1 / Sigma.sigmaBy(i, j, k);
					}
					else {
						Coeff.Bxy2(i, j, k) = dt;
						Coeff.Bzy2(i, j, k) = dt;
					}
					if (Sigma.sigmaBz(i, j, k) != (ftypePML)0.0)
					{
						Coeff.Bxz2(i, j, k) = 1.0 / Sigma.sigmaBz(i, j, k) - Bxz1 / Sigma.sigmaBz(i, j, k);
						Coeff.Byz2(i, j, k) = 1.0 / Sigma.sigmaBz(i, j, k) - Byz1 / Sigma.sigmaBz(i, j, k);
					}
					else {
						Coeff.Bxz2(i, j, k) = dt;
						Coeff.Byz2(i, j, k) = dt;
					}
				}
			}

}