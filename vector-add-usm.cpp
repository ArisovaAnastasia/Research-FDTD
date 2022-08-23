//==============================================================
// Vector Add is the equivalent of a Hello, World! sample for data parallel
// programs. Building and running the sample verifies that your development
// environment is setup correctly and demonstrates the use of the core features
// of DPC++. This sample runs on both CPU and GPU (or FPGA). When run, it
// computes on both the CPU and offload device, then compares results. If the
// code executes on both CPU and offload device, the device name and a success
// message are displayed. And, your development environment is setup correctly!
//
// For comprehensive instructions regarding DPC++ Programming, go to
// https://software.intel.com/en-us/oneapi-programming-guide and search based on
// relevant terms noted in the comments.
//
// DPC++ material used in the code sample:
// •	A one dimensional array of data shared between CPU and offload device.
// •	A device queue and kernel.
//==============================================================
// Copyright © Intel Corporation
//
// SPDX-License-Identifier: MIT
// =============================================================
#define _USE_MATH_DEFINES
#include <CL/sycl.hpp>
#include <array>
#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <chrono>
#include "type data.h"
#include "function-for-sykl.h"

#if FPGA || FPGA_EMULATOR
#include <CL/sycl/INTEL/fpga_extensions.hpp>
#endif

using namespace sycl;
using namespace std;

// Parametres for this example.
#define ftype float
#define ftypePML float
#define ftype_sigma double

#include <stdio.h>

template <typename type_data>
class data3d_sycl;

// Create an exception handler for asynchronous SYCL exceptions
static auto exception_handler = [](sycl::exception_list e_list) {
	for (std::exception_ptr const& e : e_list) {
		try {
			std::rethrow_exception(e);
		}
		catch (std::exception const& e) {
#if _DEBUG
			std::cout << "Failure" << std::endl;
#endif
			std::terminate();
		}
	}
};


//************************************
// Initialize the array from 0 to array_size - 1
//************************************
void InitializeArray(double* a, size_t size) {
	for (size_t i = 0; i < size; i++) a[i] = 0.0;
}

//************************************
// Demonstrate vector add both in sequential on CPU and in parallel on device.
//************************************



int main(int argc, char* argv[]) {
    
    cpu_selector cpu_selector;

    double n = 4;
    ftype T = 8. * M_PI, dt = 0.005;
    constexpr int delta_x = 8, delta_y = 8, delta_z = 8;
    ftype_sigma sigma_x = 46.5, sigma_y = 46.5;
    constexpr int Nx = 256, Ny = 128, Nz = 32;

    pair<ftype, ftype> ab(0.0, 4.0 * M_PI), cd(0.0, 16.0 * M_PI), fg(0.0, 16.0 * M_PI);
    ftype ax_transverse = ab.second * (ftype)3. / (ftype)4.;
    ftype ay_transverse = cd.second * (ftype)1. / (ftype)6.;
    ftype az_transverse = fg.second * (ftype)1. / (ftype)6.;

    ftype tp_x = ab.second / (ftype)6.;
    ftype tp_y = cd.second / (ftype)6.;
    ftype tp_z = fg.second / (ftype)6.;

    ftype lambda_ = (ftype)2. * (ftype)M_PI / (ftype)4.;
    ftype k_ = (ftype)2. * (ftype)M_PI / lambda_;
    ftype omega_ = (ftype)2. * (ftype)M_PI / lambda_;
    ftype A = (ftype)1.; //amplitude
    ftype x00 = (ab.second - ab.first) / (ftype)2.;
    ftype y00 = (cd.second - cd.first) / (ftype)2.;
    ftype z00 = (fg.second - fg.first) / (ftype)2.;
    ftype t0_x = (ftype)3. * tp_x;
    ftype t0_y = (ftype)3. * tp_y;
    ftype t0_z = (ftype)3. * tp_z;

    int index_start_x = delta_x + 1;
    int index_start_y = delta_y + 1;
    int index_start_z = delta_z + 1;

    int offset = 5;

	int Nt = (T / dt);
    int THREAD_COUNT = 64;
    int THREAD_COUNT_SPEC = 8;
    
    try {
        queue q(cpu_selector, exception_handler);

        // Print out the device information used for the kernel code.
        std::cout << "Running on device: "
            << q.get_device().get_info<info::device::name>() << "\n";
        std::cout << "T =  " << T << "  dt = " << dt << "  Nt = " << Nt << "\n"
            << "Nx = " << Nx << "  Ny = " << Ny << "  Nz = " << Nz << "\n"
            << "delta_x = " << delta_x << "  delta_y = " << delta_y << "  delta_z = " << delta_z << "\n"
            << "sigma_x = " << sigma_x << "  sigma_y = " << sigma_y << "  n = " << n << "\n";

        Fields<ftype> cube(q, Nx, Ny, Nz, delta_x, delta_y, delta_z);
        SplitFields<ftypePML> cube_split(q, Nx, Ny, Nz, delta_x, delta_y, delta_z);
        Sigma<ftype_sigma> sigma(q, Nx, Ny, Nz, delta_x, delta_y, delta_z);
        Coefficient<ftypePML> Coeff(q, Nx, Ny, Nz, delta_x, delta_y, delta_z);

        double* vec_energy = malloc_shared<double>(size_t(Nt + 1), q);

        //Calculating the values of some parameters

        ftype dx = (ab.second - ab.first) / (ftype)Nx, dy = (cd.second - cd.first) / (ftype)Ny,
            dz = (fg.second - fg.first) / (ftype)Nz;
        ftype dt_x = dt / (dx), dt_y = dt / (dy), dt_z = dt / (dz);
        ftypePML _1dx = (ftypePML)1.0 / dx, _1dy = (ftypePML)1.0 / dy, _1dz = (ftypePML)1.0 / dz;

        // Initialize grid 
        Initializing_cube_sycl_half_version<ftype, ftypePML>(cube, cube_split, Nx, Ny, Nz, delta_x, delta_y, delta_z);
        InitializeArray(vec_energy, Nt + 1);

        auto start = std::chrono::steady_clock::now();

        Initializing_Sigma_Coeff_sycl_half_version_2_0<ftypePML, ftype_sigma>(Coeff, sigma, Nx, Ny, Nz,
            delta_x, delta_y, delta_z, n, sigma_x, sigma_y, dt);


        //for (int i = 1; i < Nx + 2 * delta_x + 1; i++) {
        //    for (int j = 1; j < Ny + 2 * delta_y + 1; j++) {
        //        for (int k = 1; k < Nz + 2 * delta_z + 1; k++) {
        //            if (sigma.sigmaEx(i, j, k) != sigma_float.sigmaEx(i, j, k)) {
        //                cout << std::endl << "error sigma_Ex: " << sigma.sigmaEx(i, j, k) << "  " << sigma_float.sigmaEx(i, j, k) << std::endl;
        //            }
        //            if (sigma.sigmaEy(i, j, k) != sigma_float.sigmaEy(i, j, k)) {
        //                cout << std::endl << "error sigma_Ey: " << sigma.sigmaEy(i, j, k) << "  " << sigma_float.sigmaEy(i, j, k) << std::endl;
        //            }
        //            if (sigma.sigmaEz(i, j, k) != sigma_float.sigmaEz(i, j, k)) {
        //                cout << std::endl << "error sigma_Ez: " << sigma.sigmaEz(i, j, k) << "  " << sigma_float.sigmaEz(i, j, k) << std::endl;
        //            }
        //            if (sigma.sigmaBx(i, j, k) != sigma_float.sigmaBx(i, j, k)) {
        //                cout << std::endl << "error sigma_Bx: " << sigma.sigmaBx(i, j, k) << "  " << sigma_float.sigmaBx(i, j, k) << std::endl;
        //            }
        //            if (sigma.sigmaBy(i, j, k) != sigma_float.sigmaBy(i, j, k)) {
        //                cout << std::endl << "error sigma_By: " << sigma.sigmaBy(i, j, k) << "  " << sigma_float.sigmaBy(i, j, k) << std::endl;
        //            }
        //            if (sigma.sigmaBz(i, j, k) != sigma_float.sigmaBz(i, j, k)) {
        //                cout << std::endl << "error sigma_Bz: " << sigma.sigmaBz(i, j, k) << "  " << sigma_float.sigmaBz(i, j, k) << std::endl;
        //            }
        //        }
        //    }
        //}
        /*
        for (int i = 1; i < Nx + 2 * delta_x + 1; i++) {
            for (int j = 1; j < Ny + 2 * delta_y + 1; j++) {
                for (int k = 1; k < Nz + 2 * delta_z + 1; k++) {
                    if (Coeff.Exy1(i, j, k) != Coeff_float.Exy1(i, j, k)) {
                        cout << std::endl << "error coeff_Exy1: " << Coeff.Exy1(i, j, k) << "  " << Coeff_float.Exy1(i, j, k) << std::endl;
                    }
                    if (Coeff.Exz1(i, j, k) != Coeff_float.Exz1(i, j, k)) {
                        cout << std::endl << "error coeff_Exz1: " << Coeff.Exz1(i, j, k) << "  " << Coeff_float.Exz1(i, j, k) << std::endl;
                    }
                    if (Coeff.Eyx1(i, j, k) != Coeff_float.Eyx1(i, j, k)) {
                        cout << std::endl << "error coeff_Eyx1: " << Coeff.Eyx1(i, j, k) << "  " << Coeff_float.Eyx1(i, j, k) << std::endl;
                    }
                    if (Coeff.Eyz1(i, j, k) != Coeff_float.Eyz1(i, j, k)) {
                        cout << std::endl << "error coeff_Eyz1: " << Coeff.Eyz1(i, j, k) << "  " << Coeff_float.Eyz1(i, j, k) << std::endl;
                    }
                    if (Coeff.Ezx1(i, j, k) != Coeff_float.Ezx1(i, j, k)) {
                        cout << std::endl << "error coeff_Ezx1: " << Coeff.Ezx1(i, j, k) << "  " << Coeff_float.Ezx1(i, j, k) << std::endl;
                    }
                    if (Coeff.Ezy1(i, j, k) != Coeff_float.Ezy1(i, j, k)) {
                        cout << std::endl << "error coeff_Ezy1: " << Coeff.Ezy1(i, j, k) << "  " << Coeff_float.Ezy1(i, j, k) << std::endl;
                    }
                    if (Coeff.Bxy1(i, j, k) != Coeff_float.Bxy1(i, j, k)) {
                        cout << std::endl << "error coeff_Bxy1: " << Coeff.Bxy1(i, j, k) << "  " << Coeff_float.Bxy1(i, j, k) << std::endl;
                    }
                    if (Coeff.Bxz1(i, j, k) != Coeff_float.Bxz1(i, j, k)) {
                        cout << std::endl << "error coeff_Bxz1: " << Coeff.Bxz1(i, j, k) << "  " << Coeff_float.Bxz1(i, j, k) << std::endl;
                    }
                    if (Coeff.Byx1(i, j, k) != Coeff_float.Byx1(i, j, k)) {
                        cout << std::endl << "error coeff_Byx1: " << Coeff.Byx1(i, j, k) << "  " << Coeff_float.Byx1(i, j, k) << std::endl;
                    }
                    if (Coeff.Byz1(i, j, k) != Coeff_float.Byz1(i, j, k)) {
                        cout << std::endl << "error coeff_Byz1: " << Coeff.Byz1(i, j, k) << "  " << Coeff_float.Byz1(i, j, k) << std::endl;
                    }
                    if (Coeff.Bzx1(i, j, k) != Coeff_float.Bzx1(i, j, k)) {
                        cout << std::endl << "error coeff_Bzx1: " << Coeff.Bzx1(i, j, k) << "  " << Coeff_float.Bzx1(i, j, k) << std::endl;
                    }
                    if (Coeff.Bzy1(i, j, k) != Coeff_float.Bzy1(i, j, k)) {
                        cout << std::endl << "error coeff_Bzy1: " << Coeff.Bzy1(i, j, k) << "  " << Coeff_float.Bzy1(i, j, k) << std::endl;
                    }

                    if (Coeff.Exy2(i, j, k) != Coeff_float.Exy2(i, j, k)) {
                        cout << std::endl << "error coeff_Exy2: " << Coeff.Exy2(i, j, k) << "  " << Coeff_float.Exy2(i, j, k) << std::endl;
                    }
                    if (Coeff.Exz2(i, j, k) != Coeff_float.Exz2(i, j, k)) {
                        cout << std::endl << "error coeff_Exz2: " << Coeff.Exz2(i, j, k) << "  " << Coeff_float.Exz2(i, j, k) << std::endl;
                    }
                    if (Coeff.Eyx2(i, j, k) != Coeff_float.Eyx2(i, j, k)) {
                        cout << std::endl << "error coeff_Eyx2: " << Coeff.Eyx2(i, j, k) << "  " << Coeff_float.Eyx2(i, j, k) << std::endl;
                    }
                    if (Coeff.Eyz2(i, j, k) != Coeff_float.Eyz2(i, j, k)) {
                        cout << std::endl << "error coeff_Eyz2: " << Coeff.Eyz2(i, j, k) << "  " << Coeff_float.Eyz2(i, j, k) << std::endl;
                    }
                    if (Coeff.Ezx2(i, j, k) != Coeff_float.Ezx2(i, j, k)) {
                        cout << std::endl << "error coeff_Ezx2: " << Coeff.Ezx2(i, j, k) << "  " << Coeff_float.Ezx2(i, j, k) << std::endl;
                    }
                    if (Coeff.Ezy2(i, j, k) != Coeff_float.Ezy2(i, j, k)) {
                        cout << std::endl << "error coeff_Ezy2: " << Coeff.Ezy2(i, j, k) << "  " << Coeff_float.Ezy2(i, j, k) << std::endl;
                    }
                    if (Coeff.Bxy2(i, j, k) != Coeff_float.Bxy2(i, j, k)) {
                        cout << std::endl << "error coeff_Bxy2: " << Coeff.Bxy2(i, j, k) << "  " << Coeff_float.Bxy2(i, j, k) << std::endl;
                    }
                    if (Coeff.Bxz2(i, j, k) != Coeff_float.Bxz2(i, j, k)) {
                        cout << std::endl << "error coeff_Bxz2: " << Coeff.Bxz2(i, j, k) << "  " << Coeff_float.Bxz2(i, j, k) << std::endl;
                    }
                    if (Coeff.Byx2(i, j, k) != Coeff_float.Byx2(i, j, k)) {
                        cout << std::endl << "error coeff_Byx2: " << Coeff.Byx2(i, j, k) << "  " << Coeff_float.Byx2(i, j, k) << std::endl;
                    }
                    if (Coeff.Byz2(i, j, k) != Coeff_float.Byz2(i, j, k)) {
                        cout << std::endl << "error coeff_Byz2: " << Coeff.Byz2(i, j, k) << "  " << Coeff_float.Byz2(i, j, k) << std::endl;
                    }
                    if (Coeff.Bzx2(i, j, k) != Coeff_float.Bzx2(i, j, k)) {
                        cout << std::endl << "error coeff_Bzx2: " << Coeff.Bzx2(i, j, k) << "  " << Coeff_float.Bzx2(i, j, k) << std::endl;
                    }
                    if (Coeff.Bzy2(i, j, k) != Coeff_float.Bzy2(i, j, k)) {
                        cout << std::endl << "error coeff_Bzy2: " << Coeff.Bzy2(i, j, k) << "  " << Coeff_float.Bzy2(i, j, k) << std::endl;
                    }
                }
            }
        }
        */

        // FDTD method in DPC++.
        std::cout << Nt << std::endl;

        for (int it = 0; it < Nt; ++it) {

            if (it % 100 == 0)
                cout << it << std::endl;

            //bound conditions
            ftype t1 = dt * (ftype)it;
            ftype t2 = t1 + (ftype)0.5 * dt;
            ftype x1 = (ftype)(ab.first + dx * (ftype)offset + (ftype)0.5 * dx);
            ftype x2 = (ftype)(ab.first + dx * (ftype)offset);

            q.submit([&](sycl::handler& h) {
                
                //std::cout << "Success" << std::endl;
                h.parallel_for(range<2>(Ny, Nz), [=](auto index) {
                    int j = index[0];
                    int k = index[1];
                    auto kernel_cube = cube;

                    Calculate_Currents_sycl_half_version(kernel_cube, dy,  dz, cd, fg, offset, index_start_x, index_start_y, index_start_z, j, k,  A,
                         t1,  t2,  t0_x,  tp_x,  y00,  x1,  x2,  ay_transverse,  omega_,  k_);
                    });
                

                }).wait_and_throw();

            vec_energy[it] = CalculateEnergy_sycl_half_version(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
            std::cout << it << "   " << vec_energy[it] << std::endl;

            //XYZ
            //run calculate electric fields
            //internal area

            q.submit([&](sycl::handler& h) {
                //sycl::stream out(4096, 4096, h);
                    
                h.parallel_for(range<3>(THREAD_COUNT, THREAD_COUNT, THREAD_COUNT), [=](auto index) {
                    //out << sycl::scientific;
                    for (int i = index[0] + delta_x + 1; i < Nx + delta_x + 1; i += THREAD_COUNT)
                        for (int j = index[1] + delta_y + 1; j < Ny + delta_y + 1; j += THREAD_COUNT)
                            for (int k = index[2] + delta_z + 1; k < Nz + delta_z + 1; k += THREAD_COUNT) {
                                Update_electric_field_sycl_half_version<ftype>(cube, dt_x, dt_y, dt_z, i, j, k);
                            }
                    // out << "success";
                    /*vec_EX[it] = cube.Ex(169, 45, 40);
                    vec_EY[it] = cube.Ey(169, 45, 40);
                    vec_EZ[it] = cube.Ez(169, 45, 40);*/

                    });          
             
                }).wait_and_throw();

            //edge from 0,0,0 by XY
            q.submit([&](sycl::handler& h) {
                sycl::stream out(4096, 4096, h);
                h.parallel_for(range<3>(THREAD_COUNT, THREAD_COUNT, THREAD_COUNT_SPEC), [=](auto index) {
                    //#define _MM_SET_FLUSH_ZERO_MODE _MM_FLUSH_ZERO_ON
                    //#define _MM_SET_DENORMALS_ZERO_MODE _MM_DENORMALS_ZERO_ON
                    out << sycl::scientific;
                    for (int i = index[0] + 1; i < Nx + 2 * delta_x + 1; i += THREAD_COUNT)
                        for (int j = index[1] + 1; j < Ny + 2 * delta_y + 1; j += THREAD_COUNT)
                            for (int k = index[2] + 1; k < delta_z + 1; k += THREAD_COUNT_SPEC) {
                                Update_electric_field_PML_sycl_half_version<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
                            }
                    });
                }).wait_and_throw(); //std::cout<<"Success3"<<std::endl;

            //edge from 0,0,0 by XZ
            q.submit([&](sycl::handler& h) {
                sycl::stream out(4096, 4096, h);
                h.parallel_for(range<3>(THREAD_COUNT, THREAD_COUNT_SPEC, THREAD_COUNT), [=](auto index) {
                    //#define  _MM_SET_FLUSH_ZERO_MODE _MM_FLUSH_ZERO_ON
                    //#define  _MM_SET_DENORMALS_ZERO_MODE _MM_DENORMALS_ZERO_ON
                    out << sycl::scientific;
                    for (int i = index[0] + 1; i < Nx + 2 * delta_x + 1; i += THREAD_COUNT)
                        for (int j = index[1] + 1; j < delta_y + 1; j += THREAD_COUNT_SPEC)
                            for (int k = index[2] + delta_z + 1; k < Nz + 2 * delta_z + 1; k += THREAD_COUNT) {
                                Update_electric_field_PML_sycl_half_version<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
                            }
                    });
                }).wait_and_throw(); //std::cout<<"Success4"<<std::endl;

            //edge from 0,0,0 by YZ
            q.submit([&](sycl::handler& h) {
                sycl::stream out(4096, 4096, h);
                h.parallel_for(range<3>(THREAD_COUNT_SPEC, THREAD_COUNT, THREAD_COUNT), [=](auto index) {
                    //#define  _MM_SET_FLUSH_ZERO_MODE _MM_FLUSH_ZERO_ON
                    //#define  _MM_SET_DENORMALS_ZERO_MODE _MM_DENORMALS_ZERO_ON
                    out << sycl::scientific;
                    for (int i = index[0] + 1; i < delta_x + 1; i += THREAD_COUNT_SPEC)
                        for (int j = index[1] + delta_y + 1; j < Ny + 2 * delta_y + 1; j += THREAD_COUNT)
                            for (int k = index[2] + delta_z + 1; k < Nz + 2 * delta_z + 1; k += THREAD_COUNT) {
                                Update_electric_field_PML_sycl_half_version<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
                            }
                    });
                }).wait_and_throw(); //std::cout<<"Success5"<<std::endl;

            //edge from n,0,0 by YZ
            q.submit([&](sycl::handler& h) {
                sycl::stream out(4096, 4096, h);
                h.parallel_for(range<3>(THREAD_COUNT_SPEC, THREAD_COUNT, THREAD_COUNT), [=](auto index) {
                    //#define  _MM_SET_FLUSH_ZERO_MODE _MM_FLUSH_ZERO_ON
                        //#define _MM_SET_DENORMALS_ZERO_MODE _MM_DENORMALS_ZERO_ON
                    out << sycl::scientific;
                    for (int i = index[0] + Nx + delta_x + 1; i < Nx + 2 * delta_x + 1; i += THREAD_COUNT_SPEC)
                        for (int j = index[1] + delta_y + 1; j < Ny + 2 * delta_y + 1; j += THREAD_COUNT)
                            for (int k = index[2] + delta_z + 1; k < Nz + 2 * delta_z + 1; k += THREAD_COUNT) {
                                Update_electric_field_PML_sycl_half_version<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
                            }
                    });
                }).wait_and_throw(); //std::cout<<"Success6"<<std::endl;

            //edge from 0,m,0 by XZ
            q.submit([&](sycl::handler& h) {
                sycl::stream out(4096, 4096, h);
                h.parallel_for(range<3>(THREAD_COUNT, THREAD_COUNT_SPEC, THREAD_COUNT), [=](auto index) {
                    //#define _MM_SET_FLUSH_ZERO_MODE _MM_FLUSH_ZERO_ON
                    //#define _MM_SET_DENORMALS_ZERO_MODE _MM_DENORMALS_ZERO_ON
                    out << sycl::scientific;
                    for (int i = index[0] + delta_x + 1; i < Nx + delta_x + 1; i += THREAD_COUNT)
                        for (int j = index[1] + Ny + delta_y + 1; j < Ny + 2 * delta_y + 1; j += THREAD_COUNT_SPEC)
                            for (int k = index[2] + delta_z + 1; k < Nz + 2 * delta_z + 1; k += THREAD_COUNT) {
                                Update_electric_field_PML_sycl_half_version<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
                            }
                    });
                }).wait_and_throw(); //std::cout<<"Success7"<<std::endl;

            //edge from 0,0,k by XY
            q.submit([&](sycl::handler& h) {
                sycl::stream out(4096, 4096, h);
                h.parallel_for(range<3>(THREAD_COUNT, THREAD_COUNT, THREAD_COUNT_SPEC), [=](auto index) {
                    //#define _MM_SET_FLUSH_ZERO_MODE _MM_FLUSH_ZERO_ON 
                    //#define _MM_SET_DENORMALS_ZERO_MODE _MM_DENORMALS_ZERO_ON
                    out << sycl::scientific;
                    for (int i = index[0] + delta_x + 1; i < Nx + delta_x + 1; i += THREAD_COUNT)
                        for (int j = index[1] + delta_y + 1; j < Ny + delta_y + 1; j += THREAD_COUNT)
                            for (int k = index[2] + Nz + delta_z + 1; k < Nz + 2 * delta_z + 1; k += THREAD_COUNT_SPEC) {
                                Update_electric_field_PML_sycl_half_version<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
                            }
                    });
                }).wait_and_throw();
            // std::cout<<"Success8"<<std::endl;

            //time_vec[i_time++] = std::chrono::steady_clock::now();

            //run calculate magnetic fields
            //internal area
            q.submit([&](sycl::handler& h) {
                //sycl::stream out(4096, 4096, h);
                h.parallel_for(range<3>(THREAD_COUNT, THREAD_COUNT, THREAD_COUNT), [=](auto index) {
                    //#define _MM_SET_FLUSH_ZERO_MODE _MM_FLUSH_ZERO_ON
                    //#define _MM_SET_DENORMALS_ZERO_MODE _MM_DENORMALS_ZERO_ON
                    // out << sycl::scientific;
                    for (int i = index[0] + delta_x + 1; i < Nx + delta_x + 1; i += THREAD_COUNT)
                        for (int j = index[1] + delta_y + 1; j < Ny + delta_y + 1; j += THREAD_COUNT)
                            for (int k = index[2] + delta_z + 1; k < Nz + delta_z + 1; k += THREAD_COUNT) {
                                Update_magnetic_field_sycl_half_version<ftype>(cube, dt_x, dt_y, dt_z, i, j, k);
                            }
                        
                    //vec_BX[it] = cube.Bx(169, 45, 40);
                    //vec_BY[it] = cube.By(169, 45, 40);
                    //vec_BZ[it] = cube.Bz(169, 45, 40);

                    });
                //OutputDenormalMagnetic(cube, 169, 45, 40, it, vec_BX, vec_BY, vec_BZ);
                }).wait_and_throw();

            //edge from 0,0,0 by XY
            q.submit([&](sycl::handler& h) {
                //sycl::stream out(4096, 4096, h);
                h.parallel_for(range<3>(THREAD_COUNT, THREAD_COUNT, THREAD_COUNT_SPEC), [=](auto index) {
                    //#define _MM_SET_FLUSH_ZERO_MODE _MM_FLUSH_ZERO_ON
                    //#define _MM_SET_DENORMALS_ZERO_MODE _MM_DENORMALS_ZERO_ON
                    //out << sycl::scientific;
                    for (int i = index[0] + 1; i < Nx + 2 * delta_x + 1; i += THREAD_COUNT)
                        for (int j = index[1] + 1; j < Ny + 2 * delta_y + 1; j += THREAD_COUNT)
                            for (int k = index[2] + 1; k < delta_z + 1; k += THREAD_COUNT_SPEC) {
                                Update_magnetic_field_PML_sycl_half_version<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
                            }
                    });
                }).wait_and_throw();

            //edge from 0,0,0 by XZ
            q.submit([&](sycl::handler& h) {
                //sycl::stream out(4096, 4096, h);
                h.parallel_for(range<3>(THREAD_COUNT, THREAD_COUNT_SPEC, THREAD_COUNT), [=](auto index) {
                    //#define  _MM_SET_FLUSH_ZERO_MODE _MM_FLUSH_ZERO_ON
                    //#define  _MM_SET_DENORMALS_ZERO_MODE _MM_DENORMALS_ZERO_ON
                    //out << sycl::scientific;
                    for (int i = index[0] + 1; i < Nx + 2 * delta_x + 1; i += THREAD_COUNT)
                        for (int j = index[1] + 1; j < delta_y + 1; j += THREAD_COUNT_SPEC)
                            for (int k = index[2] + delta_z + 1; k < Nz + 2 * delta_z + 1; k += THREAD_COUNT) {
                                Update_magnetic_field_PML_sycl_half_version<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
                            }
                    });
                }).wait_and_throw();

            //edge from 0,0,0 by YZ
            q.submit([&](sycl::handler& h) {
                //sycl::stream out(4096, 4096, h);
                h.parallel_for(range<3>(THREAD_COUNT_SPEC, THREAD_COUNT, THREAD_COUNT), [=](auto index) {
                    //#define _MM_SET_FLUSH_ZERO_MODE _MM_FLUSH_ZERO_ON
                    //#define _MM_SET_DENORMALS_ZERO_MODE _MM_DENORMALS_ZERO_ON
                    //out << sycl::scientific;
                    for (int i = index[0] + 1; i < delta_x + 1; i += THREAD_COUNT_SPEC)
                        for (int j = index[1] + delta_y + 1; j < Ny + 2 * delta_y + 1; j += THREAD_COUNT)
                            for (int k = index[2] + delta_z + 1; k < Nz + 2 * delta_z + 1; k += THREAD_COUNT) {
                                Update_magnetic_field_PML_sycl_half_version<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
                            }
                    });
                }).wait_and_throw();

            //edge from n,0,0 by YZ
            q.submit([&](sycl::handler& h) {
                //sycl::stream out(4096, 4096, h);
                h.parallel_for(range<3>(THREAD_COUNT_SPEC, THREAD_COUNT, THREAD_COUNT), [=](auto index) {
                    //#define  _MM_SET_FLUSH_ZERO_MODE _MM_FLUSH_ZERO_ON
                        //#define _MM_SET_DENORMALS_ZERO_MODE _MM_DENORMALS_ZERO_ON
                    //out << sycl::scientific;
                    for (int i = index[0] + Nx + delta_x + 1; i < Nx + 2 * delta_x + 1; i += THREAD_COUNT_SPEC)
                        for (int j = index[1] + delta_y + 1; j < Ny + 2 * delta_y + 1; j += THREAD_COUNT)
                            for (int k = index[2] + delta_z + 1; k < Nz + 2 * delta_z + 1; k += THREAD_COUNT) {
                                Update_magnetic_field_PML_sycl_half_version<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
                            }
                    });
                }).wait_and_throw();

            //edge from 0,m,0 by XZ
            q.submit([&](sycl::handler& h) {
                //sycl::stream out(4096, 4096, h);
                h.parallel_for(range<3>(THREAD_COUNT, THREAD_COUNT_SPEC, THREAD_COUNT), [=](auto index) {
                    //#define _MM_SET_FLUSH_ZERO_MODE _MM_FLUSH_ZERO_ON
                    //#define _MM_SET_DENORMALS_ZERO_MODE _MM_DENORMALS_ZERO_ON
                    //out << sycl::scientific;
                    for (int i = index[0] + delta_x + 1; i < Nx + delta_x + 1; i += THREAD_COUNT)
                        for (int j = index[1] + Ny + delta_y + 1; j < Ny + 2 * delta_y + 1; j += THREAD_COUNT_SPEC)
                            for (int k = index[2] + delta_z + 1; k < Nz + 2 * delta_z + 1; k += THREAD_COUNT) {
                                Update_magnetic_field_PML_sycl_half_version<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
                            }
                    });
                }).wait_and_throw();

            //edge from 0,0,k by XY
            q.submit([&](sycl::handler& h) {
                //sycl::stream out(4096, 4096, h);
                h.parallel_for(range<3>(THREAD_COUNT, THREAD_COUNT, THREAD_COUNT_SPEC), [=](auto index) {
                    //out << sycl::scientific;
                    for (int i = index[0] + delta_x + 1; i < Nx + delta_x + 1; i += THREAD_COUNT)
                        for (int j = index[1] + delta_y + 1; j < Ny + delta_y + 1; j += THREAD_COUNT)
                            for (int k = index[2] + Nz + delta_z + 1; k < Nz + 2 * delta_z + 1; k += THREAD_COUNT_SPEC) {
                                Update_magnetic_field_PML_sycl_half_version<ftype, ftypePML>(cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k);
                            }
                    });
                }).wait_and_throw();

                //std::cout << it << " " << cube.Bx(14, 25, 41) << std::endl;
            //if (it == 477) {

            //    FILE* f;
            //    f = fopen("D:\\work\\Project1\\Ey_477.dat", "ab+");
            //    float temp;
            //    ofstream numbEy("Ey_test.txt");


            //    for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
            //        for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
            //            for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
            //            {
            //                temp = cl::sycl::detail::half2Float(cube.Ey(i, j, k));
            //                fwrite(&temp, sizeof(float), 1, f);
            //                numbEy << temp << "  ";
            //                if (cube.Ey(i, j, k) != 0.) {
            //                    std::cout << i << " " << j << " " << k << " = " << temp  << " == " <<cube.Ey(i, j, k) << std::endl;
            //                }
            //            }
            //    fseek(f, 0, SEEK_SET);
            //    std::cout << "FROM FILE" << std::endl;
            //    for (int i = 1; i < 49; i++) {
            //        for (int j = 1; j < 49; j++) {
            //            for (int k = 1; k < 49; k++) {
            //                fread(&temp, sizeof(float), 1, f);
            //                if (fabs(temp) > 0.) {
            //                    printf("%i %i %i = %f \n", i, j, k, temp);
            //                }
            //            }
            //        }
            //    }
            //    fclose(f);
            //    numbEy.close();
            //    std::cout << "BY" << std::endl;
            //    FILE* f1;
            //    f1 = fopen("D:\\work\\Project1\\Bz_477.dat", "ab+");

            //    ofstream numbBz("Bz_test.txt");

            //    for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
            //        for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
            //            for (int k = 1; k < Nz + 2 * delta_z + 1; k++)
            //            {
            //                temp = cl::sycl::detail::half2Float(cube.Bz(i, j, k));
            //                fwrite(&temp, sizeof(float), 1, f1);

            //                numbBz << temp << "  ";
            //                if (cube.Bz(i, j, k) != 0.) {
            //                    std::cout << i << " " << j << " " << k << " = " << temp << " == "<<cube.Bz(i, j, k) << std::endl;
            //                }
            //            }
            //    std::cout << "FROM FILE" << std::endl;
            //    for (int i = 1; i < 49; i++) {
            //        for (int j = 1; j < 49; j++) {
            //            for (int k = 1; k < 49; k++) {
            //                fread(&temp, sizeof(float), 1, f1);
            //                if (fabs(temp) > 0.) {
            //                    printf("%i %i %i = %f \n", i, j, k, temp);
            //                }
            //            }
            //        }
            //    }

            //    fclose(f1);
            //    numbBz.close();
            //    }
                
        }
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << "elapsed time = " << elapsed_seconds.count() << "\n";

        vec_energy[Nt] = CalculateEnergy_sycl_half_version(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
        double result = vec_energy[Nt] / vec_energy[Nt/2];
        std::cout << "R = " << result << std::endl;

        cube.Delete(q);
        cube_split.Delete(q);
        sigma.Delete(q);
        Coeff.Delete(q);
        sycl::free(vec_energy, q);
        
    }
	catch (sycl::exception const& e) {
		std::cout << "An exception is caught while adding two vectors.\n" << e.what();
		std::terminate();
	}

	std::cout << "Vector add successfully completed on device.\n";
	return 0;
}
