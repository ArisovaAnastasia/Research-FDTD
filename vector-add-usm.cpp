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
#define ftype half
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
void InitializeArray(float* a, size_t size) {
    for (size_t i = 0; i < size; i++) a[i] = 0.0;
}

int main(int argc, char* argv[]) {

    //cpu_selector gpu_selector;

    float n = 4;
    ftype T = ftype(8.0f * M_PI), dt = ftype(0.005f);
    constexpr int delta_x = 8, delta_y = 0, delta_z = 8;
    ftype_sigma sigma_x = 46.5, sigma_y = 46.5;
    constexpr int Nx = 256, Ny = 1, Nz = 128;

    pair<ftype, ftype> ab(ftype(0.0f), ftype(4.0f * M_PI)),
        cd(ftype(0.0), ftype(16.0f * M_PI)),
        fg(ftype(0.0), ftype(16.0f * M_PI));


    int Nt = (float(T) / float(dt));
    int THREAD_COUNT = 64;
    int THREAD_COUNT_SPEC = 7;

    try {
        sycl::queue q = sycl::queue{ sycl::gpu_selector{}, sycl::async_handler{} };

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

        float* vec_energy = malloc_shared<float>(size_t(Nt + 1), q);

        ftype dx = ftype((float(ab.second) - float(ab.first)) / float(Nx)),
            dy = ftype((float(cd.second) - float(cd.first)) / float(Ny)),
            dz = ftype((float(fg.second) - float(fg.first)) / float(Nz));
        ftype dt_x = ftype(float(dt) / float(dx)),
            dt_y = ftype(float(dt) / float(dy)),
            dt_z = ftype(float(dt) / float(dz));
        ftypePML _1dx = (ftypePML)1.0 / dx, _1dy = (ftypePML)1.0 / dy, _1dz = (ftypePML)1.0 / dz;

        // Initialize grid 
        Initializing_cube_sycl_half_version<ftype, ftypePML>(cube, cube_split, Nx, Ny, Nz, delta_x, delta_y, delta_z);
        InitializeArray(vec_energy, Nt + 1);

        auto start = std::chrono::steady_clock::now();

        Initializing_Sigma_Coeff_sycl_half_version_2_0<ftypePML, ftype_sigma>(Coeff, sigma, Nx, Ny, Nz,
            delta_x, delta_y, delta_z, n, sigma_x, sigma_y, dt);

        // FDTD method in DPC++.
        //std::cout << Nt << std::endl;

        ofstream numbEy("Split_Ey_half.csv");
        ofstream numb("Energy_float.csv");

        for (int it = 0; it < Nt; ++it) {

            if (it % 100 == 0)
                cout << it << std::endl;

            Calculate_Currents_sycl_half_version<ftype>('X', q, cube, it, Nx, Ny, Nz, delta_x, delta_y, delta_z, dx, dy, dz, dt, ab, cd, fg);
            
            if (sizeof(ftype) == 2) {
                vec_energy[it] = CalculateEnergy_sycl_half_version(q, cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
            }
            else {
                vec_energy[it] = CalculateEnergy_sycl_float_version(q, cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
            }

            std::cout << it << "   " << vec_energy[it] << std::endl;
            numb << it << ";" << vec_energy[it] << std::endl;
            numbEy << it << ";" << cube_split.Eyz(256,1,138) << std::endl;

     //XYZ
     //run calculate electric fields
     //internal area

            q.submit([&](sycl::handler& h) {
                h.parallel_for(range<3>(Nx, Ny, Nz), [=](auto index) {
                    int i = index[0] + delta_x + 1;
                    int j = index[1] + delta_y + 1;
                    int k = index[2] + delta_z + 1;

                    Update_electric_domain_half_version<ftype>(cube, dt_x, dt_y, dt_z, i, j, k);
                    });
                }).wait_and_throw();


                if (delta_x != 0) {
                    q.submit([&](sycl::handler& h) {
                        sycl::stream out(4096, 4096, h);
                        h.parallel_for(range<2>(Ny, Nz), [=](auto index) {
                            int i = delta_x;
                            int j = index[0] + delta_y + 1;
                            int k = index[1] + delta_z + 1;

                            Update_electric_PML_BOUND_float_version<ftype, ftypePML>
                                (cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k, out);
                            });
                        }).wait_and_throw();
                }

                if (delta_z != 0) {
                    q.submit([&](sycl::handler& h) {
                        sycl::stream out(4096, 4096, h);
                        h.parallel_for(range<2>(Nx, Ny), [=](auto index) {
                            int i = index[0] + delta_x + 1;
                            int j = index[1] + delta_y + 1;
                            int k = delta_z;

                            Update_electric_PML_BOUND_float_version<ftype, ftypePML>
                                (cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k, out);
                            });
                        }).wait_and_throw();
                }

                if (delta_y != 0) {
                    q.submit([&](sycl::handler& h) {
                        sycl::stream out(4096, 4096, h);
                        h.parallel_for(range<2>(Nx, Nz), [=](auto index) {
                            int i = index[0] + delta_x + 1;
                            int j = delta_y;
                            int k = index[1] + delta_z + 1;

                            Update_electric_PML_BOUND_float_version<ftype, ftypePML>
                                (cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k, out);
                            });
                        }).wait_and_throw();
                }

                Run_PML_electric<ftype, ftypePML>(cube_split, Coeff,
                    _1dx, _1dy, _1dz, Nx, Ny, Nz, delta_x, delta_y, delta_z);



                if (sizeof(ftype) == 2) {
                    if (delta_x != 0) {
                        q.submit([&](sycl::handler& h) {
                            h.parallel_for(range<2>(Ny, Nz), [=](auto index) {
                                int i = delta_x;
                                int j = index[0] + delta_y + 1;
                                int k = index[1] + delta_z + 1;

                                CubeSplit2Cube_electric_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);

                                i = Nx + delta_x + 1;

                                CubeSplit2Cube_electric_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);
                                });
                            }).wait_and_throw();
                    }

                    if (delta_z != 0) {
                        q.submit([&](sycl::handler& h) {
                            h.parallel_for(range<2>(Nx, Ny), [=](auto index) {
                                int i = index[0] + delta_x + 1;
                                int j = index[1] + delta_y + 1;
                                int k = delta_z;

                                CubeSplit2Cube_electric_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);

                                k = Nz + delta_z + 1;

                                CubeSplit2Cube_electric_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);
                                });
                            }).wait_and_throw();
                    }

                    if (delta_y != 0) {
                        q.submit([&](sycl::handler& h) {
                            h.parallel_for(range<2>(Nx, Nz), [=](auto index) {
                                int i = index[0] + delta_x + 1;
                                int j = delta_y;
                                int k = index[1] + delta_z + 1;

                                CubeSplit2Cube_electric_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);

                                j = Ny + delta_y + 1;

                                CubeSplit2Cube_electric_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);
                                });
                            }).wait_and_throw();
                    }
                }
                else {
                    if (delta_x != 0) {

                        q.submit([&](sycl::handler& h) {
                            h.parallel_for(range<2>(Ny, Nz), [=](auto index) {
                                int i = delta_x;
                                int j = index[0] + delta_y + 1;
                                int k = index[1] + delta_z + 1;

                                CubeSplit2Cube_electric_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);

                                i = Nx + delta_x + 1;

                                CubeSplit2Cube_electric_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);
                                });
                            }).wait_and_throw();
                    }

                    if (delta_z != 0) {
                        q.submit([&](sycl::handler& h) {
                            h.parallel_for(range<2>(Nx, Ny), [=](auto index) {
                                int i = index[0] + delta_x + 1;
                                int j = index[1] + delta_y + 1;
                                int k = delta_z;

                                CubeSplit2Cube_electric_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);

                                k = Nz + delta_z + 1;

                                CubeSplit2Cube_electric_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);
                                });
                            }).wait_and_throw();
                    }

                    if (delta_y != 0) {
                        q.submit([&](sycl::handler& h) {
                            h.parallel_for(range<2>(Nx, Nz), [=](auto index) {
                                int i = index[0] + delta_x + 1;
                                int j = delta_y;
                                int k = index[1] + delta_z + 1;

                                CubeSplit2Cube_electric_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);

                                j = Ny + delta_y + 1;

                                CubeSplit2Cube_electric_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);
                                });
                            }).wait_and_throw();
                    }
                }

    //run calculate magnetic fields
    //internal area
                q.submit([&](sycl::handler& h) {
                    sycl::stream out(4096, 4096, h);
                    h.parallel_for(range<3>(Nx, Ny, Nz), [=](auto index) {
                        int i = index[0] + delta_x + 1;
                        int j = index[1] + delta_y + 1;
                        int k = index[2] + delta_z + 1;
                        Update_magnetic_domain_half_version<ftype>(cube, dt_x, dt_y, dt_z, i, j, k);
                        });
                    //OutputDenormalMagnetic(cube, 169, 45, 40, it, vec_BX, vec_BY, vec_BZ);
                    }).wait_and_throw(); //std::cout << "Success internal area magnetic" << std::endl;

                    if (delta_x != 0) {
                        q.submit([&](sycl::handler& h) {
                            sycl::stream out(4096, 4096, h);
                            h.parallel_for(range<2>(Ny, Nz), [=](auto index) {
                                int i = Nx + delta_x + 1;
                                int j = index[0] + delta_y + 1;
                                int k = index[1] + delta_z + 1;

                                Update_magnetic_PML_BOUND_float_version<ftype, ftypePML>
                                    (cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k, out);
                                });
                            }).wait_and_throw();
                    }

                    if (delta_z != 0) {
                        q.submit([&](sycl::handler& h) {
                            sycl::stream out(4096, 4096, h);
                            h.parallel_for(range<2>(Nx, Ny), [=](auto index) {
                                int i = index[0] + delta_x + 1;
                                int j = index[1] + delta_y + 1;
                                int k = Nz + delta_z + 1;

                                Update_magnetic_PML_BOUND_float_version<ftype, ftypePML>
                                    (cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k, out);
                                });
                            }).wait_and_throw();
                    }

                    if (delta_y != 0) {
                        q.submit([&](sycl::handler& h) {
                            sycl::stream out(4096, 4096, h);
                            h.parallel_for(range<2>(Nx, Nz), [=](auto index) {
                                int i = index[0] + delta_x + 1;
                                int j = Ny + delta_y + 1;
                                int k = index[1] + delta_z + 1;

                                Update_magnetic_PML_BOUND_float_version<ftype, ftypePML>
                                    (cube, cube_split, Coeff, _1dx, _1dy, _1dz, i, j, k, out);
                                });
                            }).wait_and_throw();
                    }

                    Run_PML_magnetic<ftype, ftypePML>(cube_split, Coeff,
                        _1dx, _1dy, _1dz, Nx, Ny, Nz, delta_x, delta_y, delta_z);

        //Run_CubeSplit2Cube_magnetic<ftype, ftypePML>(cube, cube_split, Nx, Ny, Nz, delta_x, delta_y, delta_z);


                    if (sizeof(ftype) == 2) {
                        if (delta_x != 0) {

                            q.submit([&](sycl::handler& h) {
                                h.parallel_for(range<2>(Ny, Nz), [=](auto index) {
                                    int i = delta_x;
                                    int j = index[0] + delta_y + 1;
                                    int k = index[1] + delta_z + 1;

                                    CubeSplit2Cube_magnetic_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);

                                    i = Nx + delta_x + 1;

                                    CubeSplit2Cube_magnetic_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);
                                    });
                                }).wait_and_throw();
                        }

                        if (delta_z != 0) {
                            q.submit([&](sycl::handler& h) {
                                h.parallel_for(range<2>(Nx, Ny), [=](auto index) {
                                    int i = index[0] + delta_x + 1;
                                    int j = index[1] + delta_y + 1;
                                    int k = delta_z;

                                    CubeSplit2Cube_magnetic_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);

                                    k = Nz + delta_z + 1;

                                    CubeSplit2Cube_magnetic_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);
                                    });
                                }).wait_and_throw();
                        }

                        if (delta_y != 0) {
                            q.submit([&](sycl::handler& h) {
                                h.parallel_for(range<2>(Nx, Nz), [=](auto index) {
                                    int i = index[0] + delta_x + 1;
                                    int j = delta_y;
                                    int k = index[1] + delta_z + 1;

                                    CubeSplit2Cube_magnetic_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);

                                    j = Ny + delta_y + 1;

                                    CubeSplit2Cube_magnetic_half_version<ftype, ftypePML>(cube, cube_split, i, j, k);
                                    });
                                }).wait_and_throw();
                        }
                    }
                    else {
                        if (delta_x != 0) {

                            q.submit([&](sycl::handler& h) {
                                h.parallel_for(range<2>(Ny, Nz), [=](auto index) {
                                    int i = delta_x;
                                    int j = index[0] + delta_y + 1;
                                    int k = index[1] + delta_z + 1;

                                    CubeSplit2Cube_magnetic_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);

                                    i = Nx + delta_x + 1;

                                    CubeSplit2Cube_magnetic_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);
                                    });
                                }).wait_and_throw();
                        }

                        if (delta_z != 0) {
                            q.submit([&](sycl::handler& h) {
                                h.parallel_for(range<2>(Nx, Ny), [=](auto index) {
                                    int i = index[0] + delta_x + 1;
                                    int j = index[1] + delta_y + 1;
                                    int k = delta_z;

                                    CubeSplit2Cube_magnetic_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);

                                    k = Nz + delta_z + 1;

                                    CubeSplit2Cube_magnetic_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);
                                    });
                                }).wait_and_throw();
                        }

                        if (delta_y != 0) {
                            q.submit([&](sycl::handler& h) {
                                h.parallel_for(range<2>(Nx, Nz), [=](auto index) {
                                    int i = index[0] + delta_x + 1;
                                    int j = delta_y;
                                    int k = index[1] + delta_z + 1;

                                    CubeSplit2Cube_magnetic_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);

                                    j = Ny + delta_y + 1;

                                    CubeSplit2Cube_magnetic_float_version<ftype, ftypePML>(cube, cube_split, i, j, k);
                                    });
                                }).wait_and_throw();
                        }
                    }



                    //std::cout << "Success Run_CubeSplit2Cube_magnetic" << std::endl;

                    //Update_PeriodicBound_magnetic<ftype, ftypePML>(cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);

        }
        numbEy.close();
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << "elapsed time = " << elapsed_seconds.count() << "\n";

        if (sizeof(ftype) == 2) {
            vec_energy[Nt] = CalculateEnergy_sycl_half_version(q, cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
        }
        else {
            vec_energy[Nt] = CalculateEnergy_sycl_float_version(q, cube, Nx, Ny, Nz, delta_x, delta_y, delta_z);
        }

        numb << Nt << ";" << vec_energy[Nt] << std::endl;
        numb.close();
        float result = vec_energy[Nt] / vec_energy[2600];
        std::cout << "R = " << result << std::endl;

        cube.Delete(q);
        cube_split.Delete(q);
        sigma.Delete(q);
        Coeff.Delete(q);
        sycl::free(vec_energy, q);

    }
    //catch (std::runtime_error const& e) {
    //	std::cout << "An exception is caught while adding two vectors.\n" << e.what();
    //	std::terminate();
    //}
    catch (sycl::exception const& e) {
        std::cout << "ERROR\n" << e.get_cl_code() << std::endl << e.what();
        std::terminate();
    }

    std::cout << "Vector add successfully completed on device.\n";
    return 0;
}
