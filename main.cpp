#include "Maxwells_template.h"
#include "Maxwells_two_dimensional_space.h"
#include "FDTD_YeeGrid.h"
#include "FDTD_YeeGrid_2dimen.h"
#include "PML_2dimen.h"
#include "PML_3dimen.h"
#include "half.hpp"
#include <iostream>
using namespace std;
using half_float::half;
using namespace half_float::literal;

int main() {

	//double R1 = 2.67169e-3;// значение R для опт сигмы = 18,86 при T=14pi
	//double R2 = 4.52533e-7;// значение R для опт сигмы = 46.5 при T=8pi
	 
	double T = 8. * M_PI;
	int delta = 8;
	double sigma_x = 46.5;
	double sigma_y = 1.45312;
	
	/*FDTD_three_dimen_with_PML<double, double>(128, 16, 16, 8.0_h * M_PI, 0.005_h, delta, sigma_x, sigma_y, "No Kahan", "Energy_double-double-sigma 45.5.csv", "Coeff_graph_dd_sigma_46.5.csv", "Sigma_graph_dd_sigma_46.5.csv");
	FDTD_three_dimen_with_PML<float, double>(128, 16, 16, 8.0_h * M_PI, 0.005_h, delta, sigma_x, sigma_y, "No Kahan", "Energy_float-double-sigma 45.5.csv", "Coeff_graph_fd_sigma_46.5.csv", "Sigma_graph_fd_sigma_46.5.csv");
	FDTD_three_dimen_with_PML<float, float>(128, 16, 16, 8.0_h * M_PI, 0.005_h, delta, sigma_x, sigma_y, "No Kahan", "Energy_float-float-sigma 45.5.csv", "Coeff_graph_ff_sigma_46.5.csv", "Sigma_graph_ff_sigma_46.5.csv");

	FDTD_three_dimen_with_PML<half, double>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, sigma_x, sigma_y, "No Kahan", "Energy_hd_sigma_46.5.csv", "Coeff_graph_hd_sigma_46.5.csv", "Sigma_graph_hd_sigma_46.5.csv");
	FDTD_three_dimen_with_PML<half, float>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, sigma_x, sigma_y, "No Kahan", "Energy_hf_sigma_46.5.csv", "Coeff_graph_hf_sigma_46.5.csv", "Sigma_graph_hf_sigma_46.5.csv");
	FDTD_three_dimen_with_PML<half, half>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, sigma_x, sigma_y, "No Kahan", "Energy_hh_sigma_46.5.csv", "Coeff_graph_hh_sigma_46.5.csv", "Sigma_graph_hh_sigma_46.5.csv");

	cout << "KAHAN --- NO_KAHAN" << endl;

	FDTD_three_dimen_with_PML<half, double>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, sigma_x, sigma_y, "Kahan", "Energy_hd_sigma_46.5_Kahan.csv", "Coeff_graph_hd_sigma_46.5_Kahan.csv", "Sigma_graph_hd_sigma_46.5_Kahan.csv");
	FDTD_three_dimen_with_PML<half, float>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, sigma_x, sigma_y, "Kahan", "Energy_hf_sigma_46.5_Kahan.csv", "Coeff_graph_hf_sigma_46.5_Kahan.csv", "Sigma_graph_hf_sigma_46.5_Kahan.csv");
	FDTD_three_dimen_with_PML<half, half>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, sigma_x, sigma_y, "Kahan", "Energy_hh_sigma_46.5_Kahan.csv", "Coeff_graph_hh_sigma_46.5_Kahan.csv", "Sigma_graph_hh_sigma_46.5_Kahan.csv");*/

	//cout << "KAHAN --- KAHAN" << endl;

	//FDTD_three_dimen_with_PML<half, double>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, sigma_x, sigma_y, "Kahan", "Energy_hd_sigma_46.5_Kahan2.csv", "Coeff_graph_hd_sigma_46.5_Kahan2.csv", "Sigma_graph_hd_sigma_46.5_Kahan2.csv");
	FDTD_three_dimen_with_PML<half, float>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, sigma_x, sigma_y, "Kahan", "Energy_hf_sigma_46.5_Kahan2.csv", "Coeff_graph_hf_sigma_46.5_Kahan2.csv", "Sigma_graph_hf_sigma_46.5_Kahan2.csv");

	//FDTD_three_dimen_with_PML<half, float>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, 20., sigma_y, "Kahan", "Energy_hf_sigma_20_Kahan.csv", ".csv", ".csv");
	//FDTD_three_dimen_with_PML<half, float>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, 26., sigma_y, "Kahan", "Energy_hf_sigma_26_Kahan.csv", ".csv", ".csv");
	//FDTD_three_dimen_with_PML<half, float>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, 30., sigma_y, "Kahan", "Energy_hf_sigma_30_Kahan.csv", ".csv", ".csv");
	//FDTD_three_dimen_with_PML<half, float>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, 35., sigma_y, "Kahan", "Energy_hf_sigma_35_Kahan.csv", ".csv", ".csv");
	//FDTD_three_dimen_with_PML<half, float>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, 40., sigma_y, "Kahan", "Energy_hf_sigma_40_Kahan.csv", ".csv", ".csv");
	//FDTD_three_dimen_with_PML<half, float>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, 47., sigma_y, "Kahan", "Energy_hf_sigma_47_Kahan.csv", ".csv", ".csv");
	//FDTD_three_dimen_with_PML<half, float>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, 50., sigma_y, "Kahan", "Energy_hf_sigma_50_Kahan.csv", ".csv", ".csv");
	//FDTD_three_dimen_with_PML<half, float>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, 55., sigma_y, "Kahan", "Energy_hf_sigma_55_Kahan.csv", ".csv", ".csv");
	//FDTD_three_dimen_with_PML<half, float>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, 60., sigma_y, "Kahan", "Energy_hf_sigma_60_Kahan.csv", ".csv", ".csv");

	//FDTD_three_dimen_with_PML<half, half>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, 20., sigma_y, "Kahan", "Energy_hh_sigma_20_Kahan.csv", ".csv", ".csv");
	//FDTD_three_dimen_with_PML<half, half>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, 26., sigma_y, "Kahan", "Energy_hh_sigma_26_Kahan.csv", ".csv", ".csv");
	//FDTD_three_dimen_with_PML<half, half>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, 30., sigma_y, "Kahan", "Energy_hh_sigma_30_Kahan.csv", ".csv", ".csv");
	//FDTD_three_dimen_with_PML<half, half>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, 35., sigma_y, "Kahan", "Energy_hh_sigma_35_Kahan.csv", ".csv", ".csv");
	//FDTD_three_dimen_with_PML<half, half>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, 40., sigma_y, "Kahan", "Energy_hh_sigma_40_Kahan.csv", ".csv", ".csv");
	//FDTD_three_dimen_with_PML<half, half>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, 47., sigma_y, "Kahan", "Energy_hh_sigma_47_Kahan.csv", ".csv", ".csv");
	//FDTD_three_dimen_with_PML<half, half>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, 50., sigma_y, "Kahan", "Energy_hh_sigma_50_Kahan.csv", ".csv", ".csv");
	//FDTD_three_dimen_with_PML<half, half>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, 55., sigma_y, "Kahan", "Energy_hh_sigma_55_Kahan.csv", ".csv", ".csv");
	//FDTD_three_dimen_with_PML<half, half>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, 60., sigma_y, "Kahan", "Energy_hh_sigma_60_Kahan.csv", ".csv", ".csv");



	  /*
	double sigma = 20.0;
	 
	ofstream numb_energy("Граф коэф отражения в за-ти от сигма для half-half Kahan.csv");
	while (sigma <= 60.0)
	{
		numb_energy << sigma << ";" << FDTD_three_dimen_with_PML<half, half>(128, 16, 16, 8.0_h * 3.1415_h, 0.005_h, delta, sigma, sigma_y, "Kahan", ".csv", ".csv", ".csv") << endl;
		sigma += 1.;
	}
	
	 
	numb_energy << endl << endl;
	numb_energy.close();
*/


	// Check_Curant(2.0 * M_PI / (double)64, 2.0 * M_PI / (double)64, 2.0 * M_PI / (double)64, 0.005);
	//FDTD_three_dimen_print_to_file_task1<float>(64, 64, 64, 8.0 * M_PI, 0.005, "No_Kahan", "Error_without_PML_float_8pi.csv");
	//FDTD_three_dimen_print_to_file_task1<half>(64, 64, 64, 8.0_h * 3.1415_h, 0.005_h, "No_Kahan", "Error_without_PML_half_8pi_No_Kahan.csv");
	//FDTD_three_dimen_print_to_file_task1<half>(64, 64, 64, 8.0_h * 3.1415_h, 0.005_h, "Kahan", "Error_without_PML_half_8pi_Kahan.csv");

	return 0;
}
