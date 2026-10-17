#pragma once
#include <string>
#include <vector>
#include <map>
#include <valarray>
#include <complex>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <io.h>

using namespace std;
class Core3d {
public:
	void run();
private:
	struct D3Q19 {
		const unsigned int ndir = 19;		// lattice weights and lengths
		const double zeroWeight = 1.0 / 3.0;
		const double adjaWeight = 1.0 / 18.0;
		const double diagWeight = 1.0 / 36.0;
		const int dirx[19] = {0,1,-1,0, 0,0, 0,1,-1,1,-1,0, 0, 1,-1, 1,-1, 0, 0};
		const int diry[19] = {0,0,0, 1,-1,0, 0,1,-1,0, 0,1,-1,-1, 1, 0, 0, 1,-1};
		const int dirz[19] = {0,0,0, 0, 0,1,-1,0, 0,1,-1,1,-1, 0, 0,-1, 1,-1, 1}; // new dir to prevent static pop up
		const double dirxD[19] = { 0,1,-1,0, 0,0, 0,1,-1,1,-1,0, 0, 1,-1, 1,-1, 0, 0 };//
		const double diryD[19] = { 0,0,0, 1,-1,0, 0,1,-1,0, 0,1,-1,-1, 1, 0, 0, 1,-1 };//
		const double dirzD[19] = { 0,0,0, 0, 0,1,-1,0, 0,1,-1,1,-1, 0, 0,-1, 1,-1, 1 };//
		const double wi[19] = { zeroWeight,
			adjaWeight,adjaWeight,adjaWeight,adjaWeight,adjaWeight,adjaWeight,
			diagWeight,diagWeight,diagWeight,diagWeight,diagWeight,diagWeight, 
			diagWeight,diagWeight,diagWeight,diagWeight,diagWeight,diagWeight };
		const int no_slip[19] = { 0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17 };
	};

	struct config {
		const complex<double>IU = complex<double>(0, 1);
		string filename;
		double f0_SI, Zq_SI, dq_SI, rhoq, nu_SI, delta_f0_SI;
		double nu0, rho0;
		double target_acc_Hz;
		double rad_Sph_SI, cellheight_SI, gridreso_SI, cell_diam_by_sphere_diam, frac_above_surf;
		double nufac, tandelinv;
		int ncycles_max;
		int iovt;
		double ovtorder;
		int nx, ny, nz;
		double dq, f0, Zq;
		complex<double> nu;
		complex<double> tauinv,   one_m_tauinv_m_Iom,   tauinv_m_om2_by_2;
		complex<double> tauinv_0, one_m_tauinv_m_Iom_0, tauinv_m_om2_by_2_0;

		complex<double> tauinv_p,   one_m_tauinv_p_m_Iom, tauinv_p_m_om2_by_2;
		complex<double> tauinv_p_0, one_m_tauinv_p_m_Iom_0, tauinv_p_m_om2_by_2_0;

		complex<double> tauinv_m,   one_m_tauinv_m_m_Iom, tauinv_m_m_om2_by_2;
		complex<double> tauinv_m_0, one_m_tauinv_m_m_Iom_0, tauinv_m_m_om2_by_2_0;

		double delta_SI, delta, omega, omega2by2;
		complex<double> Zliq, Zliq0, c_theo, k_theo;
		int nsteps_max, outInterval;
		double t_computation_secs;
		double SphR, Sphx, Sphy, Sphz;
		bool B_TRT;
		double Lambda;
	};

	struct allocMem {
		valarray<valarray<valarray<valarray<complex<double>>>>> h1;
		valarray<valarray<valarray<valarray<complex<double>>>>> h2;
		valarray<valarray<valarray<complex<double>>>> ux;
		valarray<valarray<valarray<complex<double>>>> uy;
		valarray<valarray<valarray<complex<double>>>> uz;
		valarray<valarray<valarray<complex<double>>>> dr;
		valarray<valarray<complex<double>>> ZLiq_Top;
		valarray<valarray<valarray<complex<double>>>> one_m_tauinv_m_Iom_arr;
		valarray<valarray<valarray<complex<double>>>> one_m_tauinv_p_m_Iom_arr;
		valarray<valarray<valarray<complex<double>>>> one_m_tauinv_m_m_Iom_arr;
		valarray<valarray<valarray<complex<double>>>> tauinv_m_om2_by_2_arr;
		valarray<valarray<valarray<complex<double>>>> tauinv_p_m_om2_by_2_arr;
		valarray<valarray<valarray<complex<double>>>> tauinv_m_m_om2_by_2_arr;

		valarray<double> h1_lo_z_real, h1_lo_z_imag, h1_hi_z_real, h1_hi_z_imag;

		valarray<double> y_arr;
		valarray<double> y_SI_arr;
		valarray<complex<double>> ux_theo;
		
	/*	void init_h() 
		{
			for (size_t x = 0; x < h1.size(); x++) {
				for (size_t y = 0; y < h1[x].size(); y++) {
					for (size_t z = 0; z < h1[x][y].size(); z++) {
						for (size_t i = 0; i < h1[x][y][z].size(); i++) {
							h1[x][y][z][i] = 0;
						}
					}
				}
			}*/
	};
	
	void Set_Config();
	void Calc_Parameters_from_Config();
	void Write_Header_Info(string _i_sim_str);
	void Write_Header_Filename();
	void Read_Header_Filename();
	void WriteDisplacementField(string _i_sim_str);
	void Set_VE_Pars();
	void Count_and_Allocate();
	void Free_Mem();
	complex<double> Stream_LB();
	void Export_to_arrays_for_MPI();
	void Import_from_arrays_for_MPI();
	void Compute_drho_u();
	void Collide();
	void Calculate_h1_from_TRT(int x, int y, int z);
	void Simulate(string _i_sim_str);
	void init_h_as_wave();
	void WriteDisplacementField_Sphere(string _i_sim_str);
	D3Q19 vst;
	config cfg;
	allocMem dyn;
};