//06.07.2020
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
		const complex<double> IU = complex<double>(0, 1);
		string filename;
		double rad_Sph_SI, cellheight_SI, gridreso_SI, cell_diam_by_sphere_diam, frac_above_surf;
		double nu, tauinv, rho0;
		double nufac, tandelinv;
		complex<double> nu_Sph, tauinv_Sph, Zliq, Zliq_Sph;
		complex<double> one_m_tauinv_m_Iom, tauinv_m_om2b2;
		complex<double> one_m_tauinv_m_Iom_Sph, tauinv_m_om2b2_Sph;
		double target_acc_Hz, char_time_fac, ovtorder;
		int iovt, jqmax, nx, ny, nz;
		double dq, f0, Zq, rhoq, delta, omega, omega2by2;
		double dq_SI, f0_SI, Zq_SI, nu_SI, delta_SI;
		complex<double> c_theo, k_theo;
		int outInterval;
		double t_computation_secs;
		double SphR, Sphx, Sphy, Sphz;
		bool B_Advanced_Top, B_Reference;
	};
	

	struct allocMem {
		valarray<valarray<valarray<valarray<complex<double>>>>> H_to_uT;
		valarray<valarray<valarray<valarray<complex<double>>>>> h1;
		valarray<valarray<valarray<valarray<complex<double>>>>> h2;
		valarray<valarray<valarray<complex<double>>>> ux;
		valarray<valarray<valarray<complex<double>>>> uy;
		valarray<valarray<valarray<complex<double>>>> uz;
		valarray<valarray<valarray<complex<double>>>> dr;
		valarray<valarray<complex<double>>> ZLiq_Top;
		valarray<valarray<complex<double>>> uxT,uyT,uzT;
		valarray<valarray<valarray<complex<double>>>> one_m_tauinv_m_Iom_arr;
		valarray<valarray<valarray<complex<double>>>> one_m_tauinv_p_m_Iom_arr;
		valarray<valarray<valarray<complex<double>>>> one_m_tauinv_m_m_Iom_arr;
		valarray<valarray<valarray<complex<double>>>> tauinv_m_om2_by_2_arr;
		valarray<valarray<valarray<complex<double>>>> tauinv_p_m_om2_by_2_arr;
		valarray<valarray<valarray<complex<double>>>> tauinv_m_m_om2_by_2_arr;
	};

	void compute_Matrices();
	void Set_Config();
	void compute_uT();
	void Calc_Parameters_from_Config();
	void Write_Header_Info(string _i_sim_str);
	void Write_Header_Filename();
	void Read_Header_Filename();
	void WriteDisplacementField(string _i_sim_str);
	void WritePopulations(string _i_sim_str);
	void Set_VE_Pars();
	void Count_and_Allocate();
	void Free_Mem();
	complex<double> Stream_LB();
	void Compute_drho_u();
	void Collide();
	void Simulate(string _i_sim_str);
	void init_h();
	void WriteDisplacementField_Sphere(string _i_sim_str);
	D3Q19 vst;
	config cfg;
	allocMem dyn;
};