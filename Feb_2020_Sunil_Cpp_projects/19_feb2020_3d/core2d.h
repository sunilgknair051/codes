//
// FEB 26
//
// Created on: 26.2.2020
#pragma once
#include <string>
#include <vector>
#include <map>
#include <valarray>
#include <complex>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;
class Core3d {
public:
	Core3d(string _filename);
	void run(int _ncycles);
	void writeResults(string _filename);
private:
	struct D3Q19 {
		const unsigned int ndir = 19;
		// lattice weights and length
		const double zeroWeight = 1.0 / 3.0;
		const double adjacentWeight = 1.0 / 18.0;
		const double adjacentLength = 1.0;
		const double diagonalWeight = 1.0 / 36.0;
		const double diagonalLength = sqrt(2.0);
		const int dirx[19] = {0,1,-1,0,0,0,0,1,-1,1,-1,0,0,1,-1,1,-1,0,0};
		const int diry[19] = {0,0,0,1,-1,0,0,1,-1,0,0,1,-1,-1,1,0,0,1,-1};
		const int dirz[19] = {0,0,0,0,0,1,-1,0,0,1,-1,1,-1,0,0,-1,1,-1,1};
		const double wi[19] = { zeroWeight,adjacentWeight,adjacentWeight,adjacentWeight,adjacentWeight,adjacentWeight,adjacentWeight\
			diagonalWeight,diagonalWeight,diagonalWeight,diagonalWeight, diagonalWeight, diagonalWeight, diagonalWeight, diagonalWeight\
			diagonalWeight, diagonalWeight , diagonalWeight, diagonalWeight };
		const double length[19] = { 0.0, adjacentLength,adjacentLength,adjacentLength,adjacentLength,adjacentLength,adjacentLength\
			diagonalLength,diagonalLength,diagonalLength,diagonalLength, diagonalLength, diagonalLength, diagonalLength, diagonalLength\
			diagonalLength , diagonalLength , diagonalLength , diagonalLength };
		const int no_slip[19] = { 0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17 };
	};

	struct om_dependent_variables {
		complex<double> omtauinv_m_Iom = complex<double>(0, 0);
		complex<double> omtauinv_Sph_m_Iom = complex<double>(0, 0);
		complex<double> Zliq = complex<double>(0, 0);
		complex<double> tauinv_m_omega2by2 = complex<double>(0, 0);
		complex<double> tauinv_Sph_m_omega2by2 = complex<double>(0, 0);
		int nsteps_mx = 0;
	};

	struct parameter {
		const complex<double>IU = complex<double>(0, 1);
		const double f0_SI = 5.0e6;
		const double nu_SI = 1.0e-6;
		double delta_f0_SI = sqrt((2.0 * nu_SI / (2.0 * M_PI * f0_SI)));
		const double Zq_SI = 8.8e6;
		const double rho0_SI = 1e3;
		const double df_analy = 716.799;
		double nu = 1.0 / 6.0;
		double tau = 3.0 * nu + 0.5;
		const double rho0 = 1.0;
		double tauinv = 2.0 / (6.0 * nu + 1.0);
		double omtauinv = 1.0 - tauinv;
		const int n_Sph_max = 10;
		const double dq_SI = 332.075e-6;
		const double rhoq = 2.65;
		double* x_Sph = nullptr;
		double* y_Sph = nullptr;
		double* z_Sph = nullptr;
		double* xfrac_Sph = nullptr;
		double* zfrac_Sph = nullptr;
		double rad_Sph = 0.0;
		double delta_f0 = 0.0;
		bool B_Save_uxyz = false;
		complex<double> nuc_Sph = complex<double>(0.0, 0.0);
		complex<double> tauinv_Sph = complex<double>(0.0, 0.0);
		complex<double> omtauinv_Sph = complex<double>(0.0, 0.0);
		// ??! als double deklariert, aber eher als int verwendet, L108
		const int ncycles_mx = 2;
		const int nslices = 3;
		om_dependent_variables omDependent;
	};

	struct config {
		int ncycles = 20;
		double target_acc_Hz = 0.0;
		double rad_Sph_SI = 0.0;
		double cellheight_SI = 2000.0e-9;
		double gridreso_SI = 25.0e-9;
		double cell_diam_by_sphere_diam = 0.0;
		double truncation_by_diameter = 0.0;
		double nufac = 1.0;
		double tandelinv = 1.0;
		double nufac_Re_Sph = 0.0;
		double nufac_Im_Sph = 0.0;
		double Lambda = 0.0;
		int iovt = 0;
		size_t n_Sph = 0;
		vector<double> xfrac_Sph;
		vector<double> zfrac_Sph;
	};

	struct allocMem {
		valarray<valarray<valarray<complex<double>>>> h1;
		valarray<valarray<valarray<complex<double>>>> h2;
		valarray<complex<double>> uxTop;
		valarray<valarray<complex<double>>> ux;
		valarray<valarray<complex<double>>> uy;
		valarray<valarray<complex<double>>> uz;
		valarray<valarray<complex<double>>> dr;
		valarray<double> y_arr;
		valarray<double> y_SI_arr;
		void init_h() {
			for (size_t x = 0; x < h1.size(); x++) {
				for (int y = 0; y < h1[x].size(); y++) {
					for (int z = 0; y < h1[x][y].size(); z++) {
						for (int i = 0; i < h1[x][y][z].size(); i++) {
							h1[x][y][z][i] = 0;
						}
					}
				}
			}
		}
		int nx = 0;
		int ny = 0;
		int ndir = 0;
	};

	void Count_and_Allocate(int _nx, int _ny, int _nz, int _ndir, int _ntforOut);
	void freeMem();
	complex<double> streamLB(complex<double> _Zliq);
	void compute_drho_u();
	void collide(complex<double> _tauinv, double _omega, double _omega2by2);
	complex<double> calc_dfc_by_n(double _omega, complex<double> _FxLiqOnWall);
	void readConfig(string _filename);
	D2Q9 vs;
	parameter prm;
	config cfg;
	allocMem dy;
	vector<complex<double>> dfc_by_n;
};