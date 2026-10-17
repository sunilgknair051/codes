#include <iostream> 
#include <fstream> 
#include <vector>  
#include <complex> 
#include <math.h>
#include "core3d.h" 

using namespace std;

void Core3d::Count_and_Allocate(int _nx, int _ny, int _nz, int _ndir, int _ntforOut) {
	dy.uxTop.resize(_nx);
	dy.dr.resize(_nx);
	dy.ux.resize(_nx);
	dy.uy.resize(_nx);
	dy.uz.resize(_nx);
	dy.h1.resize(_nx);
	dy.h2.resize(_nx);
	for (int x = 0; x < _nx; x++) {
		dy.uxTop[x].resize(_nz);
		dy.dr[x].resize(_ny);
		dy.ux[x].resize(_ny);
		dy.uy[x].resize(_ny);
		dy.uz[x].resize(_ny);
		dy.h1[x].resize(_ny);
		dy.h2[x].resize(_ny);
		for (int y = 0; y < _ny; y++) {
			dy.dr[x][y].resize(_nz);
			dy.ux[x][y].resize(_nz);
			dy.uy[x][y].resize(_nz);
			dy.uz[x][y].resize(_nz);
			dy.h1[x][y].resize(_nz);
			dy.h2[x][y].resize(_nz);
			for (int z = 0; z < _nz; z++) {
				dy.h1[x][y][z].resize(_ndir);
				dy.h2[x][y][z].resize(_ndir);
			}
		}
	}
	dy.y_arr.resize(_ny);
	dy.y_SI_arr.resize(_ny);
	dy.ux_theo.resize(_ny);
}
void Core3d::freeMem() {
	dy.y_arr.resize(0);
	dy.y_SI_arr.resize(0);

	for (int x = 0; x < cfg.nx; x++) {
		for (int y = 0; y < cfg.ny; y++) {
			for (int z = 0; z < cfg.nz; z++) {
				dy.h1[x][y][z].resize(0);
				dy.h2[x][y][z].resize(0);
			}
			dy.h1[x][y].resize(0);
			dy.h2[x][y].resize(0);
			dy.dr[x][y].resize(0);
			dy.ux[x][y].resize(0);
			dy.uy[x][y].resize(0);
			dy.uz[x][y].resize(0);
		}
		dy.h1[x].resize(0);
		dy.h2[x].resize(0);
		dy.dr[x].resize(0);
		dy.ux[x].resize(0);
		dy.uy[x].resize(0);
		dy.uz[x].resize(0);
		dy.uxTop[x].resize(0);
	}
	dy.h1.resize(0);
	dy.h2.resize(0);
	dy.dr.resize(0);
	dy.ux.resize(0);
	dy.uy.resize(0);
	dy.uz.resize(0);
	dy.uxTop.resize(0);
	dy.y_arr.resize(0);
	dy.y_SI_arr.resize(0);
	dy.ux_theo.resize(0);
}
complex<double> Core3d::streamLB(complex<double> _Zliq) {
	complex<double> FxLiqOnWall = 0.0;

	for (int x = 0; x < cfg.nx; x++) {
		for (int y = 0; y < cfg.ny; y++) {
			for (int z = 0; z < cfg.nz; z++) {
				for (int i = 0; i < vs.ndir; i++) {
					int ibar = vs.no_slip[i];
					int xmd = (dy.h1.size() + x - vs.dirx[i]) % static_cast<int>(dy.h1.size());
					int ymd = y - vs.diry[i];
					int zmd = (dy.h1.size() + z - vs.dirz[i]) % static_cast<int>(dy.h1.size());
					if ((ymd < cfg.ny) && (ymd >= 0)) {
						dy.h2[x][y][z][i] = dy.h1[xmd][ymd][zmd][i];
						continue;
					}
					if (ymd >= cfg.ny) {
					//	dy.uxTop[x][z] = -2.0 *vs.wi[7] * (dy.h1[x][y][z][7] - dy.h1[x][y][z][14]) / (_Zliq - 1.0 / 3.0);
						dy.uxTop[x][z] = 2.0 * (dy.h1[x][y][z][7] - dy.h1[x][y][z][14]) / (_Zliq - 1.0 / 3.0);
						dy.h2[x][y][z][i] = dy.h1[x][y][z][ibar] - 6.0 * vs.wi[ibar] * vs.dirx[ibar] * dy.uxTop[x][z];
						continue;
					}
					if (ymd < 0) {
						dy.h2[x][y][z][i] = dy.h1[x][y][z][ibar] - 6.0 * vs.wi[ibar] * vs.dirx[ibar];
						FxLiqOnWall += (dy.h2[x][y][z][i] + dy.h1[x][y][z][ibar]) * static_cast<double>(vs.dirx[i]);
						continue;
					}
				}
			}
		}
	}
	return FxLiqOnWall;
}
void Core3d::compute_drho_u() {
	for (size_t x = 0; x < dy.h1.size(); x++) {
		for (size_t y = 0; y < dy.h1[x].size(); y++) {
			for (size_t z = 0; z < dy.h1[x][y].size(); z++) {
				dy.dr[x][y][z] = 0.0;
				dy.ux[x][y][z] = 0.0;
				dy.uy[x][y][z] = 0.0;
				dy.uz[x][y][z] = 0.0;
				for (size_t i = 0; i < dy.h1[x][y][z].size(); i++) {
					dy.dr[x][y][z] += dy.h2[x][y][z][i];
					dy.ux[x][y][z] += dy.h2[x][y][z][i] * static_cast<double>(vs.dirx[i]);
					dy.uy[x][y][z] += dy.h2[x][y][z][i] * static_cast<double>(vs.diry[i]);
					dy.uy[x][y][z] += dy.h2[x][y][z][i] * static_cast<double>(vs.dirz[i]);
				}
			}
		}			 
	}
}
void Core3d::collide(complex<double> _tauinv, double _omega, double _omega2by2) {
	for (size_t x = 0; x < dy.h1.size(); x++) {
		for (size_t y = 0; y < dy.h1[x].size(); y++) {
			for (size_t z = 0; z < dy.h1[x][y].size(); z++) {
				for (size_t i = 0; i < dy.h1[x][y][z].size(); i++) {
					complex<double>cidotu = static_cast<double>(vs.dirx[i])* dy.ux[x][y][z] + static_cast<double>(vs.diry[i])* dy.uy[x][y][z];
					complex<double>heq = vs.wi[i] * (dy.dr[x][y][z] + 3.0 * cidotu);
					dy.h1[x][y][z][i] = dy.h2[x][y][z][i] - _tauinv * (dy.h2[x][y][z][i] - heq) - prm.IU * _omega * dy.h2[x][y][z][i] - _omega2by2 * heq;
				}
			}
		}
	}
}
complex<double> Core3d::calc_dfc_by_n(double _omega, complex<double> _FxLiqOnWall) {
	double dq = prm.dq_SI / cfg.gridreso_SI;
	double f0 = _omega / (2.0 * M_PI);
	double Zq = 2.0 * dq * f0 * prm.rhoq;
	complex<double> dfc = prm.f0_SI * prm.IU / M_PI * _FxLiqOnWall / Zq;
	return dfc;
}
void Core3d::run() {
	complex<double> nu = 1.0 / 6.0 * cfg.nufac * (1.0 - prm.IU * cfg.tandelinv);
	complex<double> tau = 3.0 * nu + 0.5;
	double delta_SI = pow((2 * prm.nu_SI / (2 * M_PI * prm.f0_SI)), 0.5);
	double delta = delta_SI / cfg.gridreso_SI;
	double omega = 2.0 * real(nu) / (delta * delta);
	double omega2by2 = 1.0 / 2.0 * omega * omega;
	complex<double> tauinv = 2.0 / (6.0 * nu + 1.0);
	complex<double> Zliq = pow((prm.IU * prm.rho0 * omega * prm.rho0 * nu), 0.5);
	complex<double> c_theo = pow((prm.IU * omega * nu), 0.5);
	complex<double> k_theo = omega / c_theo;
	setConfig();
	int nsteps = static_cast<int>(cfg.ncycles / omega);
	cout << "nsteps = " << nsteps << endl;
	int outInterval = 100;
	int ntforOut = int(nsteps / outInterval);
	Count_and_Allocate(cfg.nx, cfg.ny, cfg.nz, vs.ndir, ntforOut);
	dy.ux_theo.resize(dy.y_arr.size());
	for (int y = 0; y < dy.y_arr.size(); y++) {
		dy.y_arr[y] = y + 1.0 / 2.0;
		dy.y_SI_arr[y] = dy.y_arr[y] * cfg.cellheight_SI / dy.y_arr.size() * 1.0e9;
		dy.ux_theo[y] = exp(-prm.IU * k_theo * dy.y_arr[y]);
	}
	dy.init_h();
	complex<double> FxLiqOnWall;
	
	for (int n = 0; n < nsteps; n++) {
		FxLiqOnWall = streamLB(Zliq);
		compute_drho_u();
		collide(tauinv, omega, omega2by2);
		--outInterval;
		if (outInterval == 0) {
			cout << "step: " << n << FxLiqOnWall << endl;
			complex<double> dfc = calc_dfc_by_n(omega, FxLiqOnWall);
			dfc_by_n.push_back(dfc);
			outInterval = 300;
		}
	}
	double dq = prm.dq_SI / cfg.gridreso_SI;
	double f0 = omega / (2.0 * M_PI);
	double Zq = 2.0 * dq * f0 * prm.rhoq;
	complex<double> dfc = prm.f0_SI * prm.IU / M_PI * FxLiqOnWall / Zq;
	cout << "dfc: " << dfc << endl;
	complex<double> dfc_theo = prm.f0_SI * prm.IU / M_PI * Zliq / Zq;
	cout << "dfc_theo: " << dfc_theo << endl;
	cout << "frac_error_df " << (dfc.real() - dfc_theo.real()) / dfc_theo.real() * 100 << endl;
	cout << "frac_error_dG " << (dfc.imag() - dfc_theo.imag()) / dfc_theo.imag() * 100 << endl;
    cout << "grid res :" << cfg.gridreso_SI << endl;
	writeResults("res");
	freeMem();
}

void Core3d::setConfig() {
	cfg.target_acc_Hz = 1;
	cout << "target_acc_Hz: " << cfg.target_acc_Hz << endl;
	cfg.rad_Sph_SI = 0;
	cfg.cellheight_SI; //= 2000e-9;
	cfg.gridreso_SI; //= 10.0e-9;
	cfg.cell_diam_by_sphere_diam = 2.;
	cfg.truncation_by_diameter = 0.5;
	cfg.nufac = 1.;
	cfg.tandelinv = 1.;
	cfg.iovt = 1;
	cfg.ncycles ;
}

void Core3d::writeResults(string _filename) {
	ofstream outfile;
	outfile.open(_filename + "_wave.dat");
	//for (size_t y = 0; y < cfg.ny;y++){
    for (size_t y = 0; y < dy.y_SI_arr.size(); y++) {
		outfile << dy.y_SI_arr[y] 
			<< "\t" << dy.ux[0][y][0].real()
			<< "\t" << dy.ux[0][y][0].imag() 
			<< "\t" << dy.ux_theo[y].real() 
			<< "\t" << dy.ux_theo[y].imag() << endl;
	}
	outfile.close();
	outfile.open(_filename + "_ring_in.dat");
    int count = 0;
	for (auto& elem : dfc_by_n) {
		outfile << count << "\t" << elem.real() << "\t" << elem.imag() << endl;
		++count;
	}
	outfile.close();
	
}