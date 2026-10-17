//
// JAN23
//
// Created on: 24.1.2020

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include "core2d.h"

using namespace std;
Core2d::Core2d(string _filename) {
	readConfig(_filename);
}

void Core2d::Count_and_Allocate(int _nx, int _ny, int _ndir, int _ntforOut) {
	dy.uxTop.resize(_nx);
	dy.dr.resize(_nx);
	dy.ux.resize(_nx);
	dy.uy.resize(_nx);
	dy.h1.resize(_nx);
	dy.h2.resize(_nx);
	for (int x = 0; x < _nx; x++) {
		dy.dr[x].resize(_ny);
		dy.ux[x].resize(_ny);
		dy.uy[x].resize(_ny);
		dy.h1[x].resize(_ny);
		dy.h2[x].resize(_ny);
		for (int y = 0; y < _ny; y++) {
			dy.h1[x][y].resize(_ndir);
			dy.h2[x][y].resize(_ndir);
		}
	}
	dy.nx = _nx;
	dy.ny = _ny;
	dy.ndir = _ndir;
	dy.y_arr.resize(_ny);
	dy.y_SI_arr.resize(_ny);
}
void Core2d::freeMem() {
	dy.uxTop.resize(0);
	dy.dr.resize(0);
	dy.ux.resize(0);
	dy.uy.resize(0);
	for (int x = 0; x < dy.nx; x++) {
		for (int y = 0; y < dy.ny; y++) {
			dy.h1[x][y].resize(0);
			dy.h2[x][y].resize(0);
		}
		dy.h1[x].resize(0);
		dy.h2[x].resize(0);
	}
	dy.h1.resize(0);
	dy.h2.resize(0);

	dy.y_arr.resize(0);
	dy.y_SI_arr.resize(0);
}
complex<double> Core2d::streamLB(complex<double> _Zliq) {
	complex<double> FxLiqOnWall = 0.0;

	for (int x = 0; x < dy.nx; x++) {
		for (int y = 0; y < dy.ny; y++) {
			for (int i = 0; i < dy.ndir; i++) {
				int ibar = vs.no_slip[i];
				int xmd = (dy.h1.size() + x - vs.dirx[i]) % static_cast<int>(dy.h1.size());
				int ymd = y - vs.diry[i];
				if ((ymd < dy.ny) && (ymd >= 0)) {
					dy.h2[x][y][i] = dy.h1[xmd][ymd][i];
					continue;
				}
				if (ymd >= dy.ny) {
					complex<double> uxTop = 2.0 * (dy.h1[x][y][5] - dy.h1[x][y][6]) / (_Zliq + 1.0 / 3.0);
					dy.h2[x][y][i] = dy.h1[x][y][ibar] - 6.0 * vs.wi[ibar] * vs.dirx[ibar] * uxTop;
					continue;
				}
				if (ymd < 0) {
					dy.h2[x][y][i] = dy.h1[x][y][ibar] - 6.0 * vs.wi[ibar] * vs.dirx[ibar];
					FxLiqOnWall += (dy.h2[x][y][i] + dy.h1[x][y][ibar]) * static_cast<double>(vs.dirx[i]);
					continue;
				}
			}
		}
	}
	return FxLiqOnWall;
}
void Core2d::compute_drho_u() {
	for (size_t x = 0; x < dy.h1.size(); x++) {
		for (size_t y = 0; y < dy.h1[x].size(); y++) {
			dy.dr[x][y] = 0.0;
			dy.ux[x][y] = 0.0;
			dy.uy[x][y] = 0.0;
			for (size_t i = 0; i < dy.h1[x][y].size(); i++) {
				dy.dr[x][y] += dy.h2[x][y][i];
				dy.ux[x][y] += dy.h2[x][y][i] * static_cast<double>(vs.dirx[i]);
				dy.uy[x][y] += dy.h2[x][y][i] * static_cast<double>(vs.diry[i]);
			}
		}
	}
}
void Core2d::collide(complex<double> _tauinv, double _omega, double _omega2by2) {
	for (size_t x = 0; x < dy.h1.size(); x++) {
		for (size_t y = 0; y < dy.h1[x].size(); y++) {
			for (size_t i = 0; i < dy.h1[x][y].size(); i++) {
				complex<double>cidotu = static_cast<double>(vs.dirx[i])* dy.ux[x][y] + static_cast<double>(vs.diry[i])* dy.uy[x][y];
				complex<double>heq = vs.wi[i] * (dy.dr[x][y] + 3.0 * cidotu);
				dy.h1[x][y][i] = dy.h2[x][y][i] - _tauinv * (dy.h2[x][y][i] - heq) - prm.IU * _omega * dy.h2[x][y][i] - _omega2by2 * heq;
			}
		}
	}
}
complex<double> Core2d::calc_dfc_by_n(double _omega, complex<double> _FxLiqOnWall) {
	double dq = prm.dq_SI / cfg.gridreso_SI;
	double f0 = _omega / (2.0 * M_PI);
	double Zq = 2.0 * dq * f0 * prm.rhoq;
	complex<double> dfc = prm.f0_SI * prm.IU / M_PI * _FxLiqOnWall / Zq;
	return dfc;
}
void Core2d::run(int _ncycles) {
	cout << "Start" << endl;
	cfg.ncycles = _ncycles;
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
	vector<complex<double>> ux_theo;
	int nsteps = static_cast<int>(cfg.ncycles / omega);
	cout << "nsteps = " << nsteps << endl;
	cout << "ncycles = " << _ncycles << endl;
	int outInterval = 100;

	int ntforOut = int(nsteps / outInterval);
	Count_and_Allocate(1, 80, vs.ndir, ntforOut);

	ux_theo.resize(dy.y_arr.size());
	for (int y = 0; y < dy.y_arr.size(); y++) {
		dy.y_arr[y] = y + 1.0 / 2.0;
		dy.y_SI_arr[y] = dy.y_arr[y] * cfg.cellheight_SI / dy.y_arr.size() * 1.0e9;
		ux_theo.at(y) = exp(-prm.IU * k_theo * dy.y_arr[y]);
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
			outInterval = 100;
		}
	}
	double dq = prm.dq_SI / cfg.gridreso_SI;
	double f0 = omega / (2.0 * M_PI);
	double Zq = 2.0 * dq * f0 * prm.rhoq;
	complex<double> dfc = prm.f0_SI * prm.IU / M_PI * FxLiqOnWall / Zq;
	cout << "dfc: " << dfc << endl;
	complex<double> dfc_theo = prm.f0_SI * prm.IU / M_PI * Zliq / Zq;
	cout << "dfc_theo: " << dfc_theo << endl;
	freeMem();
}
void Core2d::readConfig(string _filename) {
	ifstream is(_filename, ios::in);
	if (!is.is_open()) {
		cout << "file " << _filename << " not found" << endl;
		return;
	}
	string line;

	int count = 0;
	size_t countEnd = 0;

	while (getline(is, line)) {
		size_t pos;
		vector<string> splitstring;
		pos = line.find_first_of(" \t:;,");
		while (pos != std::string::npos)
		{
			line[pos] = ' ';
			pos = line.find_first_of(" \t:;,", pos + 1);
		}

		bool sepmode = false;
		string substr = "";
		for (int i = 0; i < line.length(); i++) {
			if (line.at(i) == ' ') {
				if (!sepmode) {
					sepmode = true;
					//close splitstring
					splitstring.push_back(substr);
					substr = "";
				}
			}
			else {
				substr.append(&line.at(i), 1);
				sepmode = false;
			}
		}
		if (!sepmode) {
			if (substr != "") {
				splitstring.push_back(substr);
			}
		}
		++count;
		switch (count) {
		case 1: {
			cfg.target_acc_Hz = stod(splitstring.at(1));
			break;
		}
		case 2: {
			cfg.rad_Sph_SI = stod(splitstring.at(1));
			break;
		}
		case 3: {
			cfg.cellheight_SI = stod(splitstring.at(1));
			break;
		}
		case 4: {
			cfg.gridreso_SI = stod(splitstring.at(1));
			break;
		}
		case 5: {
			cfg.cell_diam_by_sphere_diam = stod(splitstring.at(1));
			break;
		}
		case 6: {
			cfg.truncation_by_diameter = stod(splitstring.at(1));
			break;
		}
		case 7: {
			cfg.nufac = stod(splitstring.at(1));
			break;
		}
		case 8: {
			cfg.tandelinv = stod(splitstring.at(1));
			break;
		}
		case 9: {
			cfg.iovt = stoi(splitstring.at(1));
			break;
		}
		}
	}
}
void Core2d::writeResults(string _filename) {
	ofstream outfile;
	outfile.open(_filename + "_wave.dat");
	//cout << "y\treal(ux)\timag(ux)" << endl;
	for (size_t y = 0; y < dy.y_SI_arr.size(); y++) {
		outfile << dy.y_SI_arr[y] << "\t" << dy.ux[0][y].real() << "\t" << dy.ux[0][y].imag() << endl;
	}
	outfile.close();
	outfile.open(_filename + "_ring_in.dat");
	//cout << "t\tdf\tdG" << endl;
	int count = 0;
	for (auto& elem : dfc_by_n) {
		outfile << count << "\t" << elem.real() << "\t" << elem.imag() << endl;
		++count;
	}
	outfile.close();
}