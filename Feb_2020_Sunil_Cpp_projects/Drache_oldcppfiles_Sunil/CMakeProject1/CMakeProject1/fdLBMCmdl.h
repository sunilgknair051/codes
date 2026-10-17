#pragma once

#include <string>
using namespace std;
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>

struct LBParameter {
	string modelFile = "";
	string velocitySet = "D2Q9";
	bool Lattice3D = false;
	unsigned int steps = 0;
	// base: X = Y = (Z)
	int baseDimension = 0;
	int dimensionX = 0;
	int dimensionY = 0;
	int dimensionZ = 0;
	bool periodicX = true;
	bool periodicY = true;
	bool periodicZ = true;
	double x_sphere = 0.0;
	double y_sphere = 5.0;
	double z_sphere = 0.0;
	string bounceBackZ = "";
	double fluidDensity = 1.0;
	void init() {
		string vSets = "D3Q15 D3Q19 D3Q27";
		if (vSets.compare(velocitySet) != 0) {
			// 3D velocitySet identified
			Lattice3D = true;
		}
		if (dimensionX == 0) {
			// dimension not set
			dimensionX = baseDimension;
		}
		if (dimensionY == 0) {
			// dimension not set
			dimensionY = baseDimension;
		}
		if ((dimensionZ == 0) && Lattice3D) {
			// 3D model
			// dimension not set
			dimensionZ = baseDimension;
		}
	}
	size_t getDimension() {
		size_t dimension = dimensionX * dimensionY;
		if (Lattice3D) {
			dimension *= dimensionZ;
		}
		return dimension;
	}
};
struct modelParameter {
	const double nu = 1.0 / 6.0;
	const double tauinv = 2.0 / (6.0 * nu + 1.0);
	const double omtauinv = 1.0 - tauinv;
	const double delta = 10.0;
	const double omega = 2.0 * nu / (delta * delta);
	const double omega2by2 = omega * omega / 2.0;
	const double uxWall = 1.0;
	const double rho0 = 1.0;
	const double rho_Liq_SI = 1000.0;
	const double nu_SI = 1.0e-6;
	const double f0_SI = 5.0e6;
	const double delta_f0_SI = sqrt(2.0 * nu_SI / (2.0 * M_PI * f0_SI));
	const double delta_f0 = delta;
	const double Zq_SI = 8.8e6;
};
struct Fx {
	double LiqOnWall;
	double LiqOnSphe;
};
bool filenameValid(const string& _filename);
void getParameterFromModelFile(LBParameter);