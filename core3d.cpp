#include <iostream> 
#include <fstream> 
#include <vector>  
#include <complex> 
#include <math.h>
#include "core3d.h" 
#include <sstream>
#include <iomanip>
#include <ctime>
#include <cstring>
#include <string>
#include <direct.h>
#include <chrono>

using namespace std;

void Core3d::Calc_Parameters_from_Config() {
	cfg.dq = cfg.dq_SI / cfg.gridreso_SI;
	cfg.ovtorder = 2. * cfg.iovt + 1;
	cfg.delta_SI = pow((2 * cfg.nu_SI / (2 * M_PI * cfg.f0_SI * cfg.ovtorder)), 0.5);
	cfg.ny = static_cast<int>(cfg.cellheight_SI / cfg.gridreso_SI);
	cfg.nu = 1.0 / 6.0;         
	if (cfg.tandelinv == 0) { 

		cfg.nu_Sph = cfg.nu * cfg.nufac;
	}
	else {
		double tandel = 1 / cfg.tandelinv;  //tandel = inv(tandelinv)
		double cos = 1 / pow((1 + tandel * tandel), 0.5);
		double sin = tandel / pow((1 + tandel * tandel), 0.5);
		cfg.nu_Sph =  cfg.nu * cfg.nufac * (sin - cfg.IU * cos);
	}
	if (cfg.B_Reference) {
		cfg.nu_Sph = cfg.nu;
		cfg.nx = 1;
		cfg.nz = 1;
		cfg.B_Advanced_Top = false;
	}                                     //
	else {
		cfg.nu_Sph = cfg.nu * cfg.nufac * (1.0 - cfg.IU * cfg.tandelinv);  //
		cfg.nx = static_cast<int>(2 * cfg.rad_Sph_SI / cfg.gridreso_SI * cfg.cell_diam_by_sphere_diam);
		cfg.nz = cfg.nx;
		cfg.B_Advanced_Top = true;
	}
	cfg.tauinv = 2.0 / (6.0 * cfg.nu + 1.0);            
	cfg.tauinv_Sph = 2.0 / (6.0 * cfg.nu_Sph + 1.0);
	cfg.one_m_tauinv_m_Iom = 1. - cfg.tauinv - cfg.IU * cfg.omega;
	cfg.one_m_tauinv_m_Iom_Sph = 1. - cfg.tauinv_Sph - cfg.IU * cfg.omega;
	cfg.tauinv_m_om2b2 = cfg.tauinv - cfg.omega * cfg.omega / 2.;
	cfg.tauinv_m_om2b2_Sph = cfg.tauinv_Sph - cfg.omega * cfg.omega / 2.;

	cfg.delta = cfg.delta_SI / cfg.gridreso_SI;
	cfg.omega = 2.0 * cfg.nu / (cfg.delta * cfg.delta);
	cfg.omega2by2 = 1.0 / 2.0 * cfg.omega * cfg.omega;
	cfg.Zliq = pow((cfg.IU * cfg.rho0 * cfg.omega * cfg.rho0 * cfg.nu), 0.5);
	cfg.Zliq_Sph = pow((cfg.IU * cfg.rho0 * cfg.omega * cfg.rho0 * cfg.nu_Sph), 0.5);
	cfg.c_theo = pow((cfg.IU * cfg.omega * cfg.nu), 0.5);
	cfg.k_theo = cfg.omega / cfg.c_theo;
	cfg.f0 = cfg.omega / (2.0 * M_PI * cfg.ovtorder);
	cfg.dq = cfg.dq_SI / cfg.gridreso_SI;
	cfg.Zq = 2.0 * cfg.dq * cfg.f0 * cfg.rhoq;
	cfg.SphR = cfg.rad_Sph_SI / cfg.gridreso_SI;
	cfg.Sphx = static_cast<double>(cfg.nx) / 2. - 0.5;
	cfg.Sphz = static_cast<double>(cfg.nz) / 2. - 0.5;
	cfg.Sphy = 2 * cfg.SphR * (cfg.frac_above_surf - 0.5);
}

void Core3d::compute_Matrices() {
	for (int jx = 0; jx < cfg.jqmax + 1; jx++) {
		double qx = 2. * M_PI * static_cast<double>(jx) / cfg.nx;
		for (int jz = 0; jz < cfg.jqmax + 1; jz++) {
			double qz = 2. * M_PI * static_cast<double>(jz) / cfg.nz;
			double phi = atan2(qz, qx);
			double q = pow((qx * qx + qz * qz), 0.5);
			complex<double> Z = cfg.Zliq;
			double om = cfg.omega;
			complex<double> k0 = pow(cfg.IU * om / cfg.nu, 0.5);
			complex<double> xi = pow((q * q + k0 * k0), 0.5) / k0;
			// Hier wird das Koordinaten system gedreht
			double Rot_mat[3][3] = { {cos(phi),0,-sin(phi)},
									 {0,1,0},
									 {sin(phi),0,cos(phi)} };
			double Rot_mat_inv[3][3] = { {cos(phi),0,sin(phi)},
									 {0,1,0},
									 {-sin(phi),0,cos(phi)} };
			complex<double> denom = 3. * k0 * (-3. * cfg.IU * om - 4 * q) * q * q * Z + 18 * pow(q, 4) * Z * Z + \
				pow(k0, 3) * (3. * cfg.IU * om + 4. * q) * xi * (1. + 3. * xi * Z) - 2. * k0 * k0 * q * q * (2. + 9. * xi * xi * Z * Z);

			complex<double> F_mat[3][3] =
			{ { 1. / denom * 3. * k0 * k0 * (-4. * q * q + 3. * cfg.IU * k0 * om * xi + 4. * k0 * q * xi),
			    1. / denom * 9. * k0 * q * (q * q - k0 * k0 * xi * xi) * Z,
				0. },
			{  -1. / denom * 9. * k0 * q * (cfg.IU * k0 * om - 2. * q*q * Z + 2. * k0 * q * xi * Z),
			    1. / denom * 3. * k0 * q * (-q + k0 * xi) * (k0 + 3. * q * Z + 3. * k0 * xi * Z),
				0},
				{ 0,0,1. / (1. / 3. + xi * Z) } };
			complex<double> Mat1[3][3] =
				{{F_mat[0][0] * Rot_mat[0][0] + F_mat[0][1] * Rot_mat[1][0] + F_mat[0][2] * Rot_mat[2][0],
				  F_mat[0][0] * Rot_mat[0][1] + F_mat[0][1] * Rot_mat[1][1] + F_mat[0][2] * Rot_mat[2][1],
				  F_mat[0][0] * Rot_mat[0][2] + F_mat[0][2] * Rot_mat[1][2] + F_mat[0][2] * Rot_mat[2][2] },
				{ F_mat[1][0] * Rot_mat[0][0] + F_mat[1][1] * Rot_mat[1][0] + F_mat[1][2] * Rot_mat[2][0],
				  F_mat[1][0] * Rot_mat[0][1] + F_mat[1][1] * Rot_mat[1][1] + F_mat[1][2] * Rot_mat[2][1],
				  F_mat[1][0] * Rot_mat[0][2] + F_mat[1][2] * Rot_mat[1][2] + F_mat[1][2] * Rot_mat[2][2] },
				{ F_mat[2][0] * Rot_mat[0][0] + F_mat[2][1] * Rot_mat[1][0] + F_mat[2][2] * Rot_mat[2][0],
				  F_mat[2][0] * Rot_mat[0][1] + F_mat[2][1] * Rot_mat[1][1] + F_mat[2][2] * Rot_mat[2][1],
				  F_mat[2][0] * Rot_mat[0][2] + F_mat[2][2] * Rot_mat[1][2] + F_mat[2][2] * Rot_mat[2][2] } };
			dyn.H_to_uT[jx][jz] =
			{ { Rot_mat_inv[0][0] * Mat1[0][0] + Rot_mat_inv[0][1] * Mat1[1][0] + Rot_mat_inv[0][2] * Mat1[2][0],
				Rot_mat_inv[0][0] * Mat1[0][1] + Rot_mat_inv[0][1] * Mat1[1][1] + Rot_mat_inv[0][2] * Mat1[2][1],
				Rot_mat_inv[0][0] * Mat1[0][2] + Rot_mat_inv[0][1] * Mat1[1][2] + Rot_mat_inv[0][2] * Mat1[2][2] },
			  { Rot_mat_inv[1][0] * Mat1[0][0] + Rot_mat_inv[1][1] * Mat1[1][0] + Rot_mat_inv[1][2] * Mat1[2][0],
				Rot_mat_inv[1][0] * Mat1[0][1] + Rot_mat_inv[1][1] * Mat1[1][1] + Rot_mat_inv[1][2] * Mat1[2][1],
				Rot_mat_inv[1][0] * Mat1[0][2] + Rot_mat_inv[1][1] * Mat1[1][2] + Rot_mat_inv[1][2] * Mat1[2][2] },
			  { Rot_mat_inv[2][0] * Mat1[0][0] + Rot_mat_inv[2][1] * Mat1[1][0] + Rot_mat_inv[2][2] * Mat1[2][0],
				Rot_mat_inv[2][0] * Mat1[0][1] + Rot_mat_inv[2][1] * Mat1[1][1] + Rot_mat_inv[2][2] * Mat1[2][1],
				Rot_mat_inv[2][0] * Mat1[0][2] + Rot_mat_inv[2][1] * Mat1[1][2] + Rot_mat_inv[2][2] * Mat1[2][2] } };
		}
	}
};

void Core3d::compute_uT() {
	double x_Sph = (static_cast<double>(cfg.nx) / 2.) - 0.5;
	double z_Sph = (static_cast<double>(cfg.nz) / 2.) - 0.5;
	complex<double> Hxy_coscos;
	complex<double> Hyy_sincos;
	complex<double> Hzy_sinsin;

	for (int x = 0; x < cfg.nx; x++) {
		for (int z = 0; z < cfg.nz; z++) {
			dyn.uxT[x][z] = 0;
			dyn.uyT[x][z] = 0;
			dyn.uzT[x][z] = 0;
		}
	}
	for (int jx = 0; jx < cfg.jqmax + 1; jx++) {
		double qx = 2. * M_PI * static_cast<double>(jx) / cfg.nx;
		for (int jz = 0; jz < cfg.jqmax + 1; jz++) {
			double qz = 2. * M_PI * static_cast<double>(jz) / cfg.nz;
			Hxy_coscos = 0;
			Hyy_sincos = 0;
			Hzy_sinsin = 0;
			for (int x = 0; x < cfg.nx; x++) {
				for (int z = 0; z < cfg.nz; z++) {
					Hxy_coscos += (dyn.h1[x][cfg.ny - 1][z][7] - dyn.h1[x][cfg.ny - 1][z][14])
						* 2. * cos(qx * (x - x_Sph)) * cos(qz * (z - z_Sph)) / static_cast<double>(cfg.nx * cfg.nz);
					Hyy_sincos += (dyn.h1[x][cfg.ny - 1][z][7] + dyn.h1[x][cfg.ny - 1][z][14] +
						dyn.h1[x][cfg.ny - 1][z][11] + dyn.h1[x][cfg.ny - 1][z][17] +
						dyn.h1[x][cfg.ny - 1][z][3])
						* 2. * sin(qx * (x - x_Sph)) * cos(qz * (z - z_Sph)) / static_cast<double>(cfg.nx * cfg.nz);
					Hzy_sinsin += (dyn.h1[x][cfg.ny - 1][z][11] - dyn.h1[x][cfg.ny - 1][z][17])
						* 2. * sin(qx * (x - x_Sph)) * sin(qz * (z - z_Sph)) / static_cast<double>(cfg.nx * cfg.nz);
				}
			}
			complex<double> u_h_vec[3] =
			{ dyn.H_to_uT[jx][jz][0][0] * Hxy_coscos + dyn.H_to_uT[jx][jz][0][1] * Hyy_sincos + dyn.H_to_uT[jx][jz][0][2] * Hzy_sinsin,
			  dyn.H_to_uT[jx][jz][1][0] * Hxy_coscos + dyn.H_to_uT[jx][jz][1][1] * Hyy_sincos + dyn.H_to_uT[jx][jz][1][2] * Hzy_sinsin,
			  dyn.H_to_uT[jx][jz][2][0] * Hxy_coscos + dyn.H_to_uT[jx][jz][2][1] * Hyy_sincos + dyn.H_to_uT[jx][jz][2][2] * Hzy_sinsin };
			for (int x = 0; x < cfg.nx; x++) {
				for (int z = 0; z < cfg.nz; z++) {
					dyn.uxT[x][z] += u_h_vec[0] * cos(qx * (x - x_Sph)) * cos(qz * (z - z_Sph));
					dyn.uyT[x][z] += u_h_vec[1] * sin(qx * (x - x_Sph)) * cos(qz * (z - z_Sph));
					dyn.uzT[x][z] += u_h_vec[2] * sin(qx * (x - x_Sph)) * sin(qz * (z - z_Sph));
				}
			}
		}
	}
}

void Core3d::Write_Header_Filename() {
	ofstream header_filename_file;
	header_filename_file.open("header_filename.txt", std::ofstream::out | std::ofstream::trunc);
	header_filename_file << cfg.filename;
	header_filename_file.close();
}

void Core3d::Read_Header_Filename() {
	ifstream header_filename_file;
	header_filename_file.open("header_filename.txt", std::ofstream::in);
	header_filename_file >> cfg.filename;
	header_filename_file.close();
}

void Core3d::Write_Header_Info(string _i_sim_str) {
	ofstream header_info_file;
	header_info_file.open("./results/" + cfg.filename + "_HeaderInfo_" + _i_sim_str + ".dat", std::ofstream::out | std::ofstream::trunc);
	header_info_file << "[" << cfg.filename << "]" << endl;
	header_info_file << "char_time_fac" << "=" << cfg.char_time_fac << endl;
	header_info_file << "target_acc_Hz" << "=" << cfg.target_acc_Hz << endl;
	header_info_file << "rad_Sph_SI" << "=" << cfg.rad_Sph_SI << endl;
	header_info_file << "cellheight_SI" << "=" << cfg.cellheight_SI << endl;
	header_info_file << "gridreso_SI" << "=" << cfg.gridreso_SI << endl;
	header_info_file << "cell_diam_by_sphere_diam" << "=" << cfg.cell_diam_by_sphere_diam << endl;
	header_info_file << "frac_above_surf" << "=" << cfg.frac_above_surf<< endl;
	header_info_file << "nufac" << "=" << cfg.nufac << endl;
	header_info_file << "tandelinv" << "=" << cfg.tandelinv << endl;
	header_info_file << "iovt" << "=" << cfg.iovt << endl;
	header_info_file << "nx" << "=" << cfg.nx << endl;
	header_info_file << "ny" << "=" << cfg.ny << endl;
	header_info_file << "nz" << "=" << cfg.nz << endl;
	header_info_file << "outInterval" << "=" << cfg.outInterval << endl;
	header_info_file << "B_Advanced_Top" << "=" << cfg.B_Advanced_Top << endl;
	header_info_file << "jqmax" << "=" << cfg.jqmax << endl;
	header_info_file << "SphR" << "=" << cfg.SphR << endl;
	header_info_file << "Sphx" << "=" << cfg.Sphx << endl;
	header_info_file << "Sphy" << "=" << cfg.Sphy << endl;
	header_info_file << "Sphz" << "=" << cfg.Sphz << endl;
	header_info_file << "t_computation_secs" << "=" << cfg.t_computation_secs << endl;
	header_info_file.close();
}

void Core3d::Set_VE_Pars() {
	for (size_t x = 0; x < cfg.nx; x++) {
		for (size_t y = 0; y < cfg.ny; y++) {
			for (size_t z = 0; z < cfg.nz; z++) {
				double r = pow(pow(x - cfg.Sphx, 2) + pow(y - cfg.Sphy, 2) + pow(z - cfg.Sphz, 2), 0.5);
				if (r < cfg.SphR) {
					dyn.one_m_tauinv_m_Iom_arr[x][y][z] = cfg.one_m_tauinv_m_Iom_Sph;
					dyn.tauinv_m_om2_by_2_arr[x][y][z]  = cfg.tauinv_m_om2b2_Sph;
					if (y == cfg.ny - 1) {dyn.ZLiq_Top[x][z] = cfg.Zliq_Sph;}
				}
				else {
					dyn.one_m_tauinv_m_Iom_arr[x][y][z] = cfg.one_m_tauinv_m_Iom;
					dyn.tauinv_m_om2_by_2_arr[x][y][z]  = cfg.tauinv_m_om2b2;
					if (y == cfg.ny - 1) {dyn.ZLiq_Top[x][z] = cfg.Zliq;					}
				}
			}
		}
	}
}

void Core3d::WriteDisplacementField(string _i_sim_str) {
	ofstream outfile; 
	outfile.open("./results/" + cfg.filename + "_Displ_" + _i_sim_str + ".dat", std::ofstream::out | std::ofstream::trunc);
	for (int x = 0; x < cfg.nx; x++) {
		for (int y = 0; y < cfg.ny; y++) {
			outfile << x << "\t" << y << "\t" << dyn.ux[x][y][0].real() << "\t" << dyn.ux[x][y][0].imag() << 
			                             "\t" << dyn.uy[x][y][0].real() << "\t" << dyn.uy[x][y][0].imag() <<
				                         "\t" << dyn.uz[x][y][0].real() << "\t" << dyn.uz[x][y][0].imag() <<  endl;
		}
	}
	outfile.close();
}

void Core3d::WritePopulations(string _i_sim_str) {
	ofstream outfile;
	outfile.open("./results/" + cfg.filename + "_Pop_" + _i_sim_str + ".dat", std::ofstream::out | std::ofstream::trunc);
	for (int x = 0; x < cfg.nx; x++) {
		for (int y = 0; y < cfg.ny; y++) {
			for (int z = 0; z < cfg.nz; z++) {
				for (int i = 0; i < vst.ndir; i++) {
					outfile << x << "\t" << y << "\t" << "\t" << z << "\t" << i << "\t" <<
						dyn.h1[x][y][z][i].real() << "\t" << dyn.h1[x][y][z][i].imag() << endl;
				}
			}
		}
	}
	outfile.close();
};

void Core3d::Count_and_Allocate() {
	dyn.H_to_uT.resize(cfg.jqmax + 1);
	for (int jx = 0; jx < cfg.jqmax + 1; jx++) {
		dyn.H_to_uT[jx].resize(cfg.jqmax + 1);
		for (int jz = 0; jz < cfg.jqmax + 1; jz++) {
			dyn.H_to_uT[jx][jz].resize(3);
			for (int j123 = 0; j123 < 3; j123++) {
				dyn.H_to_uT[jx][jz][j123].resize(3);
			}
		}
	}


	dyn.dr.resize(cfg.nx);
	dyn.ux.resize(cfg.nx);
	dyn.uy.resize(cfg.nx);
	dyn.uz.resize(cfg.nx);
	dyn.ZLiq_Top.resize(cfg.nx);

	dyn.one_m_tauinv_m_Iom_arr.resize(cfg.nx);
	dyn.tauinv_m_om2_by_2_arr.resize(cfg.nx);
	dyn.one_m_tauinv_p_m_Iom_arr.resize(cfg.nx);
	dyn.tauinv_p_m_om2_by_2_arr.resize(cfg.nx);
	dyn.one_m_tauinv_m_m_Iom_arr.resize(cfg.nx);
	dyn.tauinv_m_m_om2_by_2_arr.resize(cfg.nx);
	dyn.h1.resize(cfg.nx);
	dyn.h2.resize(cfg.nx);

	dyn.uxT.resize(cfg.nx);
	dyn.uyT.resize(cfg.nx);
	dyn.uzT.resize(cfg.nx);


	for (int x = 0; x < cfg.nx; x++) {
		dyn.uxT[x].resize(cfg.nz);
		dyn.uyT[x].resize(cfg.nz);
		dyn.uzT[x].resize(cfg.nz);
		dyn.ZLiq_Top[x].resize(cfg.nz);
		dyn.dr[x].resize(cfg.ny);
		dyn.ux[x].resize(cfg.ny);
		dyn.uy[x].resize(cfg.ny);
		dyn.uz[x].resize(cfg.ny);
		dyn.one_m_tauinv_m_Iom_arr[x].resize(cfg.ny);
		dyn.tauinv_m_om2_by_2_arr[x].resize(cfg.ny);
		dyn.h1[x].resize(cfg.ny);
		dyn.h2[x].resize(cfg.ny);
		for (int y = 0; y < cfg.ny; y++) {
			dyn.dr[x][y].resize(cfg.nz);
			dyn.ux[x][y].resize(cfg.nz);
			dyn.uy[x][y].resize(cfg.nz);
			dyn.uz[x][y].resize(cfg.nz);
			dyn.one_m_tauinv_m_Iom_arr[x][y].resize(cfg.nz);
			dyn.tauinv_m_om2_by_2_arr[x][y].resize(cfg.nz);


			dyn.h1[x][y].resize(cfg.nz);
			dyn.h2[x][y].resize(cfg.nz);
			for (int z = 0; z < cfg.nz; z++) {
				dyn.h1[x][y][z].resize(vst.ndir);
				dyn.h2[x][y][z].resize(vst.ndir);
			}
		}
	}
}

void Core3d::Free_Mem() {
	for (int jx = 0; jx < cfg.jqmax + 1; jx++) {
		for (int jz = 0; jz < cfg.jqmax + 1; jz++) {
			for (int j123 = 0; j123 < 3; j123++) {
				dyn.H_to_uT[jx][jz][j123].resize(0);
			}
			dyn.H_to_uT[jx][jz].resize(0);
		}
		dyn.H_to_uT[jx].resize(0);
	}
	dyn.H_to_uT.resize(0);

	for (int x = 0; x < cfg.nx; x++) {
		for (int y = 0; y < cfg.ny; y++) {
			for (int z = 0; z < cfg.nz; z++) {
				dyn.h1[x][y][z].resize(0);
				dyn.h2[x][y][z].resize(0);
			}
			dyn.h1[x][y].resize(0);
			dyn.h2[x][y].resize(0);
			dyn.dr[x][y].resize(0);
			dyn.ux[x][y].resize(0);
			dyn.uy[x][y].resize(0);
			dyn.uz[x][y].resize(0);
			dyn.one_m_tauinv_m_Iom_arr[x][y].resize(0);
			dyn.tauinv_m_om2_by_2_arr[x][y].resize(0);
		}
		dyn.h1[x].resize(0);
		dyn.h2[x].resize(0);
		dyn.dr[x].resize(0);
		dyn.ux[x].resize(0);
		dyn.uy[x].resize(0);
		dyn.uz[x].resize(0);
		dyn.one_m_tauinv_m_Iom_arr[x].resize(0);
		dyn.tauinv_m_om2_by_2_arr[x].resize(0);
		dyn.ZLiq_Top[x].resize(0);
		dyn.uxT[x].resize(0);
		dyn.uyT[x].resize(0);
		dyn.uzT[x].resize(0);
	}
	dyn.h1.resize(0);
	dyn.h2.resize(0);
	dyn.dr.resize(0);
	dyn.ux.resize(0);
	dyn.uy.resize(0);
	dyn.uz.resize(0);
	dyn.ZLiq_Top.resize(0);
	dyn.uxT.resize(0);
	dyn.uyT.resize(0);
	dyn.uzT.resize(0);
	dyn.one_m_tauinv_m_Iom_arr.resize(0);
	dyn.tauinv_m_om2_by_2_arr.resize(0);
}


complex<double> Core3d::Stream_LB() {
	complex<double> FxLiqOnWall = 0.0;
	for (int x = 0; x < cfg.nx; x++) {
		for (int y = 0; y < cfg.ny; y++) {
			for (int z = 0; z < cfg.nz; z++) {
				for (int i = 0; i < vst.ndir; i++) {
					int ibar = vst.no_slip[i];
					int xmd = (cfg.nx + x - vst.dirx[i]) % cfg.nx;
					int ymd = y - vst.diry[i];
					int zmd = (cfg.nz + z - vst.dirz[i]) % cfg.nz;
					if ((ymd < cfg.ny) && (ymd >= 0)) {
						dyn.h2[x][y][z][i] = dyn.h1[xmd][ymd][zmd][i];
					}
					if (ymd >= cfg.ny) {
						if (cfg.B_Advanced_Top) {
							dyn.h2[x][y][z][i] = dyn.h1[x][y][z][ibar]
								- 6.0 * vst.wi[ibar] * vst.dirx[ibar] * dyn.uxT[x][z]
								- 6.0 * vst.wi[ibar] * vst.diry[ibar] * dyn.uyT[x][z]
								- 6.0 * vst.wi[ibar] * vst.dirz[ibar] * dyn.uzT[x][z];
						}
						else {
							complex<double> uxTop = 2.0 * (dyn.h1[x][y][z][7] - dyn.h1[x][y][z][14]) / (dyn.ZLiq_Top[x][z] + 1.0 / 3.0);
							complex<double> uzTop = 2.0 * (dyn.h1[x][y][z][11] - dyn.h1[x][y][z][17]) / (dyn.ZLiq_Top[x][z] + 1.0 / 3.0);
							dyn.h2[x][y][z][i] = dyn.h1[x][y][z][ibar] - 6.0 * vst.wi[ibar] * vst.dirx[ibar] * uxTop
								- 6.0 * vst.wi[ibar] * vst.dirz[ibar] * uzTop;
						}
					}
					if (ymd < 0) {
						dyn.h2[x][y][z][i] = dyn.h1[x][y][z][ibar] - 6.0 * vst.wi[ibar] * vst.dirx[ibar];
						FxLiqOnWall += (dyn.h2[x][y][z][i] + dyn.h1[x][y][z][ibar]) * static_cast<double>(vst.dirx[i]);
					}
				}
			}
		}
	}
	return FxLiqOnWall;
}

void Core3d::Compute_drho_u() {
	for (size_t x = 0; x < cfg.nx; x++) {
		for (size_t y = 0; y < cfg.ny; y++) {
			for (size_t z = 0; z < cfg.nz; z++) {
				dyn.dr[x][y][z] = 0.0;
				dyn.ux[x][y][z] = 0.0;
				dyn.uy[x][y][z] = 0.0;
				dyn.uz[x][y][z] = 0.0;
				for (size_t i = 0; i < dyn.h1[x][y][z].size(); i++) {
					dyn.dr[x][y][z] += dyn.h2[x][y][z][i];
					dyn.ux[x][y][z] += dyn.h2[x][y][z][i] * static_cast<double>(vst.dirx[i]);
					dyn.uy[x][y][z] += dyn.h2[x][y][z][i] * static_cast<double>(vst.diry[i]);
					dyn.uz[x][y][z] += dyn.h2[x][y][z][i] * static_cast<double>(vst.dirz[i]);
				}
			}
		}			 
	}
}

void Core3d::Collide() {
	for (size_t x = 0; x < cfg.nx; x++) {
		for (size_t y = 0; y < cfg.ny; y++) {
			for (size_t z = 0; z < cfg.nz; z++) {
				for (size_t i = 0; i < vst.ndir; i++) {
					complex<double>cidotu = static_cast<double>(vst.dirx[i])* dyn.ux[x][y][z] 
							              + static_cast<double>(vst.diry[i])* dyn.uy[x][y][z]
						                  + static_cast<double>(vst.dirz[i])* dyn.uz[x][y][z];
					complex<double>heq = vst.wi[i] * (dyn.dr[x][y][z] + 3.0 * cidotu);
					dyn.h1[x][y][z][i] = dyn.one_m_tauinv_m_Iom_arr[x][y][z] * dyn.h2[x][y][z][i]
							+ dyn.tauinv_m_om2_by_2_arr[x][y][z] * heq;
				}
			}
		}
	}
}

void Core3d::init_h(){
	for (size_t x = 0; x < cfg.nx; x++) {
		for (size_t y = 0; y < cfg.ny; y++) {
			for (size_t z = 0; z < cfg.nz; z++) {
				for (size_t i = 0; i < vst.ndir; i++) {dyn.h1[x][y][z][i] = 0;}
			}
		}
	}
}

void Core3d::WriteDisplacementField_Sphere(string _i_sim_str) {
	ofstream outfile;
	outfile.open("./results/" + cfg.filename + "_Displ_Sph_" + _i_sim_str + ".dat", std::ofstream::out | std::ofstream::trunc);
	for (int x = 0; x < cfg.nx; x++) {
		for (int y = 0; y < cfg.ny; y++) {
			for (int z = 0; z < cfg.nz; z++) {
				double r = pow(pow(x - cfg.Sphx, 2) + pow(y - cfg.Sphy, 2) + pow(z - cfg.Sphz, 2), 0.5);
				if (r <= cfg.SphR) {
					outfile << x << "\t" << y << "\t" << z <<
						"\t" << dyn.ux[x][y][z].real() << "\t" << dyn.ux[x][y][z].imag() <<
						"\t" << dyn.uy[x][y][z].real() << "\t" << dyn.uy[x][y][z].imag() <<
						"\t" << dyn.uz[x][y][z].real() << "\t" << dyn.uz[x][y][z].imag() << endl;
				}
			}
		}
	}
	outfile.close();
}

void Core3d::Simulate(string _i_sim_str) {
	ofstream ring_in_file; ring_in_file.open("./results/" + cfg.filename +  "_ring_in_" + _i_sim_str + ".dat", std::ofstream::out | std::ofstream::trunc); ring_in_file.close();
	Calc_Parameters_from_Config();
	Count_and_Allocate();
	Set_VE_Pars();
	init_h();
	cout << _i_sim_str << endl;
	compute_Matrices();
	complex<double> dfc_by_n = 0;
	double tot_time;
	int n_steps_char = int(cfg.ny * cfg.ny / cfg.nu * cfg.char_time_fac);
	cout << "n_steps_char: " << n_steps_char  << endl;
	time_t t_ini_l = time(0);
	for (int it = 0; it < n_steps_char; it++) {
		if (cfg.B_Advanced_Top) { compute_uT();};
		complex<double> FxLiqOnWall = Stream_LB();
		Compute_drho_u();
		Collide();
		if (it % cfg.outInterval == 0) {
			time_t t_current_l = time(0);
			if (it > 0) {tot_time = (t_current_l - t_ini_l) * n_steps_char / it /60.; }
			else  tot_time = 0;
			dfc_by_n = cfg.f0_SI * cfg.IU / M_PI * FxLiqOnWall / cfg.Zq / cfg.ovtorder / static_cast<double>(cfg.nx * cfg.nz);
			if (it % (1 * cfg.outInterval) == 0) {
				cout << it << " / " << n_steps_char << " " << tot_time  << " [mins] , df/n = " << dfc_by_n.real() << " Hz, " << "dG/n = " << dfc_by_n.imag() << " Hz, " << endl;
			}
			ring_in_file.open("./results/" + cfg.filename + "_ring_in_" + _i_sim_str + ".dat", std::ofstream::out | std::ofstream::app);
			ring_in_file << dfc_by_n.real() << "\t" << dfc_by_n.imag() << endl;
			ring_in_file.close();
		}
	}
	cout << "df/n = " << dfc_by_n.real() << " Hz, " << "dG/n = " << dfc_by_n.imag() << " Hz" << endl;
	WriteDisplacementField(_i_sim_str);
	WriteDisplacementField_Sphere(_i_sim_str);
	WritePopulations(_i_sim_str);
	Free_Mem();
}

void Core3d::run(){
	Set_Config();
	Read_Header_Filename();
	char date[9];       // Displays current date (mm/dd/yy)
	_strdate_s(date);   // Displays current date (mm/dd/yy)

	//time_t _tm = time(NULL);
	//struct _tm* curtime = localtime_s(&_tm);
	//cout << "The current date/time is:" << asctime(curtime);

	int strlength = cfg.filename.length();
	string numberstr = cfg.filename.substr(strlength - 3, 3);
	string descriptor_old = cfg.filename.substr(0,strlength - 3);
	int number = stoi(numberstr);
	string descriptor_new = "test";
	if (descriptor_new == descriptor_old) { number++;}
	else { number = 0; };
	std::stringstream ss; ss << std::setw(3) << std::setfill('0') << number;
	cfg.filename = descriptor_new + ss.str();
	Write_Header_Filename();
	Set_VE_Pars();
	int check = _mkdir("./results");
	int isim = 0;
	int nsim_out = 9;
	int nsim_in  = 1;
	cout << "START OF RUN:  " << date << endl;   //displays current date (mm/dd/yy)	
	for (int isim_out = 0; isim_out < nsim_out; isim_out ++) {
		//cfg.cellheight_SI = 7e-9 + (2.e-9 * static_cast<double>(isim_out));
		cout << "cell_height: " << cfg.cellheight_SI << endl;
		//cout << "gridreso_SI: " << cfg.gridreso_SI << endl;
		//cfg.tandelinv = pow(10, static_cast<double>(isim_out) / 5.);
		cout << "tandelinv: " << cfg.tandelinv << endl;	
		//cfg.jqmax = 1 + (static_cast<double>(isim_out)); 
		cout << "jq_max: " << cfg.jqmax << endl;
		//cfg.gridreso_SI = (2.0 - (0.25 * static_cast<double>(isim_out))) * 1e-9;    //SI Value
		cout << "grid resolution: " << cfg.gridreso_SI << endl;
		    if (isim_out == 0) {
				cfg.nufac = 0;
				cout << "nufac: " << cfg.nufac << endl;
			}
			if (isim_out == 1) {
				cfg.nufac = 1;
				cout << "nufac: " << cfg.nufac << endl;
			}
			if (isim_out == 2) {
				cfg.nufac = 3;
				cout << "nufac: " << cfg.nufac << endl;
			}	
			if (isim_out == 3) {
				cfg.nufac = 5;
				cout << "nufac: " << cfg.nufac << endl;
			}
			if (isim_out == 4) {
				cfg.nufac = 10;
				cout << "nufac: " << cfg.nufac << endl;
			}
			if (isim_out == 5) {
				cfg.nufac = 20;
				cout << "nufac: " << cfg.nufac << endl;
			}
			if (isim_out == 6) {
				cfg.nufac = 50;
				cout << "nufac: " << cfg.nufac << endl;
			}
			if (isim_out == 7) {
				cfg.nufac = 100;
				cout << "nufac: " << cfg.nufac << endl;
			}
			if (isim_out == 8) {
				cfg.nufac = 1000;
				cout << "nufac: " << cfg.nufac << endl;
			}
		//cfg.nufac = pow(10,static_cast<double>(isim_out)/2.);
		//cout << "nufac: " << cfg.nufac << endl;
		//cfg.cell_diam_by_sphere_diam = 1.1 + (0.2 * static_cast<double>(isim_out));
		cout << "cell_diam_by_sphere_diam: " << cfg.cell_diam_by_sphere_diam << endl;
		cout << "frac_above_surface: " << cfg.frac_above_surf << endl;
//		cout << "cfg.char_time_fac: " << cfg.char_time_fac << endl;
		for (int isim_in = 0; isim_in < nsim_in; isim_in++) {
			//cfg.jqmax = isim_in;
			//cfg.iovt = isim_in;
      		//cfg.nufac = static_cast<double>(isim_in);
			cout << "iovt: " << cfg.iovt << endl;
			//cout << "nufac: " << cfg.nufac << endl;
			cout << "tandelinv: " << cfg.tandelinv << endl;
			time_t t_ini = time(0);
//			ss << std::setw(3) << std::setfill('0') << isim;
			std::string isim_str;
			/*if (isim_in == 0) {
				cfg.nufac = 2;
				cout << "nufac: " << cfg.nufac << endl;
			}
			if (isim_in == 1) {
				cfg.nufac = 5;
				cout << "nufac: " << cfg.nufac << endl;
			}
			if (isim_in == 2) {
				cfg.nufac = 100;
				cout << "nufac: " << cfg.nufac << endl;
			}*/
			for (int i_ref = 0; i_ref < 2; i_ref++) {
				if (i_ref == 0) {
					std::stringstream ss;
					ss << std::setw(3) << std::setfill('0') << isim;
					isim_str = ss.str();
					cfg.B_Reference = false;                //water
					cout << "sphere" << endl;
				}
				if (i_ref == 1) {
					std::stringstream ss;
					ss << std::setw(3) << std::setfill('0') << isim << "_ref";
					isim_str = ss.str();
					cfg.B_Reference = true;                //water
					cout << "water" << endl;
				}
				Calc_Parameters_from_Config();
				cout << "nx: " << cfg.nx << " ny: " << cfg.ny << " nz: " << cfg.nz << endl;
				Write_Header_Info(isim_str);
				Simulate(isim_str);
				time_t t_fin = time(0);
				cfg.t_computation_secs = t_fin - t_ini;
				cout << "Computation time : " << cfg.t_computation_secs << " secs" << endl;
				Write_Header_Info(isim_str);

			}
			isim++;
		}
	}
}

void Core3d::Set_Config() {
	cfg.filename = "test";
	cfg.char_time_fac = 3;             //factor of ny^2 wrt nu
	cfg.cell_diam_by_sphere_diam = 4.43;
	cfg.gridreso_SI = 3.0e-9;         //overwritten
	cfg.cellheight_SI = 13e-9;
	cfg.frac_above_surf = 0.8;  
	cfg.rad_Sph_SI = 5e-9;
	cfg.target_acc_Hz = 1;
	cfg.nufac = 1000;                 //very very very hard sphere // overwritten by the loop
	cfg.tandelinv = 0.0;            //overwritten
	cfg.iovt = 5;
	cfg.jqmax = 3;                 //changed to 3
	cfg.B_Advanced_Top = true;    //Sticking with always true
	cfg.outInterval = 5;
	cfg.f0_SI = 5e6;
	cfg.Zq_SI = 8.8e6;
	cfg.dq_SI = 332.075e-6;
	cfg.rhoq = 2.65;
	cfg.nu_SI = 1.0e-6;
	cfg.nu = 1.0 / 6.0;
	cfg.rho0 = 1.0;
}