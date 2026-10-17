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
	cfg.ovtorder = 2. * cfg.iovt - 1; //needs to be 2*cfg.rad_Sph_SI @ line 18
	cfg.nx = static_cast<int>(2*cfg.rad_Sph_SI / cfg.gridreso_SI * cfg.cell_diam_by_sphere_diam); 
	cfg.nz = cfg.nx;
	cfg.ny = static_cast<int>(cfg.cellheight_SI / cfg.gridreso_SI);
	cfg.nu = 1.0 / 6.0 * cfg.nufac * (1.0 - cfg.IU * cfg.tandelinv);
	cfg.tauinv =   2.0 / (6.0 * cfg.nu  + 1.0);
	cfg.tauinv_0 = 2.0 / (6.0 * cfg.nu0 + 1.0);

	cfg.tauinv_p   = cfg.tauinv;
	cfg.tauinv_p_0 = cfg.tauinv_0;
	cfg.tauinv_m   = 1. / ((cfg.Lambda / (1.0 / cfg.tauinv_p   - 0.5) + 0.5));
	cfg.tauinv_m_0 = 1. / ((cfg.Lambda / (1.0 / cfg.tauinv_p_0 - 0.5) + 0.5));


	cfg.delta = cfg.delta_SI / cfg.gridreso_SI;
	cfg.omega = 2.0 * cfg.nu0 / (cfg.delta * cfg.delta);
	cfg.omega2by2 = 1.0 / 2.0 * cfg.omega * cfg.omega;
	cfg.one_m_tauinv_m_Iom   = 1. - cfg.tauinv   - cfg.IU * cfg.omega;
	cfg.one_m_tauinv_p_m_Iom = 1. - cfg.tauinv_p - cfg.IU * cfg.omega;
	cfg.one_m_tauinv_m_m_Iom = 1. - cfg.tauinv_m - cfg.IU * cfg.omega;

	cfg.one_m_tauinv_m_Iom_0 = 1. - cfg.tauinv_0 - cfg.IU * cfg.omega;
	cfg.one_m_tauinv_p_m_Iom_0 = 1. - cfg.tauinv_p_0 - cfg.IU * cfg.omega;
	cfg.one_m_tauinv_m_m_Iom_0 = 1. - cfg.tauinv_m_0 - cfg.IU * cfg.omega;
	

	cfg.tauinv_m_om2_by_2   = cfg.tauinv   - cfg.omega * cfg.omega / 2.;
	cfg.tauinv_p_m_om2_by_2 = cfg.tauinv_p - cfg.omega * cfg.omega / 2.;
	cfg.tauinv_m_m_om2_by_2 = cfg.tauinv_m - cfg.omega * cfg.omega / 2.;

	cfg.tauinv_m_om2_by_2_0   = cfg.tauinv_0   - cfg.omega * cfg.omega / 2.;
	cfg.tauinv_p_m_om2_by_2_0 = cfg.tauinv_p_0 - cfg.omega * cfg.omega / 2.;
	cfg.tauinv_m_m_om2_by_2_0 = cfg.tauinv_m_0 - cfg.omega * cfg.omega / 2.;

	cfg.Zliq  = pow((cfg.IU * cfg.rho0 * cfg.omega * cfg.rho0 * cfg.nu ), 0.5);
	cfg.Zliq0 = pow((cfg.IU * cfg.rho0 * cfg.omega * cfg.rho0 * cfg.nu0), 0.5);
	cfg.c_theo = pow((cfg.IU * cfg.omega * cfg.nu), 0.5);
	cfg.k_theo = cfg.omega / cfg.c_theo;
	cfg.f0 = cfg.omega / (2.0 * M_PI);
	cfg.dq = cfg.dq_SI / cfg.gridreso_SI;
	cfg.Zq = 2.0 * cfg.dq * cfg.f0 * cfg.rhoq;
	cfg.nsteps_max = static_cast<int>(cfg.ncycles_max / cfg.omega);
	cfg.SphR = cfg.rad_Sph_SI / cfg.gridreso_SI;
	cfg.Sphx = static_cast<double>(cfg.nx - 1) / 2.;
	cfg.Sphz = static_cast<double>(cfg.nz - 1) / 2.;
	cfg.Sphy = 2 * cfg.SphR * (cfg.frac_above_surf - 0.5);
	//cout << "cfg.Sphy " << cfg.Sphy << endl;
	//cout << "cfg.frac_above_surf " << cfg.frac_above_surf << endl;
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
	header_info_file << "ncycles_max" << "=" << cfg.ncycles_max << endl;
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
	header_info_file << "nz" << "=" << cfg.nz << endl;
	header_info_file << "outInterval" << "=" << cfg.outInterval << endl;
	header_info_file << "B_TRT" << "=" << cfg.B_TRT << endl;
	header_info_file << "Lambda" << "=" << cfg.Lambda << endl;
	header_info_file << "t_computation_secs" << "=" << cfg.t_computation_secs << endl;

	header_info_file.close();
}

void Core3d::Set_VE_Pars() {
	for (size_t x = 0; x < cfg.nx; x++) {
		for (size_t y = 0; y < cfg.ny; y++) {
			for (size_t z = 0; z < cfg.nz; z++) {
				double r = pow(pow(x - cfg.Sphx, 2) + pow(y - cfg.Sphy, 2) + pow(z - cfg.Sphz, 2), 0.5);
				if (r < cfg.SphR) {
					dyn.one_m_tauinv_m_Iom_arr[x][y][z] = cfg.one_m_tauinv_m_Iom;
					dyn.tauinv_m_om2_by_2_arr[x][y][z] = cfg.tauinv_m_om2_by_2;

					dyn.one_m_tauinv_p_m_Iom_arr[x][y][z] = cfg.one_m_tauinv_p_m_Iom;
					dyn.tauinv_p_m_om2_by_2_arr[x][y][z]  = cfg.tauinv_p_m_om2_by_2;
					dyn.one_m_tauinv_m_m_Iom_arr[x][y][z] = cfg.one_m_tauinv_m_m_Iom;
					dyn.tauinv_m_m_om2_by_2_arr[x][y][z]  = cfg.tauinv_m_m_om2_by_2;
					if (y == cfg.ny - 1) {
						dyn.ZLiq_Top[x][z] = cfg.Zliq;
					}
				}
				else {
					dyn.one_m_tauinv_m_Iom_arr[x][y][z] = cfg.one_m_tauinv_m_Iom_0;
					dyn.tauinv_m_om2_by_2_arr[x][y][z] = cfg.tauinv_m_om2_by_2_0;

					dyn.one_m_tauinv_p_m_Iom_arr[x][y][z] = cfg.one_m_tauinv_p_m_Iom_0;
					dyn.tauinv_p_m_om2_by_2_arr[x][y][z] = cfg.tauinv_p_m_om2_by_2_0;
					dyn.one_m_tauinv_m_m_Iom_arr[x][y][z] = cfg.one_m_tauinv_m_m_Iom_0;
					dyn.tauinv_m_m_om2_by_2_arr[x][y][z] = cfg.tauinv_m_m_om2_by_2_0;
					if (y == cfg.ny - 1) {
						dyn.ZLiq_Top[x][z] = cfg.Zliq0;
					}
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

void Core3d::Count_and_Allocate() {
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
	for (int x = 0; x < cfg.nx; x++) {
		dyn.ZLiq_Top[x].resize(cfg.nz);
		dyn.dr[x].resize(cfg.ny);
		dyn.ux[x].resize(cfg.ny);
		dyn.uy[x].resize(cfg.ny);
		dyn.uz[x].resize(cfg.ny);
		dyn.one_m_tauinv_m_Iom_arr[x].resize(cfg.ny);
		dyn.tauinv_m_om2_by_2_arr[x].resize(cfg.ny);
		dyn.one_m_tauinv_p_m_Iom_arr[x].resize(cfg.ny);
		dyn.tauinv_p_m_om2_by_2_arr[x].resize(cfg.ny);
		dyn.one_m_tauinv_m_m_Iom_arr[x].resize(cfg.ny);
		dyn.tauinv_m_m_om2_by_2_arr[x].resize(cfg.ny);
		dyn.h1[x].resize(cfg.ny);
		dyn.h2[x].resize(cfg.ny);
		for (int y = 0; y < cfg.ny; y++) {
			dyn.dr[x][y].resize(cfg.nz);
			dyn.ux[x][y].resize(cfg.nz);
			dyn.uy[x][y].resize(cfg.nz);
			dyn.uz[x][y].resize(cfg.nz);
			dyn.one_m_tauinv_m_Iom_arr[x][y].resize(cfg.nz);
			dyn.tauinv_m_om2_by_2_arr[x][y].resize(cfg.nz);

			dyn.one_m_tauinv_p_m_Iom_arr[x][y].resize(cfg.nz);
			dyn.tauinv_p_m_om2_by_2_arr[x][y].resize(cfg.nz);
			dyn.one_m_tauinv_m_m_Iom_arr[x][y].resize(cfg.nz);
			dyn.tauinv_m_m_om2_by_2_arr[x][y].resize(cfg.nz);

			dyn.h1[x][y].resize(cfg.nz);
			dyn.h2[x][y].resize(cfg.nz);
			for (int z = 0; z < cfg.nz; z++) {
				dyn.h1[x][y][z].resize(vst.ndir);
				dyn.h2[x][y][z].resize(vst.ndir);
			}
		}
	}
	dyn.y_arr.resize(cfg.ny);
	dyn.y_SI_arr.resize(cfg.ny);
	dyn.ux_theo.resize(cfg.ny);
	dyn.h1_lo_z_real.resize(cfg.nx * cfg.ny * vst.ndir);
	dyn.h1_lo_z_imag.resize(cfg.nx * cfg.ny * vst.ndir);
	dyn.h1_hi_z_real.resize(cfg.nx * cfg.ny * vst.ndir);
	dyn.h1_hi_z_imag.resize(cfg.nx * cfg.ny * vst.ndir);
}

void Core3d::Free_Mem() {
	dyn.y_arr.resize(0);
	dyn.y_SI_arr.resize(0);

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

			dyn.one_m_tauinv_p_m_Iom_arr[x][y].resize(0);
			dyn.tauinv_p_m_om2_by_2_arr[x][y].resize(0);
			dyn.one_m_tauinv_m_m_Iom_arr[x][y].resize(0);
			dyn.tauinv_m_m_om2_by_2_arr[x][y].resize(0);
		}
		dyn.h1[x].resize(0);
		dyn.h2[x].resize(0);
		dyn.dr[x].resize(0);
		dyn.ux[x].resize(0);
		dyn.uy[x].resize(0);
		dyn.uz[x].resize(0);
		dyn.one_m_tauinv_m_Iom_arr[x].resize(0);
		dyn.tauinv_m_om2_by_2_arr[x].resize(0);
		dyn.one_m_tauinv_p_m_Iom_arr[x].resize(0);
		dyn.tauinv_p_m_om2_by_2_arr[x].resize(0);
		dyn.one_m_tauinv_m_m_Iom_arr[x].resize(0);
		dyn.tauinv_m_m_om2_by_2_arr[x].resize(0);
		dyn.ZLiq_Top[x].resize(0);
	}
	dyn.h1.resize(0);
	dyn.h2.resize(0);
	dyn.dr.resize(0);
	dyn.ux.resize(0);
	dyn.uy.resize(0);
	dyn.uz.resize(0);
	dyn.ZLiq_Top.resize(0);
	dyn.one_m_tauinv_m_Iom_arr.resize(0);
	dyn.tauinv_m_om2_by_2_arr.resize(0);

	dyn.one_m_tauinv_p_m_Iom_arr.resize(0);
	dyn.tauinv_p_m_om2_by_2_arr.resize(0);
	dyn.one_m_tauinv_m_m_Iom_arr.resize(0);
	dyn.tauinv_m_m_om2_by_2_arr.resize(0);
	dyn.y_arr.resize(0);
	dyn.y_SI_arr.resize(0);
	dyn.ux_theo.resize(0);
	dyn.h1_lo_z_real.resize(0);
	dyn.h1_lo_z_imag.resize(0);
	dyn.h1_hi_z_real.resize(0);
	dyn.h1_hi_z_imag.resize(0);
}

void Core3d::Export_to_arrays_for_MPI() {
	int iPo = 0;
	for (int x = 0; x < cfg.nx; x++) {
		for (int y = 0; y < cfg.ny; y++) {
			for (int i = 0; i < vst.ndir; i++) {
				dyn.h1_lo_z_real[iPo] = dyn.h1[x][y][0][i].real();
				dyn.h1_lo_z_imag[iPo] = dyn.h1[x][y][0][i].imag();
				dyn.h1_hi_z_real[iPo] = dyn.h1[x][y][cfg.nz - 1][i].real();
				dyn.h1_hi_z_imag[iPo] = dyn.h1[x][y][cfg.nz - 1][i].imag();
				iPo++;
			}
		}
	}

}

void Core3d::Import_from_arrays_for_MPI(){
	int iPo = 0;
	for (int x = 0; x < cfg.nx; x++) {
		for (int y = 0; y < cfg.ny; y++) {
			for (int i = 0; i < vst.ndir; i++) {
				dyn.h1[x][y][0][i]          = dyn.h1_lo_z_real[iPo] + cfg.IU * dyn.h1_lo_z_imag[iPo];
				dyn.h1[x][y][cfg.nz - 1][i] = dyn.h1_hi_z_real[iPo] + cfg.IU * dyn.h1_hi_z_imag[iPo];
				iPo++;
			}
		}
	}
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
						continue;
					}
					if (ymd >= cfg.ny) {
						complex<double> uxTop = 2.0 * (dyn.h1[x][y][z][7 ] - dyn.h1[x][y][z][14]) / (dyn.ZLiq_Top[x][z] + 1.0 / 3.0);
						complex<double> uzTop = 2.0 * (dyn.h1[x][y][z][11] - dyn.h1[x][y][z][17]) / (dyn.ZLiq_Top[x][z] + 1.0 / 3.0);
						dyn.h2[x][y][z][i] = dyn.h1[x][y][z][ibar] - 6.0 * vst.wi[ibar] * vst.dirx[ibar] * uxTop
							                                       - 6.0 * vst.wi[ibar] * vst.dirz[ibar] * uzTop;
						continue;
					}
					if (ymd < 0) {
						dyn.h2[x][y][z][i] = dyn.h1[x][y][z][ibar] - 6.0 * vst.wi[ibar] * vst.dirx[ibar];
//						dyn.h2[x][y][z][i] = dyn.h1[x][y][z][ibar] - 6.0 * vst.wi[ibar] * vst.dirz[ibar];
						FxLiqOnWall += (dyn.h2[x][y][z][i] + dyn.h1[x][y][z][ibar]) * static_cast<double>(vst.dirx[i]);
//						FxLiqOnWall += (dyn.h2[x][y][z][i] + dyn.h1[x][y][z][ibar]) * static_cast<double>(vst.dirz[i]);
						continue;
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

void Core3d::Calculate_h1_from_TRT(int x, int y, int z) {
	complex<double> heq[19], cidotu, hp[10], hm[10], heqp[10], heqm[10],
		hp_collided[10], hm_collided[10];
	const int i_no_bar_list[10] = { 0,1,3,5,7,9,11,13,15,17 };

	for (int i = 0; i < vst.ndir; i++) {
		cidotu = static_cast<double>(vst.dirx[i])* dyn.ux[x][y][z]
			+ static_cast<double>(vst.diry[i])* dyn.uy[x][y][z]
			+ static_cast<double>(vst.dirz[i])* dyn.uz[x][y][z];
		heq[i] = vst.wi[i] * (dyn.dr[x][y][z] + 3.0 * cidotu);
	}
	for (int j = 0; j < 10; j ++) {
		int i    = i_no_bar_list[j];
		int ibar = vst.no_slip[i_no_bar_list[j]];


		hp[j] = (dyn.h2[x][y][z][i] + dyn.h2[x][y][z][ibar]) / 2.;
		hm[j] = (dyn.h2[x][y][z][i] - dyn.h2[x][y][z][ibar]) / 2.;
		heqp[j] = (heq[i] + heq[ibar]) / 2.;
		heqm[j] = (heq[i] - heq[ibar]) / 2.;
		hp_collided[j] = dyn.one_m_tauinv_p_m_Iom_arr[x][y][z] * hp[j] +
						 dyn.tauinv_p_m_om2_by_2_arr[x][y][z] * heqp[j];
		hm_collided[j] = dyn.one_m_tauinv_m_m_Iom_arr[x][y][z] * hm[j] +
						 dyn.tauinv_m_m_om2_by_2_arr[x][y][z] * heqm[j];
		dyn.h1[x][y][z][i]    = hp_collided[j] + hm_collided[j];
		dyn.h1[x][y][z][ibar] = hp_collided[j] - hm_collided[j];
	}
}

void Core3d::Collide() {
	for (size_t x = 0; x < cfg.nx; x++) {
		for (size_t y = 0; y < cfg.ny; y++) {
			for (size_t z = 0; z < cfg.nz; z++) {
				if (cfg.B_TRT == false) {
					for (size_t i = 0; i < vst.ndir; i++) {
						complex<double>cidotu = static_cast<double>(vst.dirx[i])* dyn.ux[x][y][z] 
							                  + static_cast<double>(vst.diry[i])* dyn.uy[x][y][z]
						                      + static_cast<double>(vst.dirz[i])* dyn.uz[x][y][z];
						complex<double>heq = vst.wi[i] * (dyn.dr[x][y][z] + 3.0 * cidotu);
						dyn.h1[x][y][z][i] = dyn.one_m_tauinv_m_Iom_arr[x][y][z] * dyn.h2[x][y][z][i]
								+ dyn.tauinv_m_om2_by_2_arr[x][y][z] * heq;
					}
				}
				else {
					Calculate_h1_from_TRT(x, y, z);
				}
			}
		}
	}
}

void Core3d::init_h_as_wave() 
{
	for (size_t x = 0; x < cfg.nx; x++) {
		for (size_t y = 0; y < cfg.ny; y++) {
			for (size_t z = 0; z < cfg.nz; z++) {
				complex<double> ux_loc = exp(-cfg.k_theo * static_cast<double>(y));
				complex<double> uy_loc = 0;
				complex<double> uz_loc = 0;
				complex<double> dr_loc = 0;

				for (size_t i = 0; i < vst.ndir; i++) {
					complex<double>cidotu = vst.dirxD[i] * ux_loc + vst.diryD[i] * uy_loc + vst.dirzD[i] * uz_loc;
//					dyn.h1[x][y][z][i] = vst.wi[i] * (dr_loc + 3.0 * cidotu);
					dyn.h1[x][y][z][i] = 0;

				}
			}
		}
	}
}

void Core3d::Simulate(string _i_sim_str) {
	ofstream ring_in_file; ring_in_file.open("./results/" + cfg.filename +  "_ring_in_" + _i_sim_str + ".dat", std::ofstream::out | std::ofstream::trunc); ring_in_file.close();
	Count_and_Allocate();
	Set_VE_Pars();
	//dyn.init_h(); replaced on 19_06_2020
	init_h_as_wave();
	cout << _i_sim_str << endl;
	complex<double> FxLiqOnWall;
	complex<double> dfc_by_n = 0;
	int icount = 0;
	for (int it = 0; it < cfg.nsteps_max; it++) {
		FxLiqOnWall = Stream_LB();
		Compute_drho_u();
		Collide();
		Export_to_arrays_for_MPI();
		if (icount % cfg.outInterval == 0) {
			complex<double> dfc_by_n_prev = dfc_by_n;
			dfc_by_n = cfg.f0_SI * cfg.IU / M_PI * FxLiqOnWall / cfg.Zq / cfg.ovtorder / static_cast<double>(cfg.nx * cfg.nz);
			if (icount % (1 * cfg.outInterval) == 0) {
				cout << it << ", df/n = " << dfc_by_n.real() << " Hz, " << "dG/n = " << dfc_by_n.imag() << " Hz, " <<  
					"|change| = "  << abs(dfc_by_n_prev - dfc_by_n) << endl;
			}
			ring_in_file.open("./results/" + cfg.filename + "_ring_in_" + _i_sim_str + ".dat", std::ofstream::out | std::ofstream::app);
			ring_in_file << dfc_by_n.real() << "\t" << dfc_by_n.imag() << endl;
			ring_in_file.close();
			if ((abs(dfc_by_n_prev - dfc_by_n) < cfg.target_acc_Hz) or (icount > 3000)) {
				cout << "|final change| =  " << abs(dfc_by_n_prev - dfc_by_n) << endl;
				cout << "FxLiqOnWall " << FxLiqOnWall << endl;
				break;
			}
		}
		icount++;
	}
	cout << "df/n = " << dfc_by_n.real() << " Hz, " << "dG/n = " << dfc_by_n.imag() << " Hz" << endl;
	WriteDisplacementField(_i_sim_str);
	WriteDisplacementField_Sphere(_i_sim_str);
	Free_Mem();
}



void Core3d::WriteDisplacementField_Sphere(string _i_sim_str)  //Writing displacement fields to file.
{
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

void Core3d::run() 
{
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
	int nsim_out = 10;
	int nsim_in  = 8;
	cout << "START OF RUN:  " << date << endl;   //displays current date (mm/dd/yy)	
	for (int isim_out = 0; isim_out < nsim_out; isim_out++) {
		//cfg.nufac = 2 *(1.+ static_cast<double>(isim_out));
		//cfg.cell_diam_by_sphere_diam = 2. + static_cast<double>(isim_out) * 0.5;
		//cout << "cell_diam_by_sphere_diam: " << cfg.cell_diam_by_sphere_diam << endl;

		cfg.cellheight_SI = 15e-9 + 2.e-9 * static_cast<double>(isim_out);
		cout << "cell_height: " << cfg.cellheight_SI << endl;
		cout << "cell_diam_by_sphere_diam: " << cfg.cell_diam_by_sphere_diam << endl;
		cout << "gridreso_SI: " << cfg.gridreso_SI << endl;
		cout << "tandelinv: " << cfg.tandelinv << endl;
		for (int isim_in = 0; isim_in < nsim_in; isim_in++) {
			time_t t_ini = time(0);
			cfg.nufac = 2 * (1. + static_cast<double>(isim_out));
			//cfg.frac_above_surf = 0.1 + 0.8 * static_cast<double>(isim_in) 
			//	                            / static_cast<double>(nsim_in-1);
			cout << "nufac: " << cfg.nufac << endl;
			std::stringstream ss; ss << std::setw(3) << std::setfill('0') << isim;
			std::string isim_str = ss.str();
			Calc_Parameters_from_Config();
			cout << "nx: " << cfg.nx << " ny: " << cfg.ny << " nz: " << cfg.nz << endl;
			Write_Header_Info(isim_str);
			Simulate(isim_str);
			time_t t_fin = time(0);
			cfg.t_computation_secs = t_fin - t_ini;
			cout << "Computation time : " << cfg.t_computation_secs << " secs" << endl;
			Write_Header_Info(isim_str);
			isim++;
		}
	}
}


void Core3d::Set_Config() {
	cfg.filename = "test";
	cfg.B_TRT = false;
	cfg.Lambda = 1. / 4.;
	cfg.gridreso_SI = 3e-9;
	cfg.cellheight_SI = 15e-9;
	cfg.rad_Sph_SI = 10e-9;
	cfg.target_acc_Hz = 1;
	cfg.nufac = 3.0;
	cfg.tandelinv = 5.0;
	cfg.f0_SI = 5.0e6;
	cfg.Zq_SI = 8.8e6;
	cfg.dq_SI = 332.075e-6;
	cfg.rhoq = 2.65;
	cfg.nu_SI = 1.0e-6;
	cfg.delta_f0_SI = sqrt((2.0 * cfg.nu_SI / (2.0 * M_PI * cfg.f0_SI)));
	cfg.nu0 = 1.0 / 6.0;
	cfg.rho0 = 1.0;
	cfg.ncycles_max = 100;
	cfg.dq = cfg.dq_SI / cfg.gridreso_SI;
	cfg.cell_diam_by_sphere_diam = 1.1;
	cfg.frac_above_surf= 0.1;  //what we usually call Truncation
	cfg.iovt = 1;
	cfg.nx = 15;    //Will be overwritten in Calc_Parameters_from_Config()
	cfg.nz = 15;    //Will be overwritten in Calc_Parameters_from_Config()
	cfg.delta_SI = pow((2 * cfg.nu_SI / (2 * M_PI * cfg.f0_SI)), 0.5);
	cfg.outInterval = 100;
}