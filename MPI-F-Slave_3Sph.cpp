//
// MPI-F-Slave Slave-Implementation zum Test eines MPI-Kommunikationsprotokolls, parallelisiert mit MPI
// MPI-F-Slave.cpp
//
// Created on: 20.11.2019

#include <mpi.h>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <thread>
#include <chrono>
#include <math.h>
#define _USE_MATH_DEFINES
#include "MersenneTwister.h"
#include "MPI-F-Main.h"
#include "MPI-F-Slave.h"

using namespace std;

SlaveF::SlaveF(int _rank, int _nSlaves, bool _flagMPI, bool _flagDebug) {
	Set_Config_File_Input();
	
	Calc_Parameters_from_Config();
	Count_and_Allocate();
	Set_VE_Pars(_rank, _nSlaves);
	init_h();
	complex<double> FxLiqOnWall;
	complex<double> dfc_by_n = 0;
	int icount = 0;
	dataSize = cfg.nx * cfg.ny * vst.ndir * 2;
	dimensionSendRecv = dataSize;
	if (_flagMPI) {
		dataLeftSend.dArray = new double[dataSize];
		dataRightSend.dArray = new double[dataSize];
		dataLeftRecv.dArray = new double[dataSize];
		dataRightRecv.dArray = new double[dataSize];
		for (int i = 0; i < dataSize; i++) {
			dataLeftSend.dArray[i] = static_cast<double>(_rank * 10.0);
			dataRightSend.dArray[i] = static_cast<double>(_rank * 10.0);
			dataLeftRecv.dArray[i] = static_cast<double>(_rank * 10.0);
			dataRightRecv.dArray[i] = static_cast<double>(_rank * 10.0);
		}
		randS = new MTRand();
		randS->seed(_rank);
		request = new MPI_Request[1];
		request[0] = MPI_REQUEST_NULL;
		int size = 0;
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		int left = 0;
		if (_rank == 1) {
			left = size - 1;
		}
		else {
			left = _rank - 1;
		}
		int right = 0;
		if (_rank == (size - 1)) {
			right = 1;
		}
		else {
			right = _rank + 1;
		}
		MPI_Comm noComm;
		for (int i = 1; i < size; i++) {
			if (i == _rank) {
				MPI_Comm_split(MPI_COMM_WORLD, _rank, 0, &rightComm);
				continue;
			}
			if (i == left) {
				MPI_Comm_split(MPI_COMM_WORLD, left, 1, &leftComm);
				continue;
			}
			MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, 0, &noComm);
		}
		MPI_Comm_rank(rightComm, &rightRank);
		MPI_Comm_rank(leftComm, &leftRank);
		cout << "slave: " << _rank << "init" << " r -> " << right << " l -> " << left << " dataSize: " << dataSize << endl;
		fstream logfile("./results/" + cfg.filename + "_ring_in_" + to_string(_rank) + ".dat", std::ofstream::out);
	}
}

bool SlaveF::run(int _rank, bool _flagMPI, bool _flagDebug) {
	if (_flagMPI) {
		++runCount;
		int sleepMS = static_cast<int>(randS->rand() * 1000.0);
		//this_thread::sleep_for(chrono::milliseconds(sleepMS));

		// Nachricht an Master
		int destRank = 0;
		MPIdataBlockSlaveToMaster MPIdata = {};
		MPIdata.rank = _rank;
		MPIdata.command = 1;
		if (cfg.n_steps_char == runCount) {
			MPIdata.command = 2;
		}
		MPI_Send(&MPIdata, 2, MPI_INT, destRank, 0, MPI_COMM_WORLD);
		MPI_Status srstate;
		if (_rank & 1) {
			// ungerade
			sendR(dataRightSend, _flagMPI, _flagDebug);
			//cout << "slave: " << _rank << " count: " << runCount << " ungerade 1 rightComm " << dataRightSend.dArray[0] << endl;
			MPI_Sendrecv(dataRightSend.dArray, dataSize, MPI_DOUBLE, 1, 0, dataRightRecv.dArray, dataSize, MPI_DOUBLE, 1, 0, rightComm, &srstate);

		}
		else {
			// gerade
			sendL(dataLeftSend, _flagMPI, _flagDebug);
			//cout << "slave: " << _rank << " count: " << runCount << " gerade 1 leftComm" << dataLeftSend.dArray[0] << endl;
			MPI_Sendrecv(dataLeftSend.dArray, dataSize, MPI_DOUBLE, 0, 0, dataLeftRecv.dArray, dataSize, MPI_DOUBLE, 0, 0, leftComm, &srstate);
		}
		if (_rank & 1) {
			// ungerade
			sendL(dataLeftSend, _flagMPI, _flagDebug);
			//cout << "slave: " << _rank << " count: " << runCount << " ungerade 2 leftComm" << dataLeftSend.dArray[0] << endl;
			MPI_Sendrecv(dataLeftSend.dArray, dataSize, MPI_DOUBLE, 0, 0, dataLeftRecv.dArray, dataSize, MPI_DOUBLE, 0, 0, leftComm, &srstate);
		}
		else {
			// gerade
			sendR(dataRightSend, _flagMPI, _flagDebug);
			//cout << "slave: " << _rank << " count: " << runCount << " gerade 2 rightComm" << dataRightSend.dArray[0] << endl;
			MPI_Sendrecv(dataRightSend.dArray, dataSize, MPI_DOUBLE, 1, 0, dataRightRecv.dArray, dataSize, MPI_DOUBLE, 1, 0, rightComm, &srstate);
		}
		//cout << "slave: " << _rank << " count: " << runCount << " " << dataLeftRecv.dArray[0] << "<-->" << dataRightRecv.dArray[0] << endl;
		aggregate(dataRightRecv, dataLeftRecv, _flagMPI, _flagDebug);
		MPIdataBlockMasterToSlave MPIdataToSlave = {};
		MPI_Status status;
		MPIdataToSlave.command = 1;
		MPI_Recv(&MPIdataToSlave, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		
	}

	// import arrays
	complex<double> FxLiqOnWall;
	complex<double> dfc_by_n = 0;
	FxLiqOnWall = Stream_LB();
	Compute_drho_u();
	Collide();
	int icount = 0;
	if (icount % cfg.outInterval == 0) {
		complex<double> dfc_by_n_prev = dfc_by_n;
		dfc_by_n = cfg.f0_SI * cfg.IU / M_PI * FxLiqOnWall / cfg.Zq / cfg.ovtorder / static_cast<double>(cfg.nx * cfg.nz);
		if (icount % (1 * cfg.outInterval) == 0) {
			fstream logfile("./results/" + cfg.filename + "_ring_in_" + to_string(_rank) + ".dat", std::ofstream::out | std::ofstream::app);
			logfile << FxLiqOnWall.real() << "\t" << FxLiqOnWall.imag() << endl;
			logfile.close();
		}
		
	}
	icount++;
	if (cfg.n_steps_char == runCount) {
		WritePopulations(_rank);
		cout << "terminate slave " + to_string(_rank) << endl;
		return true;
	}
	return false;
}

void SlaveF::sendR(MPIdataBlockSlaveToSlave _dataRightSend, bool _flagMPI, bool _flagDebug) {
	if (_flagMPI) {

	}
	int iPo = 0;
	for (int x = 0; x < cfg.nx; x++) {
		for (int y = 0; y < cfg.ny; y++) {
			for (int i = 0; i < vst.ndir; i++) {
				_dataRightSend.dArray[iPo] = dyn.h1[x][y][cfg.nz - 1][i].real();
				_dataRightSend.dArray[iPo + 1] = dyn.h1[x][y][cfg.nz - 1][i].imag();
				iPo += 2;
			}
		}
	}

}

void SlaveF::sendL(MPIdataBlockSlaveToSlave _dataLeftSend, bool _flagMPI, bool _flagDebug) {
	if (_flagMPI) {

	}
	int iPo = 0;
	for (int x = 0; x < cfg.nx; x++) {
		for (int y = 0; y < cfg.ny; y++) {
			for (int i = 0; i < vst.ndir; i++) {
				_dataLeftSend.dArray[iPo] = dyn.h1[x][y][0][i].real();
				_dataLeftSend.dArray[iPo + 1] = dyn.h1[x][y][0][i].imag();
				iPo += 2;
			}
		}
	}
}

void SlaveF::aggregate(MPIdataBlockSlaveToSlave _dataRightRecv, MPIdataBlockSlaveToSlave _dataLeftRecv, bool _flagMPI, bool _flagDebug) {
	if (_flagMPI) {

	}
	else {
		int iPo = 0;
		for (int x = 0; x < cfg.nx; x++) {
			for (int y = 0; y < cfg.ny; y++) {
				for (int i = 0; i < vst.ndir; i++) {
					dyn.h1[x][y][0][i] = _dataRightRecv.dArray[iPo] + cfg.IU * _dataRightRecv.dArray[iPo + 1];
					dyn.h1[x][y][cfg.nz - 1][i] = _dataLeftRecv.dArray[iPo] + cfg.IU * _dataLeftRecv.dArray[iPo + 1];
					iPo += 2;
				}
			}
		}
	}
}

void SlaveF::Calc_Parameters_from_Config() {

	cfg.dq = cfg.dq_SI / cfg.gridreso_SI;
	cfg.ovtorder = 2. * cfg.iovt + 1;
	cfg.delta_SI = pow((2 * cfg.nu_SI / (2 * M_PI * cfg.f0_SI * cfg.ovtorder)), 0.5);
	cfg.ny = static_cast<int>(cfg.cellheight_SI / cfg.gridreso_SI);
	cfg.nu = 1.0 / 6.0;
	if (cfg.tandelinv == 0) { cfg.nu_Sph = cfg.nu * cfg.nufac; }
	else {
		double tandel = 1 / cfg.tandelinv;  //tandel = inv(tandelinv)
		double cos = 1 / pow((1 + tandel * tandel), 0.5);
		double sin = tandel / pow((1 + tandel * tandel), 0.5);
		cfg.nu_Sph = cfg.nu * cfg.nufac * (sin - cfg.IU * cos);
	}
	if (cfg.B_Reference) { cfg.nu_Sph = cfg.nu; }
	else {
		cfg.nu_Sph = cfg.nu * cfg.nufac * (1.0 - cfg.IU * cfg.tandelinv);
	}
	cfg.nx = static_cast<int>(2 * cfg.rad_Sph_SI / cfg.gridreso_SI * cfg.cell_diam_by_sphere_diam);
	cfg.nz = cfg.nx;
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
	double cs = pow((1. / 3.), 0.5);
	cfg.ZComp = cfg.rho0 * cs;
	cfg.Zliq_Sph = pow((cfg.IU * cfg.rho0 * cfg.omega * cfg.rho0 * cfg.nu_Sph), 0.5);
	cfg.f0 = cfg.omega / (2.0 * M_PI * cfg.ovtorder);
	cfg.dq = cfg.dq_SI / cfg.gridreso_SI;
	cfg.Zq = 2.0 * cfg.dq * cfg.f0 * cfg.rhoq;
	cfg.SphR = cfg.rad_Sph_SI / cfg.gridreso_SI;
//	cfg.Sphx = static_cast<double>(cfg.nx) / 5. - 0.5; ////
//	cfg.Sphz = static_cast<double>(cfg.nz) / 2. - 0.5;
	cfg.Sphy = 2 * cfg.SphR * (cfg.frac_above_surf - 0.5);
	cfg.n_steps_char = int(cfg.ny * cfg.ny / cfg.nu * cfg.char_time_fac);
	cout << cfg.Sphx0;
	cout << cfg.Sphx1;
	cout << cfg.Sphx2;
}
void SlaveF::Set_Config() {

	cfg.filename = "test_000";
	cfg.char_time_fac = 2;             //factor of ny^2 wrt nu
	cfg.cell_diam_by_sphere_diam = 1.3;
	cfg.gridreso_SI = 1.0e-9;
	cfg.cellheight_SI = 8e-9;
	cfg.frac_above_surf = 0.5;
	cfg.rad_Sph_SI = 5e-9;
	cfg.nufac = 1000;                 //very very very hard sphere // overwritten by the loop
	cfg.tandelinv = 10.0;            //overwritten
	cfg.iovt = 0;
	cfg.outInterval = 20;
	cfg.f0_SI = 5e6;
	cfg.Zq_SI = 8.8e6;
	cfg.dq_SI = 332.075e-6;
	cfg.rhoq = 2.65;
	cfg.nu_SI = 1.0e-6;
	cfg.nu = 1.0 / 6.0;
	cfg.rho0 = 1.0;
}

void SlaveF::Set_Config_File_Input() {
	ifstream testOpenFile("config.dat");
	if (!testOpenFile)
	{ 
		cout << "error: config file " << "config.dat" << " not found" << endl;
		exit(1);
	}
	cfg.char_time_fac = 0;
	cfg.rad_Sph_SI = 0.0;
	cfg.cellheight_SI = 0.0;
	cfg.gridreso_SI = 0.0;
	cfg.cell_diam_by_sphere_diam = 0.0;
	cfg.frac_above_surf = 0.0;
	cfg.nufac = 0;
	cfg.tandelinv = 0;
	cfg.iovt = 0;
	cfg.outInterval = 0;
	ifstream is("config.dat", ios::in);
	string line;
	int parameterSet = 0;
	while (getline(is, line)) {
		char rmChars[] = " \t";
		// entferne Leerzeichen und tabs
		for (unsigned int i = 0; i < 2; i++)
		{
			line.erase(std::remove(line.begin(), line.end(), rmChars[i]), line.end());
		}
		string delimiter = "=";
		size_t strpos = 0;
		
		string teilstring;
		vector<string> splitstring;
		while ((strpos = line.find(delimiter)) != string::npos) {
			string teilstring = line.substr(0, strpos);
			splitstring.push_back(teilstring);
			line.erase(0, strpos + delimiter.length());
			splitstring.push_back(line);
			break;
		}

		if (splitstring.size() < 2) {
			continue;
		}
		if (splitstring.at(0) == "char_time_fac") {
			cfg.char_time_fac = stoi(splitstring.at(1));
			++parameterSet;
		}
		if (splitstring.at(0) == "rad_Sph_SI") {
			cfg.rad_Sph_SI = stod(splitstring.at(1));
			++parameterSet;
		}
		if (splitstring.at(0) == "cellheight_SI") {
			cfg.cellheight_SI = stod(splitstring.at(1));
			++parameterSet;
		}
		if (splitstring.at(0) == "gridreso_SI") {
			cfg.gridreso_SI = stod(splitstring.at(1));
			++parameterSet;
		}
		if (splitstring.at(0) == "cell_diam_by_sphere_diam") {
			cfg.cell_diam_by_sphere_diam = stod(splitstring.at(1));
			++parameterSet;
		}
		if (splitstring.at(0) == "frac_above_surf") {
			cfg.frac_above_surf = stod(splitstring.at(1));
			++parameterSet;
		}
		if (splitstring.at(0) == "nufac") {
			cfg.nufac = stoi(splitstring.at(1));
			++parameterSet;
		}
		if (splitstring.at(0) == "tandelinv") {
			cfg.tandelinv = stoi(splitstring.at(1));
			++parameterSet;
		}
		if (splitstring.at(0) == "Sphx0") {
			cfg.Sphx0 = stoi(splitstring.at(1));
			++parameterSet;
		}
		if (splitstring.at(0) == "Sphz0") {
			cfg.Sphz0 = stoi(splitstring.at(1));
			++parameterSet;
		}
		if (splitstring.at(0) == "Sphx1") {
			cfg.Sphx1 = stoi(splitstring.at(1));
			++parameterSet;
		}
		if (splitstring.at(0) == "Sphz1") {
			cfg.Sphz1 = stoi(splitstring.at(1));
			++parameterSet;
		}
		if (splitstring.at(0) == "Sphx2") {
			cfg.Sphx2 = stoi(splitstring.at(1));
			++parameterSet;
		}
		if (splitstring.at(0) == "Sphz2") {
			cfg.Sphz2 = stoi(splitstring.at(1));
			++parameterSet;
		}
		if (splitstring.at(0) == "iovt") {
			cfg.iovt = stoi(splitstring.at(1));
			++parameterSet;
		}
		if (splitstring.at(0) == "outInterval") {
			cfg.outInterval = stoi(splitstring.at(1));
			++parameterSet;
		}
	}
	if (parameterSet != 10) {
		cout << "config file incomplete" << endl;
		exit(1);
	}
}
void SlaveF::Count_and_Allocate() {
	dyn.dr.resize(cfg.nx);
	dyn.ux.resize(cfg.nx);
	dyn.uy.resize(cfg.nx);
	dyn.uz.resize(cfg.nx);

	dyn.one_m_tauinv_m_Iom_arr.resize(cfg.nx);
	dyn.tauinv_m_om2_by_2_arr.resize(cfg.nx);
	dyn.one_m_tauinv_p_m_Iom_arr.resize(cfg.nx);
	dyn.tauinv_p_m_om2_by_2_arr.resize(cfg.nx);
	dyn.one_m_tauinv_m_m_Iom_arr.resize(cfg.nx);
	dyn.tauinv_m_m_om2_by_2_arr.resize(cfg.nx);
	dyn.h1.resize(cfg.nx);
	dyn.h2.resize(cfg.nx);

	for (int x = 0; x < cfg.nx; x++) {
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


void SlaveF::Free_Mem() {
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
	}
	dyn.h1.resize(0);
	dyn.h2.resize(0);
	dyn.dr.resize(0);
	dyn.ux.resize(0);
	dyn.uy.resize(0);
	dyn.uz.resize(0);
	dyn.one_m_tauinv_m_Iom_arr.resize(0);
	dyn.tauinv_m_om2_by_2_arr.resize(0);
}
void SlaveF::Set_VE_Pars(int _rank, int _nSlaves) {
	for (size_t x = 0; x < cfg.nx; x++) {
		for (size_t y = 0; y < cfg.ny; y++) {
			for (size_t z = 0; z < cfg.nz; z++) {
				int z_shift = (_rank - 1 - _nSlaves / 2) * cfg.nz + z;
				dyn.one_m_tauinv_m_Iom_arr[x][y][z] = cfg.one_m_tauinv_m_Iom;
				dyn.tauinv_m_om2_by_2_arr[x][y][z] = cfg.tauinv_m_om2b2;

				double xd = static_cast<double>(x);
				double zd = static_cast<double>(z_shift);
				double nxd = static_cast<double>(cfg.nx);
				double nzd = static_cast<double>(cfg.nz);
				double r0 = pow(pow(xd - cfg.Sphx0, 2) + pow(y - cfg.Sphy, 2) + pow(zd - cfg.Sphz0, 2), 0.5);
				double rpx = pow(pow(xd - (cfg.Sphx0 + nxd), 2) + pow(y - cfg.Sphy, 2) + pow(zd - cfg.Sphz0, 2), 0.5);
				double rmx = pow(pow(xd - (cfg.Sphx0 - nxd), 2) + pow(y - cfg.Sphy, 2) + pow(zd - cfg.Sphz0, 2), 0.5);
				double rpz = pow(pow(xd - cfg.Sphx0, 2) + pow(y - cfg.Sphy, 2) + pow(zd - (cfg.Sphz0 + nzd), 2), 0.5);
				double rmz = pow(pow(xd - cfg.Sphx0, 2) + pow(y - cfg.Sphy, 2) + pow(zd - (cfg.Sphz0 - nzd), 2), 0.5);

				if ((r0 < cfg.SphR) || (rpx < cfg.SphR) || (rmx < cfg.SphR) || (rpz < cfg.SphR) || (rmz < cfg.SphR)) {
					dyn.one_m_tauinv_m_Iom_arr[x][y][z] = cfg.one_m_tauinv_m_Iom_Sph;
					dyn.tauinv_m_om2_by_2_arr[x][y][z] = cfg.tauinv_m_om2b2_Sph;
				}
				r0 = pow(pow(xd - cfg.Sphx1, 2) + pow(y - cfg.Sphy, 2) + pow(zd - cfg.Sphz1, 2), 0.5);
				rpx = pow(pow(xd - (cfg.Sphx1 + nxd), 2) + pow(y - cfg.Sphy, 2) + pow(zd - cfg.Sphz1, 2), 0.5);
				rmx = pow(pow(xd - (cfg.Sphx1 - nxd), 2) + pow(y - cfg.Sphy, 2) + pow(zd - cfg.Sphz1, 2), 0.5);
				rpz = pow(pow(xd - cfg.Sphx1, 2) + pow(y - cfg.Sphy, 2) + pow(zd - (cfg.Sphz1 + nzd), 2), 0.5);
				rmz = pow(pow(xd - cfg.Sphx1, 2) + pow(y - cfg.Sphy, 2) + pow(zd - (cfg.Sphz1 - nzd), 2), 0.5);

				if ((r0 < cfg.SphR) || (rpx < cfg.SphR) || (rmx < cfg.SphR) || (rpz < cfg.SphR) || (rmz < cfg.SphR)) {
					dyn.one_m_tauinv_m_Iom_arr[x][y][z] = cfg.one_m_tauinv_m_Iom_Sph;
					dyn.tauinv_m_om2_by_2_arr[x][y][z] = cfg.tauinv_m_om2b2_Sph;
				}
				r0 = pow(pow(xd - cfg.Sphx1, 2) + pow(y - cfg.Sphy, 2) + pow(zd - cfg.Sphz2, 2), 0.5);
				rpx = pow(pow(xd - (cfg.Sphx2 + nxd), 2) + pow(y - cfg.Sphy, 2) + pow(zd - cfg.Sphz2, 2), 0.5);
				rmx = pow(pow(xd - (cfg.Sphx2 - nxd), 2) + pow(y - cfg.Sphy, 2) + pow(zd - cfg.Sphz2, 2), 0.5);
				rpz = pow(pow(xd - cfg.Sphx2, 2) + pow(y - cfg.Sphy, 2) + pow(zd - (cfg.Sphz2 + nzd), 2), 0.5);
				rmz = pow(pow(xd - cfg.Sphx2, 2) + pow(y - cfg.Sphy, 2) + pow(zd - (cfg.Sphz2 - nzd), 2), 0.5);

				if ((r0 < cfg.SphR) || (rpx < cfg.SphR) || (rmx < cfg.SphR) || (rpz < cfg.SphR) || (rmz < cfg.SphR)) {
					dyn.one_m_tauinv_m_Iom_arr[x][y][z] = cfg.one_m_tauinv_m_Iom_Sph;
					dyn.tauinv_m_om2_by_2_arr[x][y][z] = cfg.tauinv_m_om2b2_Sph;
				}

			}
		}
	}
}

void SlaveF::init_h() {
	for (size_t x = 0; x < cfg.nx; x++) {
		for (size_t y = 0; y < cfg.ny; y++) {
			for (size_t z = 0; z < cfg.nz; z++) {
				for (size_t i = 0; i < vst.ndir; i++) { dyn.h1[x][y][z][i] = 0; }
			}
		}
	}
}

complex<double> SlaveF::Stream_LB() {
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
						complex<double> uxTop = 2.0 * (dyn.h1[x][y][z][7] - dyn.h1[x][y][z][14]) / (cfg.Zliq + 1.0 / 3.0);
						complex<double> uyTop = 2.0 * (dyn.h1[x][y][z][11] + dyn.h1[x][y][z][17] +
							dyn.h1[x][y][z][7] + dyn.h1[x][y][z][14] +
							dyn.h1[x][y][z][3]) / (cfg.ZComp + 1.0);
						complex<double> uzTop = 2.0 * (dyn.h1[x][y][z][11] - dyn.h1[x][y][z][17]) / (cfg.Zliq + 1.0 / 3.0);
						dyn.h2[x][y][z][i] = dyn.h1[x][y][z][ibar] - 6.0 * vst.wi[ibar] * vst.dirx[ibar] * uxTop
							- 6.0 * vst.wi[ibar] * vst.diry[ibar] * uyTop
							- 6.0 * vst.wi[ibar] * vst.dirz[ibar] * uzTop;
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

void SlaveF::Compute_drho_u() {
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
void SlaveF::Collide() {
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
void SlaveF::WritePopulations(int _rank) {
	ofstream outfile;
	outfile.open("./results/" + cfg.filename + "_Pop_" + to_string(_rank) + ".dat", std::ofstream::out);
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

int SlaveF::getDimensionSendRecv() {
	return dimensionSendRecv;
}