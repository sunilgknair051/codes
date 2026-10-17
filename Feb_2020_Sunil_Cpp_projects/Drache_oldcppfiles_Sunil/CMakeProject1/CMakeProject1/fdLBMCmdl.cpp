#include <string>
#include <fstream>
#include "cxxopts.hpp"
#include "parsing.h"
#include "fdLBMCmdl.h"
#include "velocitySets.h"
#include "node.h"
#include "lattice.h"

using namespace std;

bool filenameValid(const string& _filename) {
	if (_filename.size() < 1) {
		cout << "filename invalid" << endl;
		return false;
	}
	ifstream is;
	is.open(_filename, ifstream::in);
	if (!is.is_open()) {
		cout << "error open file " << _filename << endl;
		return false;
	}
	char c = is.get();
	if (!is.good()) {
		cout << "error reading file " << _filename << endl;
		return false;
	}
	return true;
}

LBParameter getParameterFromModelFile(const string& _filename) {
	LBParameter simulationParameter;
	ifstream is(_filename, ios::in);
	if (!is.is_open())
		return simulationParameter;
	string line;
	cout << "model file: " << _filename << endl;
	while (getline(is, line)) {
		if (line.find("velocitySet", 0) != string::npos) {
			vector<string> specLine = parsing::splitline(line);
			if (specLine.size() < 2)
				continue;
			string vSets ="D2Q9 D3Q15 D3Q19 D3Q27";
			if (vSets.compare(specLine.at(1)) != 0) {
				simulationParameter.velocitySet = specLine.at(1);
			}
			continue;
		}
		if (line.find("baseDimension", 0) != string::npos) {
			vector<string> specLine = parsing::splitline(line);
			if (specLine.size() < 2)
				continue;
			simulationParameter.baseDimension = stol(specLine.at(1));
			continue;
		}
		if (line.find("dimensionX", 0) != string::npos) {
			vector<string> specLine = parsing::splitline(line);
			if (specLine.size() < 2)
				continue;
			simulationParameter.dimensionX = stol(specLine.at(1));
			continue;
		}
		if (line.find("dimensionY", 0) != string::npos) {
			vector<string> specLine = parsing::splitline(line);
			if (specLine.size() < 2)
				continue;
			simulationParameter.dimensionY = stol(specLine.at(1));
			continue;
		}
		if (line.find("dimensionZ", 0) != string::npos) {
			vector<string> specLine = parsing::splitline(line);
			if (specLine.size() < 2)
				continue;
			simulationParameter.dimensionZ = stol(specLine.at(1));
			continue;
		}
		if (line.find("bounceBackZ", 0) != string::npos) {
			vector<string> specLine = parsing::splitline(line);
			if (specLine.size() < 2)
				continue;
			simulationParameter.bounceBackZ = specLine.at(1);
			simulationParameter.periodicZ = false;
			continue;
		}
		if (line.find("fluidDensity", 0) != string::npos) {
			vector<string> specLine = parsing::splitline(line);
			if (specLine.size() < 2)
				continue;
			simulationParameter.fluidDensity = stod(specLine.at(1));
			continue;
		}
	}

	return simulationParameter;
}

int main(int argc, char* argv[])
{
	LBParameter simulationParameter;
	if (argc == 1) {
		cout << "Lattice Boltzmann simulation program, Version 0.1" << endl;
		cout << "Copyright (C) 2018 Marco Drache, Sunil Gopalakrishna, Diethelm Johannsmann" << endl;
		cout << "Usage:" << endl;
		cout << "  lb - m modelFile [OPTION..]" << endl;
		cout << "  lb --help" << endl;
		exit(0);
	}
	try
	{
		cxxopts::Options options(argv[0], "Example:\n  lb -m basicParameter.lb --scale=10 ");
		options.positional_help("[optional args]");

		options.add_options()
			("m, model", "Model", cxxopts::value<std::vector<std::string>>(), "modelFILE")
			("help", "Print help")
			("license", "Print license")
			("scale", "relative model dimension", cxxopts::value<uint64_t>(), "unsigned int")
			;

		options.parse(argc, argv);

		if (options.count("help"))
		{
			std::cout << options.help({ "", "Group" }) << std::endl;
			exit(0);
		}

		if (options.count("m"))
		{
			auto& fileVec = options["m"].as<std::vector<std::string>>();
			if (!filenameValid(fileVec.at(0))) {
				exit(1);
			}
			simulationParameter.modelFile = fileVec.at(0);
		}
		if (options.count("steps"))
		{
			simulationParameter.steps = options["steps"].as<unsigned int>();
		}
		if (options.count("license"))
		{
			cout << "Lattice Boltzmann simulation program" << endl;
			cout << "Copyright (C) 2018 Diethelm Johannsmann, Sunil Gopalakrishna, Marco Drache" << endl;
			cout << endl;
			cout << "This program is free software: you can redistribute it and/or modify" << endl;
			cout << "it under the terms of the GNU General Public License as published by" << endl;
			cout << "the Free Software Foundation, either version 3 of the License, or" << endl;
			cout << "(at your option) any later version." << endl;
			cout << endl;
			cout << "This program is distributed in the hope that it will be useful," << endl;
			cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of" << endl;
			cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << endl;
			cout << "GNU General Public License for more details." << endl;
			cout << endl;
			cout << "You should have received a copy of the GNU General Public License" << endl;
			cout << "along with this program.  If not, see <http://www.gnu.org/licenses/>." << endl;
			exit(0);
		}
		//std::cout << "Arguments remain = " << argc << "\t" << argv[0] << std::endl;
	
	}
	

	catch (const cxxopts::OptionException& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}
	
	auto simulationParameterSet = getParameterFromModelFile(simulationParameter.modelFile);
	simulationParameterSet.init();
	cout << "model dimension: " << simulationParameterSet.getDimension() << endl;

	Lattice *lbLattice = new Lattice(simulationParameterSet);
	Fx FX;
	modelParameter MP;
	for (size_t i = 0; i < simulationParameter.steps; i++) {
		lbLattice->LB_Stream(FX);
		lbLattice->collide(MP);
		complex<double> dfc = lbLattice->calc_dfc(simulationParameterSet, MP, FX);
		cout << "step: " << i << " " << dfc << endl;
	}
	return 0;

}
