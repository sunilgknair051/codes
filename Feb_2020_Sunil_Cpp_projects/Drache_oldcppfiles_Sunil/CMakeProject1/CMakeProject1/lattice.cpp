#include <vector>
#include <string>
#include <iostream>
#include "velocitySets.h"
#include "node.h"
#include "lattice.h"

using namespace std;
size_t Lattice::nextNodeLattice(const int _xPos, const  int _yPos, const  int _zPos, const LBParameter& _lbParameter) {
	int xPos = _xPos;
	int yPos = _yPos;
	int zPos = _zPos;
	if (_yPos < 0) {
		yPos += _lbParameter.dimensionY;
	}
	if (yPos >= _lbParameter.dimensionY) {
		yPos -= _lbParameter.dimensionY;
	}
	if (xPos < 0) {
		xPos += _lbParameter.dimensionX;
	}
	if (xPos >= _lbParameter.dimensionX) {
		xPos -= _lbParameter.dimensionX;
	}
	if (zPos < 0) {
		if (_lbParameter.periodicX) {
			zPos += _lbParameter.dimensionZ;
		}
		else {
			zPos = 0;
		}
	}
	if (zPos >= _lbParameter.dimensionZ) {
		if (_lbParameter.periodicX) {
			zPos -= _lbParameter.dimensionZ;
		}
		else {
			zPos = _lbParameter.dimensionZ - 1;
		}
	}
	if (_lbParameter.Lattice3D) {
		return zPos + yPos *_lbParameter.dimensionZ + xPos * _lbParameter.dimensionZ * _lbParameter.dimensionY;
	}
	else {
		return yPos + xPos * _lbParameter.dimensionX;
	}
}
Lattice::Lattice(const LBParameter& _lbParameter) {
	// init Lattice
	if (_lbParameter.velocitySet.compare("D2Q9") == 0) {
		D2Q9 velocitySet;
		size_t id = 0;
		for (size_t i = 0; i < _lbParameter.dimensionY; i++) {
			for (size_t j = 0; j < _lbParameter.dimensionX; j++) {
				NodePosition pos;
				pos.id = id;
				pos.x = j;
				pos.y = i;
				Node *tmpNode = new Node(velocitySet, pos);
				latticeNodes.push_back(tmpNode);
				++id;
			}
		}
		// set connects 
		id = 0;
		for (size_t i = 0; i < _lbParameter.dimensionX; i++) {
			for (size_t j = 0; j < _lbParameter.dimensionY; j++) {
				for (size_t dir = 0; dir < velocitySet.numDirections; dir++) {
					int yPos = static_cast<int>(j) + velocitySet.ciY[dir];
					int xPos = static_cast<int>(i) + velocitySet.ciX[dir];
					size_t nextNode = nextNodeLattice(xPos, yPos, 0, _lbParameter);
					latticeNodes.at(id)->connect(latticeNodes.at(nextNode), dir);
					std::cout << "id: " << id << " dir: " << dir << " " << nextNode << endl;
				}
				++id;
			}
		}
	}
	if (_lbParameter.velocitySet.compare("D3Q15") == 0) {
		D3Q15 velocitySet;
		size_t id = 0;
		for (size_t i = 0; i < _lbParameter.dimensionX; i++) {
			for (size_t j = 0; j < _lbParameter.dimensionY; j++) {
				for (size_t k = 0; k < _lbParameter.dimensionZ; k++) {
					NodePosition pos;
					pos.id = id;
					pos.x = i;
					pos.y = j;
					pos.z = k;
					Node *tmpNode = new Node(velocitySet, pos);
					latticeNodes.push_back(tmpNode);
					++id;
					// connect layer 0 an 1 at position 0
					//layer0.back()->connect(layer1.back(), 0);
					//layer1.back()->connect(layer0.back(), 0);
				}
			}
		}
	}
	if (_lbParameter.velocitySet.compare("D3Q19") == 0) {
		D3Q19 velocitySet;
		size_t id = 0;
		for (size_t i = 0; i < _lbParameter.dimensionX; i++) {
			for (size_t j = 0; j < _lbParameter.dimensionY; j++) {
				for (size_t k = 0; k < _lbParameter.dimensionZ; k++) {
					NodePosition pos;
					pos.id = id;
					pos.x = i;
					pos.y = j;
					pos.z = k;
					Node *tmpNode = new Node(velocitySet, pos);
					latticeNodes.push_back(tmpNode);
					++id;
				}
			}
		}
		id = 0;
		for (size_t i = 0; i < _lbParameter.dimensionX; i++) {
			for (size_t j = 0; j < _lbParameter.dimensionY; j++) {
				for (size_t k = 0; k < _lbParameter.dimensionZ; k++) {
					for (size_t dir = 0; dir < velocitySet.numDirections; dir++) {
						int zPos = static_cast<int>(k) + velocitySet.ciZ[dir];
						int yPos = static_cast<int>(j) + velocitySet.ciY[dir];
						int xPos = static_cast<int>(i) + velocitySet.ciX[dir];

						size_t nextNode = nextNodeLattice(xPos, yPos, zPos, _lbParameter);
						latticeNodes.at(id)->connect(latticeNodes.at(nextNode), dir);
						//cout << "vX: " << velocitySet.ciX[dir] << " vY: " << velocitySet.ciY[dir] << " vZ: " << velocitySet.ciZ[dir] << endl;
						//cout << "id: " << id << " dir: " << dir << " " << nextNode << endl;
					}
					++id;
				}
			}
		}
	}

	if (_lbParameter.velocitySet.compare("D3Q27") == 0) {
		D3Q27 velocitySet;
		size_t id = 0;
		for (size_t i = 0; i < _lbParameter.dimensionX; i++) {
			for (size_t j = 0; j < _lbParameter.dimensionY; j++) {
				for (size_t k = 0; k < _lbParameter.dimensionZ; k++) {
					NodePosition pos;
					pos.id = id;
					pos.x = i;
					pos.y = j;
					pos.z = k;
					/*
					pos.layer = 0;
					Node *tmpNode0 = new Node(velocitySet, pos);
					layer0.push_back(tmpNode0);
					pos.layer = 1;
					Node *tmpNode1 = new Node(velocitySet, pos);
					layer1.push_back(tmpNode0);
					++id;
					// connect layer 0 an 1 at position 0
					layer0.back()->connect(layer1.back(), 0);
					layer1.back()->connect(layer0.back(), 0);
					*/
				}
			}
		}
	}
}
void Lattice::LB_Stream(Fx& _Fx) {
	_Fx.LiqOnSphe = 0.0;
	_Fx.LiqOnWall = 0.0;
	for (auto& elem : latticeNodes) {
		elem->LB_Stream(_Fx);
	}
}
void Lattice::collide(const modelParameter& _modelParameter) {
	for (auto& elem : latticeNodes) {
		elem->collide(_modelParameter);
	}
}

complex<double> Lattice::calc_dfc(const LBParameter& _LBParameter, const modelParameter& _modelParameter, Fx& _Fx) {
	complex<double> TanStress = (_Fx.LiqOnWall + _Fx.LiqOnSphe) / (_LBParameter.dimensionX * _LBParameter.dimensionZ);
	complex<double> ZL = TanStress / _modelParameter.uxWall;
	complex<double> Zq = (_modelParameter.rho0 / _modelParameter.rho_Liq_SI)*(_modelParameter.nu / _modelParameter.nu_SI)*(_modelParameter.delta_f0_SI / _modelParameter.delta_f0)*_modelParameter.Zq_SI;
	complex<double> c01 = (0.0, 1.0);
	complex<double> dfc = _modelParameter.f0_SI * c01 / M_PI * ZL / Zq;
	dfc = -conj(dfc);
	return dfc;
}