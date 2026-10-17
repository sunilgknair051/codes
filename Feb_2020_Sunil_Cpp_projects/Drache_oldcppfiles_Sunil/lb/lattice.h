#pragma once
#include <vector>
#include "fdLBMCmdl.h"
#include "velocitySets.h"
#include "node.h"

class Lattice
{
public:
	Lattice(const LBParameter& _lbParameter);
	void LB_Stream(Fx& _Fx);
	void collide(const modelParameter& _modelParameter);
	complex<double> calc_dfc(const LBParameter& _LBParameter, const modelParameter& _modelParameter, Fx& _Fx);
private:
	vector<Node*> latticeNodes;
	size_t nextNodeLattice(const int _xPos, const  int _yPos, const  int _zPos, const LBParameter& _lbParameter);
};