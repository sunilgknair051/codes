#pragma once
#include <complex>
#include <cmath>
#include "velocitySets.h"
#include "fdLBMCmdl.h"

using namespace std;
using namespace complex_literals;

struct NodePosition {
	// debug information
	size_t id = 0;
	size_t x = 0;
	size_t y = 0;
	size_t z = 0;
};

class Node
{
public:
	Node(const D2Q9& _velocitySet, const NodePosition& _position);
	Node(const D3Q15& _velocitySet, const NodePosition& _position);
	Node(const D3Q19& _velocitySet, const NodePosition& _position);
	Node(const D3Q27& _velocitySet, const NodePosition& _position);

	void connect(Node *_nextNode0, const size_t _position);
	NodePosition getPosition();
	void compute_dr_u();
	void collide(const modelParameter& _modelParameter);
	void LB_Stream(Fx& _Fx);

private:
	NodePosition position;
	void setPosition(const NodePosition& _position);
	void initNext(const size_t _countNext);
	Node **nextNode;
	size_t numDirections;

	// dimension numDirections
	int *ciX;
	int *ciY;
	int *ciZ;
	double *wi;
	complex<double> *h1;
	complex<double> *h2;

	complex<double> dr;
	complex<double> ux;
	complex<double> uy;
	complex<double> uz;
};
