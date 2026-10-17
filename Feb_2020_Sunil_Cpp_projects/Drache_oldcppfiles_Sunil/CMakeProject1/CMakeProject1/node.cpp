#include "velocitySets.h"
#include "node.h"
#include "fdLBMCmdl.h"

void Node::setPosition(const NodePosition& _position) {
	// register position
	position.id = _position.id;
	position.x = _position.x;
	position.y = _position.y;
	position.z = _position.z;
}
void Node::initNext(const size_t _countNext) {
	nextNode = new Node*[_countNext];
	for (size_t i = 0; i < _countNext; i++) {
		nextNode[i] = nullptr;
	}
}
Node::Node(const D2Q9& _velocitySet, const NodePosition& _position) {
	setPosition(_position);
	initNext(_velocitySet.numDirections);
}
Node::Node(const D3Q15& _velocitySet, const NodePosition& _position) {
	setPosition(_position);
	initNext(_velocitySet.numDirections);
}
Node::Node(const D3Q19& _velocitySet, const NodePosition& _position) {
	h1 = new complex<double>[_velocitySet.numDirections];
	h2 = new complex<double>[_velocitySet.numDirections];
	ciX = new int[_velocitySet.numDirections];
	ciY = new int[_velocitySet.numDirections];
	ciZ = new int[_velocitySet.numDirections];
	for (size_t i = 0; i < _velocitySet.numDirections; i++) {
		h1[i] = 0.0;
		h2[i] = 0.0;
		ciX[i] = _velocitySet.ciX[i];
		ciY[i] = _velocitySet.ciY[i];
		ciZ[i] = _velocitySet.ciZ[i];
	}
	setPosition(_position);
	initNext(_velocitySet.numDirections);
	numDirections = _velocitySet.numDirections;
}
Node::Node(const D3Q27& _velocitySet, const NodePosition& _position) {
	setPosition(_position);
	initNext(_velocitySet.numDirections);
}
void Node::connect(Node *_nextNode, const size_t _position) {
	nextNode[_position] = _nextNode;
}
NodePosition Node::getPosition() {
	return position;
}
void Node::compute_dr_u() {
	dr = (0.0, 0.0);
	ux = (0.0, 0.0);
	uy = (0.0, 0.0);
	uz = (0.0, 0.0);
	for (size_t i = 0; i < numDirections; i++) {
		dr += h2[i];
		ux += (h2[i] * static_cast<double>(ciX[i]));
		uy += (h2[i] * static_cast<double>(ciY[i]));
		uz += (h2[i] * static_cast<double>(ciZ[i]));
	}
}
void Node::collide(const modelParameter& _modelParameter) {
	for (size_t i = 0; i < numDirections; i++) {
		complex<double> cidotu = static_cast<double>(ciX[i]) * ux + static_cast<double>(ciY[i]) * uy + static_cast<double>(ciZ[i]) *uz;
		complex<double> heq = wi[i] * (dr + 3.0 * cidotu);
		complex<double> c01 = (0.0, 1.0);
		h1[i] = _modelParameter.tauinv * heq + c01 * _modelParameter.omega * h2[i] - _modelParameter.omega2by2 * heq;
	}
}
void Node::LB_Stream(Fx& _Fx) {

}