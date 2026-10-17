#pragma once
#include <math.h>
struct D2Q9 {
	const unsigned int numDirections = 9;
	// lattice weights and lenght
	const double zeroWeight = 4.0 / 9.0;
	const double adjacentWeight = 1.0 / 9.0;
	const double adjacentLenght = 1.0;
	const double diagonalWeight = 1.0 / 36.0;
	const double diagonalLenght = sqrt(2.0);
	const int ciX[9] = { 0,1,0,-1,0,1,-1,-1,1 };
	const int ciY[9] = { 0,0,1,0,-1,1,1,-1,-1 };
	const double wi[9] = { zeroWeight,adjacentWeight,adjacentWeight,adjacentWeight,adjacentWeight,diagonalWeight,diagonalWeight,diagonalWeight,diagonalWeight };
	const double length[9] = { 0.0, adjacentLenght,adjacentLenght,adjacentLenght,adjacentLenght,diagonalLenght,diagonalLenght,diagonalLenght,diagonalLenght };
};
struct D3Q15 {
	const unsigned int numDirections = 15;
	// lattice weights and lenght
	const double zeroWeight = 2.0 / 9.0;
	const double adjacentWeight = 1.0 / 9.0;
	const double adjacentLenght = 1.0;
	const double cornerWeight = 1.0 / 72.0;
	const double cornerLenght = sqrt(3.0);
	const int ciX[15] = { 0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,-1,1 };
	const int ciY[15] = { 0,0,0,1,-1,0,0,1,-1,1,-1,-1,1,1,-1 };
	const int ciZ[15] = { 0,0,0,0,0,1,-1,1,-1,-1,1,1,-1,1,-1 };
	const double wi[15] = { zeroWeight,adjacentWeight,adjacentWeight,adjacentWeight,adjacentWeight,adjacentWeight,adjacentWeight,cornerWeight,cornerWeight,cornerWeight,cornerWeight,cornerWeight,cornerWeight,cornerWeight,cornerWeight };
	const double length[15] = { 0.0, adjacentLenght,adjacentLenght,adjacentLenght,adjacentLenght,adjacentLenght,adjacentLenght,cornerLenght,cornerLenght,cornerLenght,cornerLenght,cornerLenght,cornerLenght,cornerLenght,cornerLenght };
};
struct D3Q19 {
	const unsigned int numDirections = 19;
	// lattice weights and lenght
	const double zeroWeight = 1.0 / 3.0;
	const double adjacentWeight = 1.0 / 18.0;
	const double adjacentLenght = 1.0;
	const double diagonalWeight = 1.0 / 36.0;
	const double diagonalLenght = sqrt(2.0);
	const int ciX[19] = { 0,1,-1,0,0,0,0,1,-1,1,-1,0,0,1,-1,1,-1,0,0 };
	const int ciY[19] = { 0,0,0,1,-1,0,0,1,-1,0,0,1,-1,-1,1,0,0,1,-1 };
	const int ciZ[19] = { 0,0,0,0,0,1,-1,0,0,1,-1,1,-1,0,0,-1,1,-1,1 };
	const double wi[19] = { zeroWeight,adjacentWeight,adjacentWeight,adjacentWeight,adjacentWeight,adjacentWeight,adjacentWeight,diagonalWeight,diagonalWeight,diagonalWeight,diagonalWeight,diagonalWeight,diagonalWeight,diagonalWeight,diagonalWeight,diagonalWeight,diagonalWeight,diagonalWeight,diagonalWeight };
	const double length[19] = { 0.0, adjacentLenght,adjacentLenght,adjacentLenght,adjacentLenght,adjacentLenght,adjacentLenght,diagonalLenght,diagonalLenght,diagonalLenght,diagonalLenght,diagonalLenght,diagonalLenght,diagonalLenght,diagonalLenght,diagonalLenght,diagonalLenght,diagonalLenght,diagonalLenght };
};
struct D3Q27 {
	const unsigned int numDirections = 27;
	// lattice weights and lenght
	const double zeroWeight = 8.0 / 27.0;
	const double adjacentWeight = 2.0 / 27.0;
	const double adjacentLenght = 1.0;
	const double diagonalWeight = 1.0 / 54.0;
	const double diagonalLenght = sqrt(2.0);
	const double cornerWeight = 1.0 / 216.0;
	const double cornerLenght = sqrt(3.0);
	const int ciX[27] = { 0,1,-1,0,0,0,0,1,-1,1,-1,0,0,1,-1,1,-1,0,0,1,-1,1,-1,1,-1,-1,1 };
	const int ciY[27] = { 0,0,0,1,-1,0,0,1,-1,0,0,1,-1,-1,1,0,0,1,-1,1,-1,1,-1,-1,1,1,-1 };
	const int ciZ[27] = { 0,0,0,0,0,1,-1,0,0,1,-1,1,-1,0,0,-1,1,-1,1,1,-1,-1,1,1,-1,1,-1 };
	const double wi[27] = { zeroWeight,adjacentWeight,adjacentWeight,adjacentWeight,adjacentWeight,adjacentWeight,adjacentWeight,diagonalWeight,diagonalWeight,diagonalWeight,diagonalWeight,diagonalWeight,diagonalWeight,diagonalWeight,diagonalWeight,diagonalWeight,diagonalWeight,diagonalWeight,diagonalWeight,cornerWeight,cornerWeight,cornerWeight,cornerWeight,cornerWeight,cornerWeight,cornerWeight,cornerWeight };
	const double length[27] = { 0.0, adjacentLenght,adjacentLenght,adjacentLenght,adjacentLenght,adjacentLenght,adjacentLenght,diagonalLenght,diagonalLenght,diagonalLenght,diagonalLenght,diagonalLenght,diagonalLenght,diagonalLenght,diagonalLenght,diagonalLenght,diagonalLenght,diagonalLenght,diagonalLenght,cornerLenght,cornerLenght,cornerLenght,cornerLenght,cornerLenght,cornerLenght,cornerLenght,cornerLenght };
};