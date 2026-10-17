//
// FEB 26
//
// Created on: 26.2.2020

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "core3d.h"
#include "jan23.h"

using namespace std;

int main(int argc, char* argv[]) {
	Core3d* coreModule3D = new Core3d("3Dconfig.asc");
	coreModule3D->run(20);
	coreModule3D->writeResults("res");
	return 0;
}