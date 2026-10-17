//
// MAR 06
//
// Created on: 06.03.2020

#include <iostream>
#include <fstream>
#include <vector>
#include "core3d.h"
#include "Mar06_2020.h"
#include <chrono>
#include <ctime>


using namespace std;

int main(int argc, char* argv[]) {
	Core3d* coreModule3D = new Core3d();
	coreModule3D->run();
	//CoreModule3D->writeResults("res");
	return 0;
}