//
// JAN23
//
// Created on: 24.1.2020

#include <iostream>
#include <fstream>
#include <vector>
#include "core2d.h"
#include "jan23.h"

using namespace std;

int main(int argc, char* argv[]) {
	Core2d* coreModule2D = new Core2d("2Dconfig.asc");
	coreModule2D->run(20);
	coreModule2D->writeResults("res");
	return 0;
}