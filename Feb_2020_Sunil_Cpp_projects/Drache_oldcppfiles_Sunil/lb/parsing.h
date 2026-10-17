// mcPolymer - Monte Carlo simulation program of radical polymerization
// Copyright (C) 2008 - 2017 Marco Drache, Benjamin Hosemann, Georg Drache, Katrin Kipker
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
/*
/*
*
*  Created on: 3.9.2017
*/
#pragma once
#include <string>
#include <vector>

namespace parsing {
	using namespace std;
	using FileDaten = vector<string>;
	vector<string> splitline(const string& _line);
}
