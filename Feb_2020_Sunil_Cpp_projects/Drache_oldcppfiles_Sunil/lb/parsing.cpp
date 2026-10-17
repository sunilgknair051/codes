// mcPolymer - Monte Carlo simulation program of radical polymerization
// Copyright (C) 2008 - 2018 Marco Drache, Benjamin Hosemann, Georg Drache, Katrin Brandl
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
#include <string>
#include <vector>

#include "parsing.h"

namespace parsing {
	vector<string> splitline(const string& _line) {
		vector<string> splitstring;
		const string SEP_CHAR = ",()[]=+";

		int i;
		size_t pos;
		int start;
		string clearstring, substr, linecopy;
		bool sepmode;
		linecopy = _line;

		// separate +,(,),[,],=
		for (i = 0; i < linecopy.length(); i++) {
			if (SEP_CHAR.find(linecopy.at(i) != string::npos)) {

			}
			if (((linecopy.at(i) == ',') || linecopy.at(i) == '(') || (linecopy.at(i) == ')') || (linecopy.at(i) == '[') || (linecopy.at(i) == ']') || (linecopy.at(i) == '=')) {
				if (i < (linecopy.length() - 1)) {
					linecopy.insert(i + 1, " ");
				}
				else {
					linecopy.append(" ");
				}

				linecopy.insert(i, " ");
				++i;
				++i;
			}
		}
		// replace separators
		pos = linecopy.find_first_of(" \t:;,");
		while (pos != std::string::npos)
		{
			linecopy[pos] = ' ';
			pos = linecopy.find_first_of(" \t:;,", pos + 1);
		}
		start = 0;
		// ignore leading separators 
		if (linecopy.at(0) == ' ') {
			for (i = 0; i < linecopy.length(); i++) {
				if (linecopy.at(i) == ' ') {
					++start;
				}
				else {
					break;
				}
			}
		}
		// delete duplicate separators
		sepmode = false;
		substr = "";
		for (i = start; i < linecopy.length(); i++) {
			if (linecopy.at(i) == ' ') {
				if (!sepmode) {
					sepmode = true;
					//close splitstring
					splitstring.push_back(substr);
					substr = "";
				}
			}
			else {
				substr.append(&linecopy.at(i), 1);
				sepmode = false;
			}
		}
		if (!sepmode) {
			if (substr != "") {
				splitstring.push_back(substr);
			}
		}

		return splitstring;
	}
}
