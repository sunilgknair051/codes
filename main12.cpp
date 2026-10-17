//=================================================================================================
//
//  Advanced Programming Techniques (AdvPT)
//  Winter Term 2010
//  Assignment 1 - Task 12
//
//  Write a program, that reads in several words at once and stores them in a 'std::vector' of
//  'std::string's. Process all words by converting them to upper case strings and afterwards
//  print them in reverse order on the screen. Use both the subscript operator and iterators to
//  traverse the vector!
//
//=================================================================================================


#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <ctype.h>

using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::string;
	
int main(int argc, char *argv[]) {
	
	cout << "Please specify a couple of words (terminate with EOF):" << endl;
		
	vector<string> vecStr;
	string strTmp;
	
	while(std::getline(cin, strTmp, ' '))
		vecStr.push_back(strTmp);
		
	for (unsigned int i = 0; i < vecStr.size(); ++i)
		transform(vecStr[i].begin(), vecStr[i].end(), vecStr[i].begin(), toupper);
	
	cout << endl;
	
	for (vector<string>::const_iterator c_it = vecStr.end() - 1; c_it >= vecStr.begin(); --c_it)
		cout << *c_it << " ";
	
	cout << endl;
	
	return 0;	
}

