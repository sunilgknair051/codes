//=================================================================================================
//
//  Advanced Programming Techniques (AdvPT)
//  Winter Term 2010
//  Assignment 1 - Task 10
//
//  Write a program to strip the punctuation from a 'std::string'. The input to the program
//  should be a string of characters including punctuation; the output should be a 'std::string'
//  in which the punctuation is removed.
//
//=================================================================================================

#include <iostream>
#include <string>
#include <algorithm>
#include <ctype.h>

using std::cout;
using std::cin;
using std::endl;
using std::getline;
using std::string;
	
int main(int argc, char *argv[]) {
	
	cout << "Please specify a string:" << endl;
	
	string str;
	
	getline(cin, str);
	
	str.erase(std::remove_if(str.begin(), str.end(), ispunct), str.end());
	
	cout << str << endl;
	
	return 0;	
}

