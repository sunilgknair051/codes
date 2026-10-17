//=================================================================================================
//
//  Advanced Programming Techniques (AdvPT)
//  Winter Term 2010
//  Assignment 1 - Task 3
//
//  Write a program that prompts the user to enter a number and then uses a 'for' loop to
//  calculate the faculty of the given number and writes it to the standard output.
//
//=================================================================================================

#include <iostream>

using std::cout;
using std::cin;
using std::endl;

int main(int argc, char *argv[]) {
	cout << "Please specify a number" << endl;

	int nNumber;
	long nFaculty = 1;
	
	cin >> nNumber;
	
	for (int nIt = nNumber; nIt > 1; --nIt)
		nFaculty *= nIt;
		
	cout << nFaculty << endl;

	return 0;	
}


