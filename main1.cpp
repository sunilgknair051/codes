//=================================================================================================
//
//  Advanced Programming Techniques (AdvPT)
//  Winter Term 2010
//  Assignment 1 - Task 1
//
//  Write a program that prompts the user for two numbers and writes each number in the range
//  specified by the two numbers to the standard output.
//
//=================================================================================================

#include <iostream>
#include <algorithm>

using std::cout;
using std::cin;
using std::endl;

int main(int argc, char *argv[]) {
	cout << "Please specify two numbers" << endl;
	
	int nStart, nEnd;
	
	cin >> nStart;
	cin >> nEnd;
	
	if (nStart > nEnd) std::swap(nStart,nEnd);
	
	for (int nIt = nStart; nIt <= nEnd; ++nIt)
		cout << nIt << " ";
		
	cout << endl;

	return 0;	
}
