//=================================================================================================
//
//  Advanced Programming Techniques (AdvPT)
//  Winter Term 2010
//  Assignment 1 - Task 2
//
//  Write a program to sum the numbers in a user-specified range, omitting the 'if' test that
//  sets the upper and lower bounds. Predict what happens if the input is the numbers 7 and 3,
//  in that order. Now run the program giving it the numbers 7 and 3, and see if the results
//  match your expectation. If not, restudy the behavior of 'for' and 'while' loops until you
//  understand what happened.
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
	int nIn1, nIn2;
	
	cin >> nIn1;
	cin >> nIn2;
	
	nStart = std::min(nIn1,nIn2);
	nEnd   = std::max(nIn1,nIn2);
	
	int nSum = 0;
	
	for (int nIt = nStart; nIt <= nEnd; ++nIt)
		nSum += nIt;
		
	cout << nSum << endl;

	return 0;	
}
