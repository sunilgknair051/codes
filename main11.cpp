//=================================================================================================
//
//  Advanced Programming Techniques (AdvPT)
//  Winter Term 2010
//  Assignment 1 - Task 11
//
//  Write a program to create a 'std::vector<int>' with 10 elements of your choice. Print the
//  initial elements. Using an iterator, assign each element a value that is twice its current
//  value. Afterwards, print the content of the vector again (using 'const_iterator's).
//
//=================================================================================================

#include <cstdlib>
#include <ctime>
#include <vector>
#include <iostream>

using std::cout;
using std::cin;
using std::endl;
using std::vector;

int main(int argc, char *argv[]) {
	
	srand ( time(NULL) );
	
	vector<int> vec(10,0);
	
	for (vector<int>::iterator it = vec.begin(); it != vec.end(); ++it)
		*it = rand() % 100 + 1;
		
	for (vector<int>::const_iterator c_it = vec.begin(); c_it != vec.end(); ++c_it)
		cout << *c_it << " ";
	
	cout << endl;
	
	for (vector<int>::iterator it = vec.begin(); it != vec.end(); ++it)
		*it = *it * 2;
		
	for (vector<int>::const_iterator c_it = vec.begin(); c_it != vec.end(); ++c_it)
		cout << *c_it << " ";
	
	cout << endl;
	
	return 0;
}
