/* Program to compute Pi using Monte Carlo methods for accurcy of 0.001 or itterration of 200000 what ever is reached first */

#include "montecarlosim.h"

int main()
{
	montecarlo(1.2652);   // change paramter for different radius of circle. These value also works as intial seed for randon number engine genrator.
	system("pause");
	optmontecarlo(1.2652);  // change parameter for different radius cirlce. These value also works as intial seed for randon number engine genrator...only first quadrant considered.
	system("pause");
	return 0;
}

