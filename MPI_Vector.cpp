//*************************************************************************************************
// MPI_Vector.cpp
//
// The MPI example demonstrates the setup of a new MPI data type "vector". This new data type
// allows "vector"-messages that may be columns of a matrix.
//*************************************************************************************************


// Includes
#include <mpi.h>
#include <iostream>

using namespace std;

// Global constants
const int SIZE = 4;


//*************************************************************************************************
int main( int argc, char* argv[] )
{
   // Definition and initialization of the variables
   int size;  //The total number of processes
   int rank;  //The rank/number of this process

   // Setting the sending process to 0
   const int source = 0;

   // Setting the message tag (offers the possibility to give
   //  each message an explicit number)
   const int tag = 1;

   // Definition and initialization of the data field a
   float a[SIZE][SIZE] = {  1.0,  2.0,  3.0,  4.0,
                            5.0,  6.0,  7.0,  8.0,
                            9.0, 10.0, 11.0, 12.0,
                           13.0, 14.0, 15.0, 16.0 };

   // Definition of the receive buffer b
   float b[SIZE];

   // MPI_Status structure: the structure contains information
   //  about the source and message tag (MPI_SOURCE and MPI_TAG)
   MPI_Status stat;

   // The new MPI datatype
   MPI_Datatype columntype;

   // Initialization of MPI
   // TODO

   // Initialization of each CPU: determine the total number
   //  of processes and the rank of each process
   // TODO


   // Creating a new derived data type
   MPI_Type_vector( SIZE, 1, SIZE, MPI_FLOAT, &columntype );
   MPI_Type_commit( &columntype );

   // Let the CPU 0 send a column of 'a' to each of the 4 CPUs,
   //  which receive this column in the receive array 'b'
   if( size == SIZE ) {

      // TODO






      // Output
      cout << " rank= " << rank << " , b= " << b[0] << ", " << b[1] << ", " << b[2] << ", " << b[3] << "\n";
   }
   else {
      cout << " Must specify " << SIZE << " processors. Terminating.\n";
   }

   // MPI finalizations
   // TODO


   return 0;
}
//*************************************************************************************************
