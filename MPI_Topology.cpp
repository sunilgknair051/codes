//*************************************************************************************************
// MPI_Topology.cpp
//
// This MPI example demonstrates the setup of a MPI topology. The topology will consist of
// four processes arranged in a 2x2 fashion:
//
//      | 0  1 |
//      | 2  3 |
//
// After the setup of the MPI topology, each process is supposed to send his own rank to the
// neighboring processes. The neighbors will receive the messages and compare the received
// ranks to the expected ranks.
//*************************************************************************************************


// Includes
#include <mpi.h>
#include <iostream>

using namespace std;

// Global constants
const int SIZE  = 4;
const int UP    = 0;
const int DOWN  = 1;
const int LEFT  = 2;
const int RIGHT = 3;


//*************************************************************************************************
int main( int argc, char* argv[] )
{
   // Definition of the variables
   int size;    //The total number of processes
   int rank;    //The rank/number of this process
   int source;  //The source of a MPI message
   int dest;    //The destination for a MPI message
   int outbuf;  //Send buffer

   // Setting the message tag (offers the possibility to give
   //  each message an explicit number)
   const int tag = 1;

   // Definition and initialization of the receive buffer 'inbuf', where
   //  the received messages are stored and which will in this case
   //  contain the ranks of the neighboring CPUs
   int inbuf[4] = { MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL };

   // Definition of the neighbors array 'nrbs'
   //  Will contain the ranks of the neighbors in the order UP, DOWN, LEFT and RIGHT
   int nbrs[4];

   // Definition and initialization of the dimension array 'dims'
   //  In this case, we will construct a 2D-grid with 2 processes in x-
   //  and 2 processes in y-direction
   int dims[2] = {2,2};

   // Definition and initialization of the behavior of the Cartesian topology
   int periods[2] = {0,0};
   const int reorder = 0;

   // Definition of the coordinates array 'coords'
   int coords[2];

   // Definition of MPI request structures for the non-blocking communication
   MPI_Request reqs[8];

   // MPI_Status structure: the structure contains information
   //  about the source and message tag (MPI_SOURCE and MPI_TAG)
   MPI_Status stats[8];

   // The new MPI communicator with Cartesian topology
   MPI_Comm cartcomm;

   // Initialization of MPI
   // TODO

   // Determining the total number of CPUs
   // TODO

   // Creating the Cartesian topology and sending the own rank
   //  to the neighboring CPUs
   if( size == SIZE )
   {
      // Creating the Cartesian topology:
      //  - Creating a 2D grid with 2 processes in x- and 2 processes in y-direction
      //  - Determining the coordinates of the processes
      //  - Determining the neighbors in UP, DOWN, LEFT and RIGHT direction

      // TODO




      // Setting the output buffer
      outbuf = rank;

      // Sending the own rank to all four neighboring CPUs
      //  and receiving the rank of all four neighbors

      // TODO





      // Syncronisation of the CPUs
      // TODO

      // Output
      cout << " rank= " << rank << " , coords= " << coords[0] << "," << coords[1]
           << " , neighbors(u,d,l,r)= " << nbrs[UP] << "," << nbrs[DOWN] << ","
           << nbrs[LEFT] << "," << nbrs[RIGHT] << "\n";
      cout << " rank= " << rank << " ,                   inbuf(u,d,l,r)= "
           << inbuf[UP] << "," << inbuf[DOWN] << "," << inbuf[LEFT] << "," << inbuf[RIGHT] << "\n";
   }
   else {
      cout << " Must specify " << SIZE << " processors. Terminating.\n";
   }

   // MPI finalizations
   // TODO

   return 0;
}
//*************************************************************************************************
