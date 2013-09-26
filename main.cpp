/**
 * @file	main.cpp 
 * @author	Jiri Jaros \n
 *   	 	Brno University of Technology \n
 *              Faculty of Information Technology \n
 *              
 *              and			\n
 * 
 *              The Australian National University	\n
 *              ANU College of Engineering & Computer Science	\n
 *
 * 		jarosjir@fit.vutbr.cz
 * 	        www.fit.vutbr.cz/~jarosjir
 * 
 * @brief 	Efficient MPI implementation of the Island-Based, \n
 *            	Genetic Algorithm solving the Knapsack problem.
 *
 * @version	1.0
 * @date	06 June      2012, 00:00 (created)
		20 September 2012, 13:00 (revised)

 * @mainpage 	MPI-GA-knapsack
 * 
 * @section     License 
 * 		This source code is distribute under OpenSource GNU GPL license
 *                
 *              If using this code, please consider citation of related papers
 *              at http://www.fit.vutbr.cz/~jarosjir/pubs.php        
 *      
 *
 *
 * @section	Usage
		Efficient MPI implementation of the Island-Based, 
 *            	Genetic Algorithm solving the Knapsack problem.
 * 		
 * \verbatim	
 ----------------------------------------------------------------
Parameters: 
  -p <integer value>	: Population size       (default 128)
  -g <integer value>	: Number of generations (default 100)

  -m <float value> 	: Mutation rate 	(default 0.01)
  -c <float value> 	: Crossover rate	(default 0.7)
  -o <float value> 	: Offspring rate	(default 0.5)

  -e <float value>	: Emigrant rate		(default 0.1)
  -n <integer value>	: Migration interval	(default 1)
  -s <integer value>	: Statistics interval	(default 1)

  -b 			: Print best individual
  -f <file_name>	: Benchmark_file_name	(default knapsack_data.txt)
\endverbatim
 */




#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <mpi.h>

#include "CPU_Evolution.h"
#include "Parameters.h"


using namespace std;



/** 
 *  The main function
 * @param argc
 * @param argv
 * @return error code
 */
int main(int argc, char **argv)
{
        
      // MPI initialization
    MPI_Init(&argc, &argv);

       // Create CPU evolution class
    TCPU_Evolution CPU_Evolution(argc,argv);    
    
    double AlgorithmStartTime;
    AlgorithmStartTime = MPI_Wtime();
    
      // Run evolution
    CPU_Evolution.Run();
    
    
    double AlgorithmStopTime = MPI_Wtime();
    
      // Master process prints execution time
    if (CPU_Evolution.IsMaster()) printf("Execution time: %0.3f s.\n",  AlgorithmStopTime - AlgorithmStartTime);    
                                  
      // MPI finalization
    MPI_Finalize();

    return EXIT_SUCCESS;;
}// end of main
//------------------------------------------------------------------------------

