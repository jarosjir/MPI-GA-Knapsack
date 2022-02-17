/**
 * @file        main.cpp
 * @author      Jiri Jaros
 *              Brno University of Technology
 *              Faculty of Information Technology
 *
 *              and
 *              The Australian National University
 *              ANU College of Engineering & Computer Science
 *
 *              jarosjir@fit.vutbr.cz
 *              www.fit.vutbr.cz/~jarosjir
 *
 * @brief       Efficient MPI implementation of the Island-Based,
 *            	Genetic Algorithm solving the Knapsack problem.
 *
 * @version    	1.1
 * @date	      06 June      2012, 00:00 (created)
 *              15 February  2022, 10:58 (revised)
 *
 * @mainpage 	MPI-GA-knapsack
 *
 * @section	Usage
 *		Efficient MPI implementation of the Island-Based, Genetic Algorithm solving the Knapsack problem.
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
 *
 * @copyright   Copyright (C) 2012 - 2022 Jiri Jaros.
 *
 * This source code is distribute under OpenSouce GNU GPL license.
 * If using this code, please consider citation of related papers
 * at http://www.fit.vutbr.cz/~jarosjir/pubs.php
 *
 */

#include <mpi.h>

#include "CPU_Evolution.h"
#include "Parameters.h"


/**
 * The main function.
 * @param [in] argc - Number of command line arguments.
 * @param [in] argv - Command line arguments.
 * @return error code
 */
int main(int argc, char **argv)
{
  // MPI initialization
  MPI_Init(&argc, &argv);

  // Create evolution class on CPU.
  Evolution evolution(argc,argv);

  double algorithmStartTime = MPI_Wtime();

  // Run evolution
  evolution.run();

  double algorithmStopTime = MPI_Wtime();

  // Master process prints execution time
  if (evolution.isMaster())
  {
    printf("Execution time: %0.3f s.\n",  algorithmStopTime - algorithmStartTime);
  }

  // MPI finalization
  MPI_Finalize();

  return EXIT_SUCCESS;;
}// end of main
//----------------------------------------------------------------------------------------------------------------------

