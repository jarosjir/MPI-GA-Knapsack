/* 
 * File:        main.cpp 
 * Author:      Jiri Jaros
 * Affiliation: Brno University of Technology
 *              Faculty of Information Technology
 *              
 *              and
 * 
 *              The Australian National University
 *              ANU College of Engineering & Computer Science
 *
 * Email:       jarosjir@fit.vutbr.cz
 * Web:         www.fit.vutbr.cz/~jarosjir
 * 
 * Comments:    Efficient MPI implementation of the Island-Based, 
 *              Genetic Algorithm solving the Knapsack problem.
 *
 * 
 * License:     This source code is distribute under OpenSource GNU GPL license
 *                
 *              If using this code, please consider citation of related papers
 *              at http://www.fit.vutbr.cz/~jarosjir/pubs.php        
 *      
 *
 * 
 * Created on 06 June 2012, 00:00 PM
 */


#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <mpi.h>

#include "CPU_Evolution.h"
#include "Parameters.h"


using namespace std;



/*
 * Main Function
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

    return 0;
}

