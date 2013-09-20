/* 
 * File:        CPU_Population.h
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
 * Comments:    Header file of the GA population
 *              This class maintains and GA populations
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


#ifndef CPU_POPULATION_H
#define	CPU_POPULATION_H


#include <string>

#include "GlobalKnapsackData.h"

using namespace std;

typedef unsigned int TGene;
typedef float        TFitness;



/*
 * CPU version of GA Population
 * 
 */
class TCPU_Population{
public:
    int PopulationSize;        // Number of chromosomes
    int ChromosomeSize;        // Size of chromosome in INTs

    
    TGene    * Population;    // 1D array of genes (chromosome-based encoding)
    TFitness * Fitness;       // 1D array of fitness values

        
    TCPU_Population(const int PopulationSize, const int ChromosomeSize);
        
        // Get string representation of chromosome
    string GetStringOfChromosome(const int Idx);
    
    // Calculate fitness values of all chromosomes
    void CalculateFitness(TGlobalKnapsackData * GlobalKnapsackData);
    
    // Calculate fitenss values of all chromosome, SSE4.1 version (works on Nehalem and newer cores)
    void CalculateFitness_SSE(TGlobalKnapsackData * GlobalKnapsackData);

    
    
    virtual ~TCPU_Population();

protected:
    
        // Memory allocation
    void AllocateMemory();
    
    // Memory dealocation
    void FreeMemory();

    
    
private:
        // default constructor not allowed
    TCPU_Population();            
    
        // copy constructor not allowed
    TCPU_Population(const TCPU_Population& orig);
        
};// end of TCPU_Population
//------------------------------------------------------------------------------



#endif	/* CPU_POPULATION_H */

