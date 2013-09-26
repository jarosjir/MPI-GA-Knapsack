/**
 * @file:       CPU_Population.h
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
 * 
 * @brief 	Header file of the GA population
 *              
 *
 * 
 * @section	License
 *		This source code is distribute under OpenSource GNU GPL license
 *                
 *              If using this code, please consider citation of related papers
 *              at http://www.fit.vutbr.cz/~jarosjir/pubs.php        
 *      
 *
 * 
 * @version	1.0
 * @date	06 June      2012, 00:00 (created)
		26 September 2013, 10:50 (revised)
 */



#ifndef CPU_POPULATION_H
#define	CPU_POPULATION_H


#include <string>

#include "GlobalKnapsackData.h"

using namespace std;

/**
 * @typedef TGene 
 * @brief Datatype for one gene 
 */
typedef unsigned int TGene;

/**
 * @typedef TFitness
 * @brief Datatype fitness value
 */
typedef float        TFitness;



/**
 * @class TCPU_Population
 * @brief CPU version of GA Population
 * 
 */
class TCPU_Population{
public:
    /// Number of chromosomes
    int PopulationSize;      
    /// Size of chromosome in INTs  
    int ChromosomeSize;      

    /// 1D array of genes (chromosome-based encoding)
    TGene    * Population;    
    /// 1D array of fitness values
    TFitness * Fitness;       

    /// Constructor  
    TCPU_Population(const int PopulationSize, const int ChromosomeSize);
        
    /// Get string representation of chromosome
    string GetStringOfChromosome(const int Idx);
    
    /// Calculate fitness values of all chromosomes
    void CalculateFitness(const TGlobalKnapsackData * GlobalKnapsackData);
    
    /// Calculate fitenss values of all chromosome, SSE4.1 version (works on Nehalem and newer cores)
    void CalculateFitness_SSE(const TGlobalKnapsackData * GlobalKnapsackData);

    
    /// Destructor
    virtual ~TCPU_Population();

protected:
    
    /// Memory allocation
    void AllocateMemory();
    
    /// Memory dealocation
    void FreeMemory();
   
    
private:
    /// Default constructor not allowed
    TCPU_Population();            
    
    /// Copy constructor not allowed
    TCPU_Population(const TCPU_Population& orig);
        
};// end of TCPU_Population
//------------------------------------------------------------------------------



#endif	/* CPU_POPULATION_H */

