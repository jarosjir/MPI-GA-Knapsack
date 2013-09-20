/* 
 * File:        CPU_Statistics.h
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
 * Comments:    Header file of the GA statistics
 *              This class maintains and collects GA statistics
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

#ifndef CPU_STATISTICS_H
#define	CPU_STATISTICS_H


#include "Parameters.h"
#include "CPU_Population.h"
#include "GlobalKnapsackData.h"

/*
 * Statistic structure for exchange between nodes
 */
struct TStatDataToExchange{
  TFitness MinFitness;   // Minimum fitness value in population
  TFitness MaxFitness;   // Minimum fitness value in population     
  
  TFitness SumFitness;   // Sum of all fitness values over population
  TFitness Sum2Fitness;  // Sum of squares of all fitness values over poopulation             
};
//------------------------------------------------------------------------------

/*
 * Derived statistics structure
 */
struct TDerivedStats{
    float AvgFitness;
    float Divergence;
};

/*
 * CPU statistics class
 * 
 */
class TCPU_Statistics {
public:
    
    //-- Only master process can read these values
    TFitness GetMaxFitness()    const {return StatDataBuffer.MaxFitness;};
    TFitness GetMinFitness()    const {return StatDataBuffer.MinFitness;};
    float    GetAvgFitness()    const {return DerivedStats.AvgFitness;};
    float    GetDivergence()    const {return DerivedStats.Divergence;};
    int      GetBestIslandIdx() const {return BestIslandIdx;};
    int      GetBestLocalSolutionIdx() const {return LocalBestSolutionIdx;};
    
    // Only master can call this function
    string   GetBestIndividualStr(TGlobalKnapsackData * GlobalKnapsackData);

    // All threads MUST call this function
    void     Calculate(TCPU_Population * Population);       
    
    void     SetBestLocalIndividualAndMaxFintess(TCPU_Population * Population);
    
             TCPU_Statistics();
    virtual ~TCPU_Statistics();

protected:
    
    TStatDataToExchange   StatDataBuffer;       // statistics data to exchange
    TStatDataToExchange * StatReceiveBuffer;    // buffer to receive statistics
    
    TDerivedStats DerivedStats;                 // derived statistics
    int           LocalBestSolutionIdx;         // Index to the local best solution
    int           BestIslandIdx;                // Best island (island with the best solution)
       
    TGene    *    ReceiveIndividualBuffer;      // Buffer for best solution coming over network
        
    
    void AllocateMemory();                      // Allocate Memory
    void FreeMemory();                          // Free Memory
             
    void CalculateLocalStatistics(TCPU_Population * Population); // Calculate local stats
    void CalculateGlobalStatistics();           // Calculate global stats 
    
private:
    TCPU_Statistics(const TCPU_Population& orig);

};// end of TGPU_Statistics
//------------------------------------------------------------------------------






#endif	/* CPU_STATISTICS_H */

