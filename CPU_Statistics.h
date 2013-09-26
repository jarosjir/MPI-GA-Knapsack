/**
 * @file:       CPU_Statistics.h
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
 * @brief 	Header file of the GA statistics over islands
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
		26 September 2013, 10:00 (revised)
 */

#ifndef CPU_STATISTICS_H
#define	CPU_STATISTICS_H


#include "Parameters.h"
#include "CPU_Population.h"
#include "GlobalKnapsackData.h"

/**
 * @struct TStatDataToExchange
 * @brief  Statistic structure for exchange between nodes
 */
struct TStatDataToExchange{
  /// Minimum fitness value in population
  TFitness MinFitness;
  /// Minimum fitness value in population        
  TFitness MaxFitness;
  /// Sum of all fitness values over population
  TFitness SumFitness;   
  /// Sum of squares of all fitness values over poopulation
  TFitness Sum2Fitness;               
};
//------------------------------------------------------------------------------

/**
 * @class TDerivedStats
 * @brief Derived statistics structure
 */
struct TDerivedStats{
    /// Average Fintess
    float AvgFitness;
    /// Divergence
    float Divergence;
};

/**
 * @class TCPU_Statistics
 * @brief Class that maintains local and global statistic of the island GA
 * 
 */
class TCPU_Statistics {
public:
    
    /// Get maximum fitness value over all nodes, only master can read!
    TFitness GetMaxFitness()    const {return StatDataBuffer.MaxFitness;};
    /// Get minimum fitness value over all nodes, only master can read!
    TFitness GetMinFitness()    const {return StatDataBuffer.MinFitness;};
    /// Get average fitness value over all nodes, only master can read!
    float    GetAvgFitness()    const {return DerivedStats.AvgFitness;};
    /// Get fitness divergence value over all nodes, only master can read!
    float    GetDivergence()    const {return DerivedStats.Divergence;};
    /// Get the ID of the best island over all nodes, only master can read!
    int      GetBestIslandIdx() const {return BestIslandIdx;};
    /// Get the best local solution index 
    int      GetBestLocalSolutionIdx() const {return LocalBestSolutionIdx;};
    
    /// Get best individual chromosome as a string, only master can call this function!
    string   GetBestIndividualStr(const TGlobalKnapsackData * GlobalKnapsackData);

    /// Calculate statistics - all islands must call this routine
    void     Calculate(const TCPU_Population * Population);       
    /// Set best individual and max fitness for the local node
    void     SetBestLocalIndividualAndMaxFintess(const TCPU_Population * Population);
    
    /// Constructor
             TCPU_Statistics();
    /// Destructor
    virtual ~TCPU_Statistics();

protected:
    /// Statistics data to exchange
    TStatDataToExchange   StatDataBuffer;       
    /// Buffer to receive statistics
    TStatDataToExchange * StatReceiveBuffer;    
    
    /// Derived statistics
    TDerivedStats DerivedStats;                 
    /// Index to the local best solution
    int           LocalBestSolutionIdx;         
    /// Best island (island with the best solution)
    int           BestIslandIdx;                
       
    /// Buffer for best solution coming over network
    TGene    *    ReceiveIndividualBuffer;      
        
    /// Allocate Memory    
    void AllocateMemory();                      
    /// Free Memory
    void FreeMemory();                          
             
    /// Calculate local stats
    void CalculateLocalStatistics(const TCPU_Population * Population); 
    /// Calculate global stats 
    void CalculateGlobalStatistics();           
    
private:
    /// copy constructor not allowed
    TCPU_Statistics(const TCPU_Population& orig);

};// end of TGPU_Statistics
//------------------------------------------------------------------------------






#endif	/* CPU_STATISTICS_H */

