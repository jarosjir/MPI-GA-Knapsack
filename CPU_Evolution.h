/**
 * @file:       CPU_Evolution.h
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
 * @brief 	Header file  of the island based GA
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
		26 September 2013, 11:05 (revised)
 */
 

#ifndef CPU_EVOLUTION_H
#define	CPU_EVOLUTION_H


#include "Parameters.h"
#include "CPU_Population.h"
#include "CPU_Statistics.h"
#include "GlobalKnapsackData.h"


/**
 * @struct r123_seed  
 * @brief struct of seed
 */
struct r123_seed{
    /// first seed
    unsigned long seed1;
    /// second seed
    unsigned long seed2;
};// end of r123_seed
//------------------------------------------------------------------------------



/**
 * @class TCPU_Evolution
 * @brief CPU evolution process
 * 
 */
class TCPU_Evolution{
public:
   
    /// Constructor  
    TCPU_Evolution(int argc, char **argv);
    /// Destructor
    virtual ~TCPU_Evolution();
    
    /// Run evolution
    void Run();
    
    /// Is this a master process?
    bool IsMaster() const {return Params->IslandIdx() == 0;};
    
protected:    
    /// Parameters of evolution
    TParameters * Params;                 
    /// Actual generation            
    int ActGeneration;                    
    
    /// Master GA population  
    TCPU_Population* MasterPopulation;    
    /// Population of offsprings
    TCPU_Population* OffspringPopulation; 

    /// Population of emigrants to Send
    TCPU_Population* EmigrantsToSend;     
    /// Population of immigrants to recieve
    TCPU_Population* EmigrantsToReceive;  
    
    /// Statistics over GA process
    TCPU_Statistics *   CPUStatistics;    
    /// Global data of knapsack    
    TGlobalKnapsackData GlobalData;       
    
    /// Get Random Seed for Random123
    r123_seed  GetSeed();
    
    
    /// Initialize
    void Initialize();        
    /// Run evolution cycle
    void RunEvolutionCycle();
    
    /// Generate first population
    void GenerateFirstPopulation(); 
    
    /// Genetic manipulation
    void GeneticManipulation();
    
    /// Replacement phase
    void Replacement();
    
    /// Migrate phase
    void Migrate();
        
    /// Pack Emigrants to a buffer for dispatch
    void PackEmigrantsToSend();
    
    /// Unpack immigrants form buffer
    void UnpackReceivedEmigrants();
    
    /// Selection of two chromosomes
    inline unsigned int Selection(const TCPU_Population * ParentsData, const unsigned int Random1, const unsigned int Random2);

    /// Uniform Crossover    
    inline void CrossoverUniformFlip(TGene& GeneOffspring1, TGene& GeneOffspring2,
                                     const TGene& GeneParent1, const TGene& GeneParent2,
                                     const unsigned int RandomValue);
    /// Bit Flip mutaion Mutation 
    inline void MutationBitFlip(TGene& GeneOffspring1, TGene& GeneOffspring2,
                                const unsigned int RandomValue1, const unsigned int RandomValue2, const int BitID);
        
    /// Copy constructor
    TCPU_Evolution(const TCPU_Evolution& orig);
};



#endif	/* CPU_EVOLUTION_H */

