/* 
 * File:        CPU_Evolution.h
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
 * Comments:    Header file of the GA evolution
 *              This class controls the evolution process on multicore CPU
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
 

#ifndef CPU_EVOLUTION_H
#define	CPU_EVOLUTION_H


#include "Parameters.h"
#include "CPU_Population.h"
#include "CPU_Statistics.h"
#include "GlobalKnapsackData.h"


/*
 * struct of seed
 */
struct r123_seed{
    unsigned long seed1;
    unsigned long seed2;
};// end of r123_seed
//------------------------------------------------------------------------------



/*
 * CPU evolution process
 * 
 */
class TCPU_Evolution{
public:
    
    TCPU_Evolution(int argc, char **argv);
    virtual ~TCPU_Evolution();
    
    // Run evolution
    void Run();
    
    // Is this a master process?
    bool IsMaster() {return Params->IslandIdx() == 0;};
    
protected:    
    TParameters * Params;                 // Parameters of evolution
    int ActGeneration;                    // Actual generation            
    
    
    TCPU_Population* MasterPopulation;    // Master GA population  
    TCPU_Population* OffspringPopulation; // Population of offsprings

    
    TCPU_Population* EmigrantsToSend;     // Population of emigrants to Send
    TCPU_Population* EmigrantsToReceive;  // Population of immigrants to recieve
    
    
    TCPU_Statistics *   CPUStatistics;    // Statistics over GA process    
    TGlobalKnapsackData GlobalData;       // Global data of knapsack
    
    // Get Random Seed for Random123
    r123_seed  GetSeed();
    
    
    // Initialize
    void Initialize();        
    // Run evolution cycle
    void RunEvolutionCycle();
    
    // Generate first population
    void GenerateFirstPopulation(); 
    
    // Genetic manipulation
    void GeneticManipulation();
    
    // Replacement phase
    void Replacement();
    
    // Migrate phase
    void Migrate();
        
    // Pack Emigrants to a buffer for dispatch
    void PackEmigrantsToSend();
    
    // Unpack immigrants form buffer
    void UnpackReceivedEmigrants();
    
    // Selection of two chromosomes
    inline unsigned int Selection(TCPU_Population * ParentsData, unsigned int Random1, unsigned int Random2);
    // Uniform Crossover    
    inline void CrossoverUniformFlip(TGene& GeneOffspring1, TGene& GeneOffspring2,
                                     TGene& GeneParent1   , TGene& GeneParent2,
                                     unsigned int RandomValue);
    // Bit Flip mutaion Mutation 
    inline void MutationBitFlip(TGene& GeneOffspring1, TGene& GeneOffspring2,
                                unsigned int RandomValue1,unsigned int RandomValue2, int BitID);
        
    
    TCPU_Evolution(const TCPU_Evolution& orig);
};



#endif	/* CPU_EVOLUTION_H */

