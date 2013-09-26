/**
 * @file:       CPU_Statistics.cpp
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
 * @brief 	Implementaion of  the GA statistics over islands
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


#include <sstream>
#include <malloc.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>

#include "CPU_Statistics.h"
#include "CPU_Statistics.h"



//----------------------------------------------------------------------------//
//                              Definitions                                   //
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
//                       TGPU_Statistics Implementation                       //
//                              public methods                                //
//----------------------------------------------------------------------------//


/**
 * Constructor of the class
 * 
 */
TCPU_Statistics::TCPU_Statistics()
{
        
    StatReceiveBuffer       = NULL;        
    ReceiveIndividualBuffer = NULL;
        
    BestIslandIdx = 0;
    
    AllocateMemory();
            
}// end of TGPU_Population
//------------------------------------------------------------------------------


/**
 * Destructor of the class
 * 
 */
TCPU_Statistics::~TCPU_Statistics()
{
    
    FreeMemory();
           
}// end of ~TGPU_Population
//------------------------------------------------------------------------------


/**
 * Calculate Statistics
 * 
 * @param [in]  Population - Pointer to population to calculate statistics over 
 * 
 */    
void TCPU_Statistics::Calculate(const TCPU_Population * Population)
{
    
    
    // Calculate local statistics
    CalculateLocalStatistics(Population);
   
    const TParameters * Params = TParameters::GetInstance();
        
    
    // Collect statistical data 
    MPI_Gather(&StatDataBuffer   ,sizeof(StatDataBuffer),MPI_BYTE,
                StatReceiveBuffer,sizeof(StatDataBuffer),MPI_BYTE, 0, MPI_COMM_WORLD);
    
    
    // Collect Individuals (to have the best one)
    MPI_Gather(&(Population->Population[LocalBestSolutionIdx*Params->ChromosomeSize()]), Params->ChromosomeSize(),MPI_UNSIGNED,
                ReceiveIndividualBuffer                                                , Params->ChromosomeSize(),MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    
    
    // Only master calculates the global statistics
    if (Params->IslandIdx() == 0) CalculateGlobalStatistics();
            
    
}// end of Calculate
//------------------------------------------------------------------------------


/**
 * Print best individual as a string
 *
 * @param [in] GlobalKnapsackData
 * @return Best individual in from of a sting 
 */   
string TCPU_Statistics::GetBestIndividualStr(const TGlobalKnapsackData * GlobalKnapsackData)
{

    stringstream  S; 
    
    TParameters * Params = TParameters::GetInstance();
    
    TGene * BestIndividualStart = ReceiveIndividualBuffer + BestIslandIdx * Params->ChromosomeSize();
    
    // Convert by eight bits
    for (int BlockID=0; BlockID < Params->ChromosomeSize() -1; BlockID++){
     
     for (int BitID = 0; BitID < Params->IntBlockSize() -1; BitID++ ) {         
         char c = ((BestIndividualStart[BlockID] & (1 << BitID)) == 0) ? '0' : '1';         
         S << c;
         if (BitID % 8 ==7) S << " ";
     }    
         
     S << "\n";
     
  }
 
 
     // Convert the remainder
    for (int BitID = 0; BitID < Params->IntBlockSize() - (GlobalKnapsackData->NumberOfItems - GlobalKnapsackData->OriginalNumberOfItems); BitID++) {
         char c =  ((BestIndividualStart[Params->ChromosomeSize() -1] & (1 << BitID)) == 0) ? '0' : '1';
         S << c;
         if (BitID % 8 ==7) S << " ";
    }
          
 
 return S.str();   
} // end of GetBestIndividualStr
//------------------------------------------------------------------------------

//----------------------------------------------------------------------------//
//                       TLocalStatistics Implementation                       //
//                              protected methods                             //
//----------------------------------------------------------------------------//

/**
 * Allocate memory
 */
void TCPU_Statistics::AllocateMemory()
{
    
    
    TParameters*  Params = TParameters::GetInstance();
    
    // Allocate MPI buffers for stat data and best individual
    if (Params->IslandIdx() == 0) {
      StatReceiveBuffer          = (TStatDataToExchange *) memalign(16, sizeof(TStatDataToExchange) * Params->IslandCount());          
      ReceiveIndividualBuffer    = (TGene *)               memalign(16, sizeof(TGene)  * Params->ChromosomeSize()* Params->IslandCount());
    }
        
    
}// end of AllocateMemory
//------------------------------------------------------------------------------

/**
 * Free GPU memory
 */
void TCPU_Statistics::FreeMemory()
{
       
    if (StatReceiveBuffer){
        free(StatReceiveBuffer);       
        StatReceiveBuffer = NULL;
    }
    if (ReceiveIndividualBuffer){
        free(ReceiveIndividualBuffer);
        ReceiveIndividualBuffer = NULL;
    }
    
}// end of FreeMemory
//------------------------------------------------------------------------------
    
    
/**
 * Calculate local statistics
 * 
 * @param [in] Population - Calculate local statistics over population
 * 
 */
void TCPU_Statistics::CalculateLocalStatistics(const TCPU_Population * Population)
{
   
    // initialize values
    StatDataBuffer.MinFitness = TFitness(UINT_MAX);
    StatDataBuffer.MaxFitness = TFitness(0);
    
    StatDataBuffer.SumFitness = 0.0f;
    StatDataBuffer.Sum2Fitness= 0.0f;   
    
    LocalBestSolutionIdx  = 0;         
    
    // Calculate statistics
    for (int i = 0; i < Population->PopulationSize; i++){
        if (StatDataBuffer.MaxFitness < Population->Fitness[i]){
            StatDataBuffer.MaxFitness = Population->Fitness[i];
            LocalBestSolutionIdx = i;
        }
        if (StatDataBuffer.MinFitness > Population->Fitness[i]){
            StatDataBuffer.MinFitness = Population->Fitness[i];            
        }
        StatDataBuffer.SumFitness  +=  Population->Fitness[i];
        StatDataBuffer.Sum2Fitness += (Population->Fitness[i]) * (Population->Fitness[i]);                
    }
    
    
}// end of CalculateLocalStatistics
//------------------------------------------------------------------------------
    
/**
 * Calculate global statistics
 */
void TCPU_Statistics::CalculateGlobalStatistics()
{
     
 BestIslandIdx = 0;
 
 
 // initialize values
 StatDataBuffer.MaxFitness  = StatReceiveBuffer[0].MaxFitness;
 StatDataBuffer.MinFitness  = StatReceiveBuffer[0].MinFitness;
 StatDataBuffer.SumFitness  = StatReceiveBuffer[0].SumFitness;
 StatDataBuffer.Sum2Fitness = StatReceiveBuffer[0].Sum2Fitness;
 
 TParameters * Params = TParameters::GetInstance();
 
  // Calculate numeric stats over receiving buffer
  for (int i = 1 ;  i < Params->IslandCount(); i++){      
      

      if (StatDataBuffer.MaxFitness < StatReceiveBuffer[i].MaxFitness) {
          StatDataBuffer.MaxFitness = StatReceiveBuffer[i].MaxFitness;
          BestIslandIdx = i;
      }
      
      if (StatDataBuffer.MinFitness > StatReceiveBuffer[i].MinFitness) {
          StatDataBuffer.MinFitness = StatReceiveBuffer[i].MinFitness;
      }
      
      StatDataBuffer.SumFitness  +=  StatReceiveBuffer[i].SumFitness;
      StatDataBuffer.Sum2Fitness +=  StatReceiveBuffer[i].Sum2Fitness;

  } // for
    

 
 // Calculate derived statistics
 DerivedStats.AvgFitness = StatDataBuffer.SumFitness / (Params->PopulationSize() * Params->IslandCount());
 DerivedStats.Divergence = sqrtf(fabsf( (StatDataBuffer.Sum2Fitness / (Params->PopulationSize() * Params->IslandCount()) -
                                         DerivedStats.AvgFitness * DerivedStats.AvgFitness))
                           );
 
 
 
}// end of CalculateGlobalStatistics
//------------------------------------------------------------------------------


/**
 * Set ONLY best individual Idx and its fitness to statistics
 * @param [in] Population to find the best solution and set max index
 * 
 */
void  TCPU_Statistics::SetBestLocalIndividualAndMaxFintess(const TCPU_Population * Population)
{
    
    StatDataBuffer.MaxFitness = TFitness(0);
        
    LocalBestSolutionIdx  = 0;         
    
    // Find the best insladn
    for (int i = 0; i < Population->PopulationSize; i++){
        if (StatDataBuffer.MaxFitness < Population->Fitness[i]){
            StatDataBuffer.MaxFitness = Population->Fitness[i];
            LocalBestSolutionIdx = i;
        }        
    }
    
 
}// end of SetBestLocalIndividualAndMaxFintess
//------------------------------------------------------------------------------
