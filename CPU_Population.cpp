/**
 * @file:       CPU_Population.cpp
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
 * @brief 	Implementation file of the GA population
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

#include <sstream>
#include <smmintrin.h>

#include "malloc.h"
#include "string.h"
#include "CPU_Population.h"
#include "Parameters.h"




//----------------------------------------------------------------------------//
//                              Definitions                                   //
//----------------------------------------------------------------------------//

    
//----------------------------------------------------------------------------//
//                       TCPU_Population Implementation                       //
//                              public methods                                //
//----------------------------------------------------------------------------//

/**
 * Constructor of the class
 * 
 * @param [in] PopulationSize
 * @param [in] ChromosomeSize
 * 
 */
TCPU_Population::TCPU_Population(const int PopulationSize, const int ChromosomeSize)
{    
    this->ChromosomeSize = ChromosomeSize;
    this->PopulationSize = PopulationSize;
    
    AllocateMemory();    
    
}// end of TCPU_Population
//------------------------------------------------------------------------------
    

/**
 * Destructor of the class
 */
TCPU_Population::~TCPU_Population()
{    
    FreeMemory();
       
}// end of TCPU_Population
//------------------------------------------------------------------------------
    


/**
 * Print chromosome to string
 * 
 * @param [in] Idx Idx of chromosome in population
 * @return String representation of chromosome
 */
string TCPU_Population::GetStringOfChromosome(const int Idx)
{
       
 stringstream  S;    

  // simple print of chromosome 
 for (int BlockID=0; BlockID<ChromosomeSize; BlockID++){
     
     for (int BitID = 0; BitID < 32; BitID++ ) {
         char c = ((Population[Idx*ChromosomeSize + BlockID] & (1 << BitID)) == 0) ? '0' : '1';
         S << c;
     }    
         
     S << "\n";
     
  }
 
 return S.str();   
    
}// end of GetStringOfChromozome
//------------------------------------------------------------------------------


/**
 * Calculate fitness using SSE 4.1 instructions
 * 
 * @param [in] GlobalKnapsackData
 * 
 */
void TCPU_Population::CalculateFitness_SSE(const TGlobalKnapsackData * GlobalKnapsackData)
{
 
  TParameters * Params = TParameters::GetInstance();  
    
  
  for ( int ChromosomeIdx = 0; ChromosomeIdx < PopulationSize; ChromosomeIdx++){
      int StartIdx = ChromosomeIdx * ChromosomeSize;
  
        int TempArray[4] __attribute__((aligned(16))) ;
      
      __m128i Part_Weight_SSE = _mm_setzero_si128 ();
      __m128i Part_Price_SSE  = _mm_setzero_si128 ();      
      
      __m128i Mask = _mm_set1_epi32(1);

      
      // Calculate prices and weights using SSE 
      for (int BlockID = 0; BlockID < ChromosomeSize; BlockID++) {          
          
          unsigned int Gene = Population[StartIdx + BlockID];
          __m128i GeneValue_SSE = _mm_set_epi32 (Gene >>3, Gene >> 2, Gene >> 1, Gene >>0);
          
          for (int i = 0; i < Params->IntBlockSize(); i+=4){
             size_t index = BlockID * Params->IntBlockSize() + i ;
   
             __m128i Weights_SSE    = _mm_load_si128((__m128i *)&(GlobalKnapsackData->ItemWeight[index]));             
             __m128i Prices_SSE     = _mm_load_si128((__m128i *)&(GlobalKnapsackData->ItemPrice[index]));             
                                                                

             __m128i ActualGene_SSE = _mm_srli_epi32(GeneValue_SSE, i);
             __m128i ItemValue_SSE  = _mm_and_si128 (ActualGene_SSE,Mask);    
                                      
             
             Weights_SSE     = _mm_mullo_epi32(Weights_SSE, ItemValue_SSE);
             Prices_SSE      = _mm_mullo_epi32(Prices_SSE, ItemValue_SSE);
             
             Part_Weight_SSE = _mm_add_epi32(Weights_SSE, Part_Weight_SSE);
             Part_Price_SSE  = _mm_add_epi32(Prices_SSE,  Part_Price_SSE);             
          }
      }  
  
      _mm_store_si128((__m128i *) TempArray, Part_Weight_SSE);
      int TotalWeight = (TempArray[0] + TempArray[1]) + (TempArray[2] + TempArray[3]);
      
      _mm_store_si128((__m128i *) TempArray, Part_Price_SSE);
      int TotalPrice = (TempArray [0] + TempArray [1]) + (TempArray [2] + TempArray [3]);
      
  
      // penalization by best price/weight ratio    
      TFitness result = (float) TotalPrice;

      if (TotalWeight > GlobalKnapsackData->KnapsackCapacity){
          int Penalty = TotalWeight - GlobalKnapsackData->KnapsackCapacity;

          result = result  - GlobalKnapsackData->MaxPriceWightRatio * (Penalty);            
          if (result < 0 ) result = TFitness(0);                  
      }                           
      
      Fitness[ChromosomeIdx] = result;
      
  }// end master for
  
    
}// end of CalculateFitness
//-----------------------------------------------------------------------------


/**
 * Calculate fitness, normal version
 * 
 * @param [in] GlobalKnapsackData
 * 
 */
void TCPU_Population::CalculateFitness(const TGlobalKnapsackData * GlobalKnapsackData)
{
 
  TParameters * Params = TParameters::GetInstance();  
    
  static const unsigned int ONE = 1; 
  
  
  for (int ChromosomeIdx = 0; ChromosomeIdx < PopulationSize; ChromosomeIdx++){
       int StartIdx = ChromosomeIdx * ChromosomeSize;
  
      
      int weight[4]={0,0,0,0};
      int price [4]={0,0,0,0};
      int Value [4]={0,0,0,0};

      
      // Calculate price and weight, loop is unrolled 
      for (int BlockID = 0; BlockID < ChromosomeSize; BlockID++) {
          for (int i = 0; i < Params->IntBlockSize(); i+=4){
             size_t index = BlockID * Params->IntBlockSize() + i ;

             Value[0] = (Population[StartIdx + BlockID] >> (i)     ) & ONE;
             Value[1] = (Population[StartIdx + BlockID] >> (i+1)   ) & ONE;
             Value[2] = (Population[StartIdx + BlockID] >> (i+2)   ) & ONE;
             Value[3] = (Population[StartIdx + BlockID] >> (i+3)   ) & ONE;                          

             weight[0] += Value[0] * GlobalKnapsackData->ItemWeight[index];
             price [0] += Value[0] * GlobalKnapsackData->ItemPrice [index];

             weight[1] += Value[1] * GlobalKnapsackData->ItemWeight[index + 1];
             price [1] += Value[1] * GlobalKnapsackData->ItemPrice [index + 1];

             weight[2] += Value[2] * GlobalKnapsackData->ItemWeight[index + 2];
             price [2] += Value[2] * GlobalKnapsackData->ItemPrice [index + 2];

             weight[3] += Value[3] * GlobalKnapsackData->ItemWeight[index + 3];
             price [3] += Value[3] * GlobalKnapsackData->ItemPrice [index + 3];
          }
      }  
  
      int TotalPrice = (price [0] + price [1]) + (price [2] + price [3]);
      int TotalWidth = (weight[0] + weight[1]) + (weight[2] + weight[3]);
  
            
      //-- penalization --//   
      TFitness result = (float) TotalPrice;

      if (TotalWidth > GlobalKnapsackData->KnapsackCapacity){
          unsigned int Penalty = TotalWidth - GlobalKnapsackData->KnapsackCapacity;

          result = result  - GlobalKnapsackData->MaxPriceWightRatio * (Penalty);            
          if (result < 0 ) result = TFitness(0);                  
      }                           
      
      Fitness[ChromosomeIdx] = result;
      
  }// end master for
      
}// end of CalculateFitness
//-----------------------------------------------------------------------------






//----------------------------------------------------------------------------//
//                       TCPU_Population Implementation                       //
//                            protected methods                               //
//----------------------------------------------------------------------------//

/**
 * Allocate memory
 */
void TCPU_Population::AllocateMemory()
{
    
    Population = (TGene    *) memalign(16, sizeof(TGene)    * ChromosomeSize * PopulationSize );
    Fitness    = (TFitness *) memalign(16, sizeof(TFitness) * PopulationSize );
    
}// end of AllocateMemory
//------------------------------------------------------------------------------

/** 
 * Free memory
 */
void TCPU_Population::FreeMemory()
{

    free(Population);
    free(Fitness);
}// end of FreeMemory
//------------------------------------------------------------------------------
