/** 
 * @file        GlobalKnapsackData.cpp
 * @author      Jiri Jaros
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
 * @brief       Implementation file of the knapsack global data class. 
 *              
 * @section     This source code is distribute under OpenSource GNU GPL license
 *                
 *              If using this code, please consider citation of related papers
 *              at http://www.fit.vutbr.cz/~jarosjir/pubs.php        
 *      
 *
 * @version	1.0
 * @date	06 June      2012, 00:00 (created)
 *		26 September 2013, 09:58 (revised)
 */




#include <malloc.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <mpi.h>

#include "GlobalKnapsackData.h"
#include "Parameters.h"


//----------------------------------------------------------------------------//
//                              Definitions                                   //
//----------------------------------------------------------------------------//




//----------------------------------------------------------------------------//
//                              public methods                                //
//----------------------------------------------------------------------------//

/**
 * Constructor of the class
 */
TGlobalKnapsackData::TGlobalKnapsackData()
{
    
    NumberOfItems               = 0;
    OriginalNumberOfItems       = 0;
    KnapsackCapacity            = 0;
    MaxPriceWightRatio          = 0;
    
    ItemPrice                   = NULL;
    ItemWeight                  = NULL;
       
}// end of constructor
//------------------------------------------------------------------------------


/**
 * Destructor of the class
 */   
TGlobalKnapsackData::~TGlobalKnapsackData()
{
    
    FreeMemory();
    
}// end of TGlobalKnapsackData
//------------------------------------------------------------------------------
    


/**
 * Load data from file, filename given in Parameter class
 */
void TGlobalKnapsackData::LoadFromFile()
{
    
    
    // Get instance of Parameter class    
    TParameters * Params = TParameters::GetInstance();

    // Open file with benchmark data
    ifstream fr(Params->BenchmarkFileName().c_str()); 
    if (!fr.is_open()) {
         cerr << "Global Benchmark Data: File not found" << endl;       
         Params->PrintUsageAndExit();
    }
    
    
    // Read number of items
    NumberOfItems = 0;    
    fr >>NumberOfItems;
    
    OriginalNumberOfItems = NumberOfItems;
    
    // Calculate padding    
    int Overhead = NumberOfItems % Params->IntBlockSize();
    if (Overhead != 0) NumberOfItems = NumberOfItems + (Params->IntBlockSize() - Overhead);
    
    // Allocate memory for arrays    
    AllocateMemory(NumberOfItems);
        
    
    // Load weights 
    for (int i = 0; i < OriginalNumberOfItems; i++){
        fr >> ItemPrice[i];    
                                   
    }    
    for (int i = OriginalNumberOfItems; i < NumberOfItems; i++){
        ItemPrice[i] = TPriceType(0);                       
    }

    
    
    // Load weights         
    for (int i = 0; i < OriginalNumberOfItems; i++){
        fr >> ItemWeight[i];
    }

    for (int i = OriginalNumberOfItems; i < NumberOfItems; i++){
        ItemWeight[i] = TPriceType(0);        
    }
    
    
    // Get max Price/Weight ratio
    MaxPriceWightRatio = 0.0f;
    
    for (int i = 0; i < OriginalNumberOfItems; i++){
        if (ItemWeight[i] != 0) {
                float Ratio = ItemPrice[i] / ItemWeight[i];
                if (Ratio > MaxPriceWightRatio)  MaxPriceWightRatio = Ratio;
        }
        
    }
    
    //Read Knapsack capacity    
    fr >> KnapsackCapacity;
      
    // Update chromosome size in parameters    
    Params->SetChromosomeSize(NumberOfItems/Params->IntBlockSize());
    
    
    
}// end of LoadFromFile
//------------------------------------------------------------------------------
    
    
    
//----------------------------------------------------------------------------//
//                           protected methods                                //
//----------------------------------------------------------------------------//

/**
 * Allocate memory
 * 
 * @param       NumberOfItems - Number of Items in Knapsack with padding
 */
void TGlobalKnapsackData::AllocateMemory(const int NumberOfItems)
{    
    
    // Allocate class memory with 16B with alignment
    
    ItemPrice  = (TPriceType *)  memalign(16, sizeof(TPriceType)  * NumberOfItems );
    ItemWeight = (TWeightType *) memalign(16, sizeof(TWeightType) * NumberOfItems);
    
    
}// end of AllocateMemory
//------------------------------------------------------------------------------


/**
 * Free Memory
 */
void TGlobalKnapsackData::FreeMemory()
{

    // Free class memory
    free(ItemPrice);
    free(ItemWeight);    
}// end of AllocateMemory
//------------------------------------------------------------------------------


