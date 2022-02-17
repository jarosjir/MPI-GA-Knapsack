/**
 * @file        GlobalKnapsackData.cpp
 * @author      Jiri Jaros
 *              Brno University of Technology
 *              Faculty of Information Technology
 *
 *              and
 *
 *              The Australian National University
 *              ANU College of Engineering & Computer Science
 *
 *              jarosjir@fit.vutbr.cz
 *              www.fit.vutbr.cz/~jarosjir
 *
 * @brief       Implementation file of the knapsack global data class.
 *              This class maintains the benchmark data.
 *
 * @date        06 June      2012, 00:00 (created)
 *              15 February  2022, 12:33 (revised)
 *
 * @copyright   Copyright (C) 2012 - 2022 Jiri Jaros.
 *
 * This source code is distribute under OpenSouce GNU GPL license.
 * If using this code, please consider citation of related papers
 * at http://www.fit.vutbr.cz/~jarosjir/pubs.php
 *
 */

#include <mpi.h>
#include <malloc.h>
#include <fstream>


#include "GlobalKnapsackData.h"
#include "Parameters.h"

//--------------------------------------------------------------------------------------------------------------------//
//--------------------------------------------------- Definitions ----------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

static const char * const ERROR_FILE_NOT_FOUND = "Global Benchmark Data: File not found\n";

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Public methods ---------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Constructor of the class
 */
GlobalKnapsackData::GlobalKnapsackData()
  : numberOfItems(0),
    originalNumberOfItems(0),
    knapsackCapacity(0),
    maxPriceWightRatio(0),
    itemPrice(nullptr),
    itemWeight(nullptr)
{

}// end of constructor
//----------------------------------------------------------------------------------------------------------------------

/**
 * Destructor of the class
 */
GlobalKnapsackData::~GlobalKnapsackData()
{
  freeMemory();
}// end of destructor.
//----------------------------------------------------------------------------------------------------------------------

/**
 * Load data from file, filename given in Parameter class.
 */
void GlobalKnapsackData::loadFromFile()
{
    using namespace std;

    // Get instance of Parameter class
    Parameters& params = Parameters::getInstance();

    // Open file with benchmark data
    ifstream fr(params.getBenchmarkFileName().c_str());
    if (!fr.is_open())
    {
      fprintf(stderr, ERROR_FILE_NOT_FOUND);
      params.printUsageAndExit();
    }

    // Read number of items
    numberOfItems = 0;
    fr >> numberOfItems;

    originalNumberOfItems = numberOfItems;

    // Calculate padding
    int Overhead = numberOfItems % params.getIntBlockSize();
    if (Overhead != 0)
    {
      numberOfItems = numberOfItems + (params.getIntBlockSize() - Overhead);
    }

    // Allocate memory for arrays
    allocateMemory(numberOfItems);


    // Load weights
    for (int i = 0; i < originalNumberOfItems; i++)
    {
      fr >> itemPrice[i];
    }
    for (int i = originalNumberOfItems; i < numberOfItems; i++)
    {
      itemPrice[i] = PriceType(0);
    }


    // Load weights
    for (int i = 0; i < originalNumberOfItems; i++)
    {
      fr >> itemWeight[i];
    }

    for (int i = originalNumberOfItems; i < numberOfItems; i++)
    {
      itemWeight[i] = PriceType(0);
    }


    // Get max Price/Weight ratio
    maxPriceWightRatio = 0.0f;

    for (int i = 0; i < originalNumberOfItems; i++)
    {
      if (itemWeight[i] != 0)
      {
        float Ratio = itemPrice[i] / itemWeight[i];
        if (Ratio > maxPriceWightRatio)
        {
          maxPriceWightRatio = Ratio;
        }
      }
    }

    //Read Knapsack capacity
    fr >> knapsackCapacity;

    // Update chromosome size in parameters
    params.setChromosomeSize(numberOfItems / params.getIntBlockSize());
}// end of loadFromFile
//----------------------------------------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------------------------------//
//----------------------------------------------- Protected methods --------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Allocate memory.
 *
 * @param NumberOfItems - Number of Items in Knapsack with padding.
 */
void GlobalKnapsackData::allocateMemory(const int NumberOfItems)
{
  // Allocate class memory with 64B with alignment

  itemPrice  = (PriceType *)  memalign(64, sizeof(PriceType)  * NumberOfItems);
  itemWeight = (WeightType *) memalign(64, sizeof(WeightType) * NumberOfItems);
}// end of allocateMemory
//----------------------------------------------------------------------------------------------------------------------

/**
 * Free Memory
 */
void GlobalKnapsackData::freeMemory()
{
  // Free class memory
  free(itemPrice);
  free(itemWeight);
}// end of freeMemory
//----------------------------------------------------------------------------------------------------------------------
