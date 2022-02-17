/**
 * @file        Statistics.cpp
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
 * @brief       Implementaion of  the GA statistics over islands
 *
 * @date        06 June      2012, 00:00 (created)
 *              15 February  2022, 15:00 (revised)
 *
 * @copyright   Copyright (C) 2012 - 2022 Jiri Jaros.
 *
 * This source code is distribute under OpenSouce GNU GPL license.
 * If using this code, please consider citation of related papers
 * at http://www.fit.vutbr.cz/~jarosjir/pubs.php
 */


#include <mpi.h>
#include <malloc.h>
#include <limits.h>
#include <cmath>

#include "Statistics.h"

//--------------------------------------------------------------------------------------------------------------------//
//--------------------------------------------------- Definitions ----------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Public methods ---------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Constructor of the class
 */
Statistics::Statistics()
  : mStatReceiveBuffer(nullptr),
    mLocalBestSolutionIdx(0),
    mBestIslandIdx(0),
    mReceiveIndividualBuffer(nullptr)
{
  allocateMemory();
}// end of Statistics
//----------------------------------------------------------------------------------------------------------------------

/**
 * Destructor of the class
 *
 */
Statistics::~Statistics()
{
  freeMemory();
}// end of ~Statistics
//----------------------------------------------------------------------------------------------------------------------

/**
 * Calculate Statistics.
 *
 * @param [in] population - Pointer to population to calculate statistics over
 */
void Statistics::calculate(const Population* population)
{
  // Calculate local statistics
  calculateLocalStatistics(population);

  const Parameters& params = Parameters::getInstance();

  // Collect statistical data
  MPI_Gather(&mStatDataBuffer   , sizeof(mStatDataBuffer), MPI_BYTE,
              mStatReceiveBuffer, sizeof(mStatDataBuffer), MPI_BYTE, 0, MPI_COMM_WORLD);


  // Collect Individuals (to have the best one)
  MPI_Gather(&(population->population[mLocalBestSolutionIdx * params.getChromosomeSize()]),
             params.getChromosomeSize(),
             MPI_UNSIGNED,
             mReceiveIndividualBuffer,
             params.getChromosomeSize(),
             MPI_UNSIGNED,
             0,
             MPI_COMM_WORLD);

  // Only master calculates the global statistics
  if (params.getIslandIdx() == 0)
  {
    calculateGlobalStatistics();
  }
}// end of calculate
//----------------------------------------------------------------------------------------------------------------------

/**
 * Print best individual as a string.
 */
std::string Statistics::getBestIndividualStr(const GlobalKnapsackData* globalKnapsackData)
{
 // Lambda function to convert 1 int into a bit string
  auto convertIntToBitString= [] (Gene value, int nValidDigits) -> std::string
  {
    std::string str = "";

    for (int bit = 0; bit < nValidDigits; bit++)
    {
      str += ((value & (1 << bit)) == 0) ? "0" : "1";
      str += (bit % 8 == 7) ? " " : "";
    }

    for (int bit = nValidDigits; bit < 32; bit++)
    {
      str += 'x';
      str += (bit % 8 == 7) ? " " : "";
    }

    return str;
  };// end of convertIntToBitString


  const Parameters& params = Parameters::getInstance();

  const Gene* bestIndividual = mReceiveIndividualBuffer + mBestIslandIdx * params.getChromosomeSize();

  std::string bestChromozome = "";

  const int nBlocks = globalKnapsackData->originalNumberOfItems / 32;

  for (int blockId = 0; blockId < nBlocks; blockId++)
  {
    bestChromozome += convertIntToBitString(bestIndividual[blockId], 32) + "\n";
  }

  // Reminder
  if (globalKnapsackData->originalNumberOfItems % 32 > 0 )
  {
    bestChromozome += convertIntToBitString(bestIndividual[nBlocks],
                                            globalKnapsackData->originalNumberOfItems % 32);
  }

 return bestChromozome;
} // end of getBestIndividualStr
//----------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------//
//----------------------------------------------- Protected methods --------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Allocate memory
 */
void Statistics::allocateMemory()
{
  const Parameters& params = Parameters::getInstance();

  // Allocate MPI buffers for stats data and best individual
  if (params.getIslandIdx() == 0)
  {
    mStatReceiveBuffer       = (StatsDataToExchange *) memalign(64,
                                                                sizeof(StatsDataToExchange) * params.getIslandCount());
    mReceiveIndividualBuffer = (Gene *)                memalign(64,
                                                                sizeof(Gene) * params.getChromosomeSize() *
                                                                params.getIslandCount());
  }
}// end of allocateMemory
//----------------------------------------------------------------------------------------------------------------------

/**
 * Free GPU memory
 */
void Statistics::freeMemory()
{
  if (mStatReceiveBuffer)
  {
    free(mStatReceiveBuffer);
    mStatReceiveBuffer = nullptr;
  }
  if (mReceiveIndividualBuffer)
  {
    free(mReceiveIndividualBuffer);
    mReceiveIndividualBuffer = nullptr;
  }

}// end of freeMemory
//----------------------------------------------------------------------------------------------------------------------


/**
 * Calculate local statistics
 */
void Statistics::calculateLocalStatistics(const Population* population)
{
  // initialize values
  mStatDataBuffer.maxFitness = Fitness(0);
  mStatDataBuffer.minFitness = Fitness(UINT_MAX);

  mStatDataBuffer.sumFitness = 0.0f;
  mStatDataBuffer.sum2Fitness= 0.0f;

  mLocalBestSolutionIdx  = 0;

  // Calculate statistics
  for (int individualIdx = 0; individualIdx < population->populationSize; individualIdx++)
  {
    if (mStatDataBuffer.maxFitness < population->fitness[individualIdx])
    {
      mStatDataBuffer.maxFitness = population->fitness[individualIdx];
      mLocalBestSolutionIdx = individualIdx;
    }

    mStatDataBuffer.minFitness = std::min(mStatDataBuffer.minFitness, population->fitness[individualIdx]);

    mStatDataBuffer.sumFitness  +=  population->fitness[individualIdx];
    mStatDataBuffer.sum2Fitness += (population->fitness[individualIdx]) * (population->fitness[individualIdx]);
  }
}// end of calculateLocalStatistics
//----------------------------------------------------------------------------------------------------------------------

/**
 * Calculate global statistics
 */
void Statistics::calculateGlobalStatistics()
{
  mBestIslandIdx = 0;

  // initialize values
  mStatDataBuffer.maxFitness  = mStatReceiveBuffer[0].maxFitness;
  mStatDataBuffer.minFitness  = mStatReceiveBuffer[0].minFitness;
  mStatDataBuffer.sumFitness  = mStatReceiveBuffer[0].sumFitness;
  mStatDataBuffer.sum2Fitness = mStatReceiveBuffer[0].sum2Fitness;

   const Parameters& params = Parameters::getInstance();

  // Calculate numeric stats over receiving buffer
  for (int islandIdx = 1;  islandIdx < params.getIslandCount(); islandIdx++)
  {
    if (mStatDataBuffer.maxFitness < mStatReceiveBuffer[islandIdx].maxFitness)
    {
      mStatDataBuffer.maxFitness = mStatReceiveBuffer[islandIdx].maxFitness;
      mBestIslandIdx = islandIdx;
    }

    mStatDataBuffer.minFitness = std::min(mStatDataBuffer.minFitness, mStatReceiveBuffer[islandIdx].minFitness);

    mStatDataBuffer.sumFitness  +=  mStatReceiveBuffer[islandIdx].sumFitness;
    mStatDataBuffer.sum2Fitness +=  mStatReceiveBuffer[islandIdx].sum2Fitness;
  } // for

 // Calculate derived statistics
 const int nIndividuals = params.getPopulationSize() * params.getIslandCount();
 mDerivedStats.avgFitness = mStatDataBuffer.sumFitness / nIndividuals;
 mDerivedStats.divergence = sqrt(fabs((mStatDataBuffer.sum2Fitness / nIndividuals) -
                                      (mDerivedStats.avgFitness * mDerivedStats.avgFitness)));
}// end of calculateGlobalStatistics
//----------------------------------------------------------------------------------------------------------------------

/**
 * Set ONLY best individual Idx and its fitness to statistics
 */
void  Statistics::setBestLocalIndividualAndMaxFintess(const Population* population)
{
  mStatDataBuffer.maxFitness = Fitness(0);

  mLocalBestSolutionIdx = 0;

  // Find the best island
  for (int individualIdx = 0; individualIdx < population->populationSize; individualIdx++)
  {
    if (mStatDataBuffer.maxFitness < population->fitness[individualIdx])
    {
      mStatDataBuffer.maxFitness = population->fitness[individualIdx];
      mLocalBestSolutionIdx      = individualIdx;
    }
  }
}// end of setBestLocalIndividualAndMaxFintess
//----------------------------------------------------------------------------------------------------------------------
