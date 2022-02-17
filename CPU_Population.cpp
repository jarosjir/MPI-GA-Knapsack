/**
 * @file        Population.cpp
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
 * @brief       Implementation file of the GA population.
 *
 * @date        06 June      2012, 00:00 (created)
 *              15 February  2022, 16:37 (revised)
 *
 * @copyright   Copyright (C) 2012 - 2022 Jiri Jaros.
 *
 * This source code is distribute under OpenSouce GNU GPL license.
 * If using this code, please consider citation of related papers
 * at http://www.fit.vutbr.cz/~jarosjir/pubs.php
 */

#include <sstream>
#include <smmintrin.h>
#include <malloc.h>
#include <string.h>

#include "CPU_Population.h"
#include "Parameters.h"


//--------------------------------------------------------------------------------------------------------------------//
//--------------------------------------------------- Definitions ----------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Public methods ---------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Constructor of the class.
 *
 * @param [in] populationSize
 * @param [in] chromosomeSize
 *
 */
Population::Population(const int populationSize,
                       const int chromosomeSize)
  : populationSize(populationSize),
    chromosomeSize(chromosomeSize)
{
  allocateMemory();
}// end of Population
//----------------------------------------------------------------------------------------------------------------------

/**
 * Destructor of the class
 */
Population::~Population()
{
  freeMemory();
}// end of Population
//----------------------------------------------------------------------------------------------------------------------

/**
 * Calculate fitness, normal version
 */
void Population::calculateFitness(const GlobalKnapsackData& globalKnapsackData)
{
  const Parameters& params = Parameters::getInstance();

  constexpr unsigned int ONE = 1;

  for (int chromosomeIdx = 0; chromosomeIdx < populationSize; chromosomeIdx++)
  {
    const int startIdx = chromosomeIdx * chromosomeSize;

    int weight[4] = {0, 0, 0, 0};
    int price [4] = {0, 0, 0, 0};
    int value [4] = {0, 0, 0, 0};

    // Calculate price and weight, loop is unrolled
    for (int blockIdx = 0; blockIdx < chromosomeSize; blockIdx++)
    {
      for (int i = 0; i < params.getIntBlockSize(); i += 4)
      {
        size_t index = blockIdx * params.getIntBlockSize() + i;

        value[0] = (population[startIdx + blockIdx] >> (i)    ) & ONE;
        value[1] = (population[startIdx + blockIdx] >> (i + 1)) & ONE;
        value[2] = (population[startIdx + blockIdx] >> (i + 2)) & ONE;
        value[3] = (population[startIdx + blockIdx] >> (i + 3)) & ONE;

        weight[0] += value[0] * globalKnapsackData.itemWeight[index];
        price [0] += value[0] * globalKnapsackData.itemPrice [index];

        weight[1] += value[1] * globalKnapsackData.itemWeight[index + 1];
        price [1] += value[1] * globalKnapsackData.itemPrice [index + 1];

        weight[2] += value[2] * globalKnapsackData.itemWeight[index + 2];
        price [2] += value[2] * globalKnapsackData.itemPrice [index + 2];

        weight[3] += value[3] * globalKnapsackData.itemWeight[index + 3];
        price [3] += value[3] * globalKnapsackData.itemPrice [index + 3];
      }
    }

    const int totalPrice  = (price [0] + price [1]) + (price [2] + price [3]);
    const int totalWeight = (weight[0] + weight[1]) + (weight[2] + weight[3]);

    //-- penalization --//
    Fitness result = (float) totalPrice;

    if (totalWeight > globalKnapsackData.knapsackCapacity)
    {
      unsigned int penalty = totalWeight - globalKnapsackData.knapsackCapacity;

      result = result - globalKnapsackData.maxPriceWightRatio * penalty;
    }

    fitness[chromosomeIdx] = (result < 0) ? Fitness(0) : result;
  }// end master for

}// end of calculateFitness
//----------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------//
//----------------------------------------------- Protected methods --------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Allocate memory.
 */
void Population::allocateMemory()
{
  population = (Gene    *) memalign(64, sizeof(Gene)    * chromosomeSize * populationSize);
  fitness    = (Fitness *) memalign(64, sizeof(Fitness) * populationSize);
}// end of allocateMemory
//----------------------------------------------------------------------------------------------------------------------

/**
 * Free memory
 */
void Population::freeMemory()
{
  free(population);
  free(fitness);
}// end of freeMemory
//----------------------------------------------------------------------------------------------------------------------
