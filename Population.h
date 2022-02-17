/**
 * @file        Population.h
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
 * @brief       Header file of the GA population
 *
 * @date        06 June      2012, 00:00 (created)
 *              15 February  2022, 16:23 (revised)
 *
 * @copyright   Copyright (C) 2012 - 2022 Jiri Jaros.
 *
 * This source code is distribute under OpenSouce GNU GPL license.
 * If using this code, please consider citation of related papers
 * at http://www.fit.vutbr.cz/~jarosjir/pubs.php
 */

#ifndef POPULATION_H
#define	POPULATION_H

#include <string>

#include "GlobalKnapsackData.h"

/**
 * @typedef Gene
 * @brief Datatype for one gene.
 */
using Gene = unsigned int ;

/**
 * @typedef Fitness
 * @brief Datatype fitness value.
 */
using Fitness = float;


/**
 * @class Population
 * @brief CPU version of GA Population
 */
class Population
{
  public:
    /// Default constructor not allowed.
    Population() = delete;
    /// Constructor
    Population(const int PopulationSize, const int ChromosomeSize);
    /// Copy constructor not allowed.
    Population(const Population& orig) = delete;
    /// Destructor.
    virtual ~Population();

    /// Assignment operator not allowed.
    Population& operator=(const Population&) = delete;

    /**
     * Calculate fitness values of all chromosomes.
     * @param [in] globalKnapsackData - Knapsack data.
     */
    void calculateFitness(const GlobalKnapsackData& globalKnapsackData);
    
    /// Number of chromosomes.
    int populationSize;
    /// Size of chromosome in INTs.
    int chromosomeSize;

    /// 1D array of genes (chromosome-based encoding).
    Gene*    population;
    /// 1D array of fitness values.
    Fitness* fitness;

  protected:

    /// Memory allocation.
    void allocateMemory();

    /// Memory deallocation.
    void freeMemory();

};// end of Population
//----------------------------------------------------------------------------------------------------------------------

#endif	/* POPULATION_H */
