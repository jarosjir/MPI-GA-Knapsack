/**
 * @file        Statistics.h
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
 * @brief       Header file of the GA statistics over islands
 *
 * @date        06 June      2012, 00:00 (created)
 *              15 February  2022, 14:38 (revised)
 *
 * @copyright   Copyright (C) 2012 - 2022 Jiri Jaros.
 *
 * This source code is distribute under OpenSouce GNU GPL license.
 * If using this code, please consider citation of related papers
 * at http://www.fit.vutbr.cz/~jarosjir/pubs.php
 */

#ifndef STATISTICS_H
#define	STATISTICS_H


#include "Parameters.h"
#include "CPU_Population.h"
#include "GlobalKnapsackData.h"

/**
 * @struct StatsDataToExchange
 * @brief  Statistic structure for exchange between nodes.
 */
struct StatsDataToExchange
{
  /// Minimum fitness value in population.
  Fitness minFitness;
  /// Minimum fitness value in population.
  Fitness maxFitness;
  /// Sum of all fitness values over population.
  Fitness sumFitness;
  /// Sum of squares of all fitness values over population.
  Fitness sum2Fitness;
};// end of StatsDataToExchange
//----------------------------------------------------------------------------------------------------------------------

/**
 * @class DerivedStats
 * @brief Derived statistics structure.
 */
struct DerivedStats
{
  /// Average Fitness
  double avgFitness;
  /// Divergence
  double divergence;
};// end of DerivedStats
//----------------------------------------------------------------------------------------------------------------------

/**
 * @class Statistics
 * @brief Class that maintains local and global statistic of the island GA.
 */
class Statistics
{
  public:
    /// Constructor
    Statistics();
    /// copy constructor not allowed
    Statistics(const Population& orig) = delete;
    /// Destructor
    virtual ~Statistics();

    /// Disable assignment operator.
    Statistics& operator=(const Statistics&) = delete;


    /// Get maximum fitness value over all nodes, only master can read!
    Fitness getMaxFitness()    const { return mStatDataBuffer.maxFitness; };
    /// Get minimum fitness value over all nodes, only master can read!
    Fitness getMinFitness()    const { return mStatDataBuffer.minFitness; };
    /// Get average fitness value over all nodes, only master can read!
    double  getAvgFitness()    const { return mDerivedStats.avgFitness; };
    /// Get fitness divergence value over all nodes, only master can read!
    double  getDivergence()    const { return mDerivedStats.divergence; };
    /// Get the ID of the best island over all nodes, only master can read!
    int     getBestIslandIdx() const { return mBestIslandIdx; };
    /// Get the best local solution index.
    int     getBestLocalSolutionIdx() const {return mLocalBestSolutionIdx;};

    /**
     * Get best individual chromosome as a string, only master can call this function!
     * @param [in] globalKnapsackData
     * @return Best individual in from of a sting
     */
    std::string getBestIndividualStr(const GlobalKnapsackData* globalKnapsackData);

    /**
     * Calculate statistics - all islands must call this routine.
     * @param [in] population - Population to calculate statistics at.
     */
    void    calculate(const Population* population);

    /**
     * Set best individual and max fitness for the local node
     * @param [in] population - Population to find the best solution and set max index
     */
    void    setBestLocalIndividualAndMaxFintess(const Population* population);

  protected:
    /// Allocate Memory.
    void allocateMemory();
    /// Free Memory.
    void freeMemory();

    /**
     * Calculate local stats.
     * @param [in] population - Population to calculate statistics at.
     */
    void calculateLocalStatistics(const Population* population);
    /// Calculate global stats.
    void calculateGlobalStatistics();

    /// Statistics data to exchange.
    StatsDataToExchange  mStatDataBuffer;
    /// Buffer to receive statistics.
    StatsDataToExchange* mStatReceiveBuffer;

    /// Derived statistics.
    DerivedStats mDerivedStats;
    /// Index to the local best solution.
    int          mLocalBestSolutionIdx;
    /// Best island (island with the best solution).
    int          mBestIslandIdx;

    /// Buffer for best solution coming over network.
    Gene*        mReceiveIndividualBuffer;
};// end of Statistics
//----------------------------------------------------------------------------------------------------------------------

#endif	/* STATISTICS_H */
