/**
 * @file        Evolution.h
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
 * @brief       Header file  of the island based GA
 *
 * @date        06 June      2012, 00:00 (created)
 *              15 February  2022, 16:51 (revised)
 *
 * @copyright   Copyright (C) 2012 - 2022 Jiri Jaros.
 *
 * This source code is distribute under OpenSouce GNU GPL license.
 * If using this code, please consider citation of related papers
 * at http://www.fit.vutbr.cz/~jarosjir/pubs.php
 */

#ifndef EVOLUTION_H
#define	EVOLUTION_H


#include "Parameters.h"
#include "CPU_Population.h"
#include "CPU_Statistics.h"
#include "GlobalKnapsackData.h"


/**
 * @struct r123Seed
 * @brief struct of seed
 */
struct r123Seed
{
  /// first seed
  uint32_t seed1;
  /// second seed
  uint32_t seed2;
};// end of r123_seed
//----------------------------------------------------------------------------------------------------------------------

/**
 * @class Evolution
 * @brief CPU evolution process.
 *
 */
class Evolution
{
  public:
    /// Disable default constructor
    Evolution() = delete;
    /// Constructor
    Evolution(int argc, char **argv);
        /// Copy constructor
    Evolution(const Evolution&) = delete;
    /// Destructor
    virtual ~Evolution();

    /// Disable assignment operator.
    Evolution& operator=(const Evolution&) = delete;


    /// Run evolution
    void run();
    /// Is this a master process?
    bool isMaster() const { return mParams.getIslandIdx() == 0; };

  protected:
    /// Get Random Seed for Random123.
    r123Seed getSeed();


    /// Initialize.
    void initialize();
    /// Run evolution cycle.
    void runEvolutionCycle();

    /// Generate first population.
    void generateFirstPopulation();

    /// Genetic manipulation.
    void geneticManipulation();

    /// Replacement phase.
    void replacement();

    /// Migrate phase.
    void migrate();

    /// Pack Emigrants to a buffer for dispatch.
    void packEmigrantsToSend();

    /// Unpack immigrants form buffer.
    void unpackReceivedEmigrants();

    /**
     * Selection of two chromosomes.
     * @param [in] parentPopulation - Parent population.
     * @param [in] randomValue1     - First random number.
     * @param [in] randomValue2     - Second random number.
     * @return Selected chromosome idx.
     */
    inline unsigned int selection(const Population*  parentPopulation,
                                  const unsigned int randomValue1,
                                  const unsigned int randomValue2);

    /**
     * Uniform Crossover.
     * @param [out] offspring1  - First generated offspring.
     * @param [out] offspring2  - Second generated offspring.
     * @param [in]  parent1     - First parent
     * @param [in]  parent2     - Second parent
     * @param [in]  randomValue - Random value
     */
    inline void crossoverUniformFlip(Gene&              offspring1,
                                     Gene&              offspring2,
                                     const Gene&        parent1,
                                     const Gene&        parent2,
                                     const unsigned int randomValue);

    /**
     * Bit Flip Mutation.
     * @param [out] offspring1   - First generated offspring.
     * @param [out] offspring2   - Second generated offspring.
     * @param [in]  randomValue1 - First random number.
     * @param [in]  randomValue2 - Second random number.
     * @param [in]  bitIdx       - Which bit to mutate.
     */
    inline void mutationBitFlip(Gene&              offspring1,
                                Gene&              offspring2,
                                const unsigned int randomValue1,
                                const unsigned int randomValue2,
                                const int          bitIdx);

    /// Parameters of evolution.
    Parameters& mParams;
    /// Actual generation.
    int         mActGeneration;

    /// Master GA population.
    Population* mMasterPopulation;
    /// Population of offsprings.
    Population* mOffspringPopulation;

    /// Population of emigrants to send.
    Population* mEmigrantsToSend;
    /// Population of immigrants to receive.
    Population* mEmigrantsToReceive;

    /// Statistics over GA process
    Statistics*        mStatistics;
    /// Global data of knapsack
    GlobalKnapsackData mGlobalData;

    /// MPI tag for sending data on the ring to the left
    static constexpr int kMpiTagDataLeft  = 100;
    /// MPI tag for sending data on the ring to the right
    static constexpr int kMpiTagDataRight = 101;

};// end of Evolution
//----------------------------------------------------------------------------------------------------------------------

#endif	/* EVOLUTION_H */

