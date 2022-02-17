/**
 * @file        Parameters.h
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
 * @brief       Header file of the parameter class. This class maintains all the parameters of evolution.
 *
 * @date        06 June      2012, 00:00 (created)
 *              15 February  2022, 11:11 (revised)
 *
 * @copyright   Copyright (C) 2012 - 2022 Jiri Jaros.
 *
 * This source code is distribute under OpenSouce GNU GPL license.
 * If using this code, please consider citation of related papers
 * at http://www.fit.vutbr.cz/~jarosjir/pubs.php
 */

#ifndef PARAMETERS_H
#define	PARAMETERS_H

#include <string>

/**
 * @struct EvolutionParameters
 * @brief  Parameters of the evolutionary process.
 */
struct EvolutionParameters
{
  /// Population size - number of chromosomes in the population.
  int   populationSize;
    /// Offspring population size - number of newly generated chromosomes.
  int   offspringPopulationSize;

  /// Length of binary chromosome in integer granularity (32b).
  int   chromosomeSize;
  /// Total number of generations to evolve.
  int   numOfGenerations;

  /// Crossover probability per individual (float number).
  float crossoverPst;
  /// Mutation probability per gene (float number).
  float mutationPst;
  /// Crossover rate converted to int for faster comparison and random number generation.
  unsigned int crossoverUintBoundary;
  /// Mutation rate converted to int for faster comparison and random number generation.
  unsigned int mutationUintBoundary;

  /// Number of migrating individuals between islands.
  int emigrantCount;
  /// Migration interval (how often to migrate).
  int migrationInterval;
  /// Index of the island.
  int islandIdx;
  /// Number of independent islands.
  int islandCount;
  /// How often to print statistics
  int statisticsInterval;

  /// size of int block (32 bin genes)
  int intBlockSize;
};// end of EvolutionParameters
//----------------------------------------------------------------------------------------------------------------------

/**
 * @class Parameters
 * @brief Singleton class with Parameters maintaining them in CPU and GPU constant memory.
 */
class Parameters
{
  public:
    /// Get instance of the singleton class.
    static Parameters& getInstance();

    // /Prevent copy-construction.
    Parameters(const Parameters&) = delete;
    ///Prevent assignment.
    Parameters& operator=(const Parameters&) = delete;

    /// Destructor
    virtual ~Parameters() { sInstanceFlag = false; };

    /// Parse command line and populate the class.
    void  parseCommandline(int argc, char **argv);

    //--------------------------------------------------- Getters ----------------------------------------------------//

    /// Get number of chromosomes in the population
    int   getPopulationSize()               const { return mEvolutionParameters.populationSize; };
    /// Get size of the chromosome (including padding)
    int   getChromosomeSize()               const { return mEvolutionParameters.chromosomeSize; };
    /// Set size of the chromosome
    void  setChromosomeSize(unsigned int value)   { mEvolutionParameters.chromosomeSize = value; };
    /// Get number of generations to evolve
    int   getNumOfGenerations()             const { return mEvolutionParameters.numOfGenerations; };

    /// Get crossover probability for two individuals.
    float        getCrossoverPst()          const { return mEvolutionParameters.crossoverPst; };
    /// Get per gene mutation probability.
    float        getMutationPst()           const { return mEvolutionParameters.mutationPst; };
    /// Get crossover probability in scaled to uint.
    unsigned int getCrossoverUintBoundary() const { return mEvolutionParameters.crossoverUintBoundary; };
    /// Get mutation probability in scaled to uint.
    unsigned int getMutationUintBoundary()  const { return mEvolutionParameters.mutationUintBoundary; };

    /// Get offspring population size.
    int   getOffspringPopulationSize()      const {return mEvolutionParameters.offspringPopulationSize; };
    /// Get number of emigrants.
    int   getEmigrantCount() 	            	const {return mEvolutionParameters.emigrantCount; };
    /// Get migration interval.
    int   getMigrationInterval()    		    const {return mEvolutionParameters.migrationInterval; };
    /// Get island index.
    int   getIslandIdx()            		    const {return mEvolutionParameters.islandIdx;};
    /// Get total number of islands.
    int   getIslandCount()          		    const {return mEvolutionParameters.islandCount; };
    /// Get how often to print statistics.
    int  getStatisticsInterval()            const { return mEvolutionParameters.statisticsInterval; };
    /// Get the integer block size.
    int  getIntBlockSize()                  const { return mEvolutionParameters.intBlockSize; };

    /// Get filename with global data.
    std::string getBenchmarkFileName()    	const {return mGlobalDataFileName;};

    /// Print best solution?
    bool getPrintBest()                    	const {return mPrintBest;};

    /// Print usage end exit
    void printUsageAndExit();
    /// print parameters to stdout.
    void printAllParameters();

  private:
    /// Singleton constructor.
    Parameters();

    /// Evolution  parameters
    EvolutionParameters mEvolutionParameters;
    /// Global data filename
    std::string         mGlobalDataFileName;

    /// Shall it print the best solution?
    bool                mPrintBest;

    /// Singleton instance
    static bool        	sInstanceFlag;
    /// Singleton instance
    static Parameters*  sSingletonInstance;
};// end of Parameters
//----------------------------------------------------------------------------------------------------------------------

#endif	/* PARAMETERS_H */

