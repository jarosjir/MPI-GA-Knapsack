/**
 * @file        Evolution.cpp
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
 * @brief       Implementation file  of the island based GA
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

#include <limits.h>
#include <string.h>
#include <mpi.h>
#include <sys/time.h>

#include "Random123/philox.h"

#include "Evolution.h"
#include "Statistics.h"
#include "Parameters.h"

using namespace std;
using namespace r123;

/**
  * @typedef RNG_4x32
  * @brief four 32b values given by the Philox random generator.
  */
using RNG_4x32 = r123::Philox4x32 ;

/**
  * @typedef RNG_2x32
  * @brief two 32b values given by the Philox random generator.
  */
using RNG_2x32 = r123::Philox2x32;


//--------------------------------------------------------------------------------------------------------------------//
//--------------------------------------------------- Definitions ----------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------- Public methods ---------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Constructor of the class.
 */
Evolution::Evolution(int argc, char **argv)
  : mParams(Parameters::getInstance())
{
  // Read parameters from command line
  mParams.parseCommandline(argc,argv);

  // Load data from input file
  mGlobalData.loadFromFile();
  mParams.printAllParameters();

  // Create populations
  mMasterPopulation    = new Population(mParams.getPopulationSize(), mParams.getChromosomeSize());
  mOffspringPopulation = new Population(mParams.getOffspringPopulationSize(), mParams.getChromosomeSize());

  mEmigrantsToSend     = new Population(mParams.getEmigrantCount(), mParams.getChromosomeSize());
  mEmigrantsToReceive  = new Population(mParams.getEmigrantCount(), mParams.getChromosomeSize());

  // create Statistics
  mStatistics          = new Statistics();

  mActGeneration = 0;
}// end of Evolution
//----------------------------------------------------------------------------------------------------------------------

/**
 * Destructor of the class
 */
Evolution::~Evolution()
{
  delete mMasterPopulation;
  delete mOffspringPopulation;

  delete mEmigrantsToSend;
  delete mEmigrantsToReceive;

  delete mStatistics;
}// end of Destructor
//----------------------------------------------------------------------------------------------------------------------

/**
 * Run Evolution.
 */
void Evolution::run()
{
  initialize();

  runEvolutionCycle();
}// end of run
//----------------------------------------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------------------------------//
//----------------------------------------------- Protected methods --------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------//

/**
 * Get seed for random generator.
 * @return Seed for r123 generator.
 */
r123Seed Evolution::getSeed()
{
  struct r123Seed seed;
  struct timeval  tp1;

  gettimeofday(&tp1, nullptr);

  // use time in seconds divided by IslandID and
  seed.seed1 = static_cast<uint32_t>(tp1.tv_sec);
  seed.seed2 = static_cast<uint32_t>(tp1.tv_usec);

  /*seed.seed1 = static_cast<unsigned long>(0);
  seed.seed2 = static_cast<unsigned long>(0);*/

  return seed;
}// end of getSeed
//----------------------------------------------------------------------------------------------------------------------

/**
 * Initialization evolution.
 */
void Evolution::initialize()
{
   generateFirstPopulation();

   mMasterPopulation->calculateFitness(mGlobalData);
}// end of initialize
//----------------------------------------------------------------------------------------------------------------------

/**
 * Run evolutionary cycle for defined number of generations.
 */
void Evolution::runEvolutionCycle()
{
  // Execute N generations
  for (mActGeneration = 1; mActGeneration < mParams.getNumOfGenerations(); mActGeneration++)
  {
    // If it's time for migration then migrate
    if (mActGeneration % mParams.getMigrationInterval() == 0)
    {
      migrate();
    }

    // Create new population
    geneticManipulation();

    // Evaluate fitness
    mOffspringPopulation->calculateFitness(mGlobalData);

    // Merge populations together
    replacement();

    // If necessary, calculate and print statistics
    if (mActGeneration % mParams.getStatisticsInterval() == 0)
    {
      mStatistics->calculate(mMasterPopulation);

      if (isMaster())
      {
        printf("Generation %6d, BestIsland %d, MaxFitness %6f, MinFitness %6f, AvgFitness %6f, Diver %6f \n",
               mActGeneration, mStatistics->getBestIslandIdx(),
               mStatistics->getMaxFitness(), mStatistics->getMinFitness(),
               mStatistics->getAvgFitness(), mStatistics->getDivergence());

        if (mParams.getPrintBest())
        {
          printf("%s\n", mStatistics->getBestIndividualStr(&mGlobalData).c_str());
        }
      }
    } // print stats
  }// for actGeneration


  // After evolution has finished
  mStatistics->calculate(mMasterPopulation);

  if (isMaster())
  {
    printf("---------------------------------------------------------------------------------------------\n");
    printf("BestIsland %d, FinalMaxFitness %6f, FinalMinFitness %6f, FinalAvgFitness %6f, FinalDiver %6f \n",
           mStatistics->getBestIslandIdx(),
           mStatistics->getMaxFitness(), mStatistics->getMinFitness(),
           mStatistics->getAvgFitness(), mStatistics->getDivergence());


    printf("Best solution: \n");
    printf("%s\n", mStatistics->getBestIndividualStr(&mGlobalData).c_str());
  }
}// end of runEvolutionCycle
//----------------------------------------------------------------------------------------------------------------------

/**
 * Generate first population
 */
void Evolution::generateFirstPopulation()
{
  // Total number of genes
  const int nGenes = mMasterPopulation->chromosomeSize * mMasterPopulation->populationSize;

  // Init Random generator
  RNG_4x32  rng_4x32;
  RNG_4x32::key_type key     = {{uint32_t(mParams.getIslandIdx()), 0xcaffe123}};
  RNG_4x32::ctr_type counter = {{0, getSeed().seed1 ,getSeed().seed2 ,0xbeeff000}};
  RNG_4x32::ctr_type randomValues;


  // Randomly init  genes, 4x unroll
  for (int i = 0; i < (nGenes >> 2) << 2; i+=4)
  {
    counter.incr();
    randomValues = rng_4x32(counter, key);

    mMasterPopulation->population[i]     = randomValues.v[0];
    mMasterPopulation->population[i + 1] = randomValues.v[1];
    mMasterPopulation->population[i + 2] = randomValues.v[2];
    mMasterPopulation->population[i + 3] = randomValues.v[3];
   }

   // Init the rest
    counter.incr();
    randomValues = rng_4x32(counter, key);

    int idx = 0;
    for (int i = (nGenes >> 2) << 2; i < nGenes; i++)
    {
      mMasterPopulation->population[i] = randomValues.v[idx];
      idx++;
    }

   // init fitness values
   for (int i = 0; i < mMasterPopulation->populationSize; i++)
   {
      mMasterPopulation->fitness[i] = Fitness(0);
   }
}// end of generateFirstPopulation
//----------------------------------------------------------------------------------------------------------------------

/**
 * Genetic manipulation.
 */
void Evolution::geneticManipulation()
{
  // Init random generator
  RNG_4x32  rng_4x32;
  RNG_4x32::key_type key     = {{uint32_t(mParams.getIslandIdx()), 0}};
  RNG_4x32::ctr_type counter = {{0, 0xabcd1234, getSeed().seed1, getSeed().seed2 }};
  RNG_4x32::ctr_type randomValues;

  // Go over master population
  for (int chromosomeIdx = 0; chromosomeIdx < mParams.getOffspringPopulationSize(); chromosomeIdx += 2)
  {
    //------------------------------------------------- Selection ----------------------------------------------------//
    counter.incr();
    randomValues = rng_4x32(counter, key);

    unsigned int parent1Idx = selection(mMasterPopulation, randomValues.v[0], randomValues.v[1]);
    unsigned int parent2Idx = selection(mMasterPopulation, randomValues.v[2], randomValues.v[3]);


    counter.incr();
    randomValues = rng_4x32(counter, key);

   if ( randomValues.v[0] < mParams.getCrossoverUintBoundary())
   {
    //-------------------------------------------------- Crossover ---------------------------------------------------//
    int geneIdxMod4 = 0;

    for (int geneIdx = 0; geneIdx < mParams.getChromosomeSize(); geneIdx++ )
    {
      const Gene parent1Genes = mMasterPopulation->population[parent1Idx * mParams.getChromosomeSize() + geneIdx];
      const Gene parent2Genes = mMasterPopulation->population[parent2Idx * mParams.getChromosomeSize() + geneIdx];

      Gene offspring1Genes = 0;
      Gene offspring2Genes = 0;

      RNG_4x32::ctr_type randomValueForMutation1;
      RNG_4x32::ctr_type randomValueForMutation2;

      if (geneIdxMod4 == 4)
      {
        geneIdxMod4 = 0;
        counter.incr();
        randomValues = rng_4x32(counter, key);
      }

      // crossover
      crossoverUniformFlip(offspring1Genes,
                           offspring2Genes,
                           parent1Genes,
                           parent2Genes,
                           randomValues.v[geneIdxMod4]);
      geneIdxMod4++;

      //-------------------------------------------------- Mutation --------------------------------------------------//
      for (int bitIdx = 0; bitIdx < mParams.getIntBlockSize(); bitIdx += 4)
      {
        counter.incr();
        randomValueForMutation1 = rng_4x32(counter, key);

        counter.incr();
        randomValueForMutation2 = rng_4x32(counter, key);

        for (int i = 0; i < 4 ; i++)
        {
          // mutation
          mutationBitFlip(offspring1Genes,
                          offspring2Genes,
                          randomValueForMutation1.v[i],
                          randomValueForMutation2.v[i],
                          bitIdx + i);
          }
        }// mutation

        mOffspringPopulation->population[chromosomeIdx       * mParams.getChromosomeSize() + geneIdx] = offspring1Genes;
        mOffspringPopulation->population[(chromosomeIdx + 1) * mParams.getChromosomeSize() + geneIdx] = offspring2Genes;
      }
    }// crossover
    else
    //--------------------------------------------- Clone individuals ----------------------------------------------//
    {
      for (int GeneIdx = 0; GeneIdx < mParams.getChromosomeSize(); GeneIdx++)
      {
        Gene offspring1Genes = mMasterPopulation->population[parent1Idx * mParams.getChromosomeSize() + GeneIdx];
        Gene offspring2Genes = mMasterPopulation->population[parent2Idx * mParams.getChromosomeSize() + GeneIdx];

        RNG_4x32::ctr_type randomValueForMutation1;
        RNG_4x32::ctr_type randomValueForMutation2;

        //-------------------------------------------------- Mutation --------------------------------------------------//
        for (int bitIdx = 0; bitIdx < mParams.getIntBlockSize(); bitIdx+=4)
        {
          counter.incr();
          randomValueForMutation1 = rng_4x32(counter, key);

          counter.incr();
          randomValueForMutation2 = rng_4x32(counter, key);

          for (int i = 0; i < 4 ; i++)
          {
            // mutation
            mutationBitFlip(offspring1Genes,
                            offspring2Genes,
                            randomValueForMutation1.v[i],
                            randomValueForMutation2.v[i],
                            bitIdx + i);
          }
        }// one cloning

        mOffspringPopulation->population[chromosomeIdx       * mParams.getChromosomeSize() + GeneIdx] = offspring1Genes;
        mOffspringPopulation->population[(chromosomeIdx + 1) * mParams.getChromosomeSize() + GeneIdx] = offspring2Genes;
      }// for 1 id
    }// cloning

    // Reset fitness
    mOffspringPopulation->fitness[chromosomeIdx]     = 0.0f;
    mOffspringPopulation->fitness[chromosomeIdx +1 ] = 0.0f;
  }// for individuals
}// end of geneticManipulation
//----------------------------------------------------------------------------------------------------------------------

/**
 * Binary tournament selection.
 */
unsigned int Evolution::selection(const Population * parentPopulation,
                                  const unsigned int randomValue1,
                                  const unsigned int randomValue2)
{
  const unsigned int idx1 = randomValue1 % (parentPopulation->populationSize);
  const unsigned int idx2 = randomValue2 % (parentPopulation->populationSize);

  return (parentPopulation->fitness[idx1] > parentPopulation->fitness[idx2]) ? idx1 : idx2;
}// selection
//----------------------------------------------------------------------------------------------------------------------

/**
 * Flip bites of parents to produce offsprings.
 */
void Evolution::crossoverUniformFlip(Gene& offspring1,
                                     Gene& offspring2,
                                     const Gene& parent1,
                                     const Gene& parent2,
                                     const unsigned int randomValue)
{
  offspring1 =  (~randomValue & parent1) | ( randomValue & parent2);
  offspring2 =  ( randomValue & parent1) | (~randomValue & parent2);
}// end of crossoverUniformFlip
//----------------------------------------------------------------------------------------------------------------------

/**
 * Bit Flip Mutation.
 */
 void Evolution::mutationBitFlip(Gene& offspring1,
                                 Gene& offspring2,
                                 const unsigned int randomValue1,
                                 const unsigned int randomValue2,
                                 const int bitIdx)
{
  if (randomValue1 < mParams.getMutationUintBoundary()) offspring1 ^= (1 << bitIdx);
  if (randomValue2 < mParams.getMutationUintBoundary()) offspring2 ^= (1 << bitIdx);
}// end of mutationBitFlip
//----------------------------------------------------------------------------------------------------------------------

/**
 * Replacement kernel.
 */
 void Evolution::replacement()
{
  // Init Random Number Generator
  RNG_2x32  rng_4x32;
  RNG_2x32::key_type key     = {{uint32_t(mParams.getIslandIdx())}};
  RNG_2x32::ctr_type counter = {{getSeed().seed1 ,getSeed().seed2}};
  RNG_2x32::ctr_type randomValues;

  for (int chromosomeIdx = 0; chromosomeIdx < mParams.getOffspringPopulationSize(); chromosomeIdx += 2)
  {
    counter.incr();
    randomValues = rng_4x32(counter, key);

    // Inline selection
    const unsigned int idx1 = randomValues.v[0] % mMasterPopulation->populationSize;
    const unsigned int idx2 = randomValues.v[1] % mMasterPopulation->populationSize;

    //  Replace individuals if necessary
    if (mOffspringPopulation->fitness[chromosomeIdx] > mMasterPopulation->fitness[idx1])
    {
      memcpy(&mMasterPopulation->population   [idx1          * mParams.getChromosomeSize()],
             &mOffspringPopulation->population[chromosomeIdx * mParams.getChromosomeSize()],
             mParams.getChromosomeSize() * sizeof(Gene));

      mMasterPopulation->fitness[idx1] = mOffspringPopulation->fitness[chromosomeIdx];
    }

    if (mOffspringPopulation->fitness[chromosomeIdx + 1] > mMasterPopulation->fitness[idx2])
    {
      memcpy(&mMasterPopulation->population   [idx2                * mParams.getChromosomeSize()],
             &mOffspringPopulation->population[(chromosomeIdx + 1) * mParams.getChromosomeSize()],
             mParams.getChromosomeSize() * sizeof(Gene));
      mMasterPopulation->fitness[idx2] = mOffspringPopulation->fitness[chromosomeIdx + 1];

    }
  }// chromosomeIdx
}// end of replacement

/**
 * Unidirectional migration.
 */
void Evolution::migrate()
{
  // Four messages running in parallel
  constexpr int NumOfMessages = 4;

  // MPI status and request
  MPI_Status  status [NumOfMessages];
  MPI_Request request[NumOfMessages];

  // Get source and target
  const int target = (mParams.getIslandIdx() + 1) % mParams.getIslandCount();
  const int source = (mParams.getIslandIdx() == 0) ? mParams.getIslandCount() - 1 : mParams.getIslandIdx() - 1;

  // receive new fitnesses and chromosomes
  MPI_Irecv(mEmigrantsToReceive->fitness,
            mParams.getEmigrantCount(),
            MPI_FLOAT,
            source,
            kMpiTagDataRight,
            MPI_COMM_WORLD,
            &request[2]);
  MPI_Irecv(mEmigrantsToReceive->population,
            mParams.getEmigrantCount() * mParams.getChromosomeSize(),
            MPI_UNSIGNED,
            source,
            kMpiTagDataRight,
            MPI_COMM_WORLD,
            &request[3]);

  packEmigrantsToSend();

  // dispatch new emigrants
  MPI_Isend(mEmigrantsToSend->fitness,
            mParams.getEmigrantCount(),
            MPI_FLOAT,
            target,
            kMpiTagDataRight,
            MPI_COMM_WORLD,
            &request[0]);

  MPI_Isend(mEmigrantsToSend->population,
            mParams.getEmigrantCount() * mParams.getChromosomeSize(),
            MPI_UNSIGNED,
            target,
            kMpiTagDataRight,
            MPI_COMM_WORLD,
            &request[1]);

  MPI_Waitall(NumOfMessages, request, status);

  unpackReceivedEmigrants();
}// end of migrate
//----------------------------------------------------------------------------------------------------------------------

/**
 * Create Population to send
 */
void Evolution::packEmigrantsToSend()
{
  mStatistics->setBestLocalIndividualAndMaxFintess(mMasterPopulation);

  // Best individual
  mEmigrantsToSend->fitness[0] = mStatistics->getMaxFitness();
  memcpy(&mEmigrantsToSend->population[0],
         &mMasterPopulation->population[mStatistics->getBestLocalSolutionIdx() * mParams.getChromosomeSize()],
         mParams.getChromosomeSize() * sizeof(Gene));


  // other emigrants
  RNG_4x32  rng_4x32;
  RNG_4x32::key_type key     = {{uint32_t(mParams.getIslandIdx()), 0xcaffe123}};
  RNG_4x32::ctr_type counter = {{0, getSeed().seed1, getSeed().seed2, 0xbeeff00d}};
  RNG_4x32::ctr_type randomValues;

  // Pack all emigrants into a buffer
  for (int emigrantIdx = 1; emigrantIdx < mParams.getEmigrantCount(); emigrantIdx += 2)
  {
    counter.incr();
    randomValues = rng_4x32(counter, key);

    unsigned int parent1Idx = selection(mMasterPopulation, randomValues.v[0], randomValues.v[1]);
    unsigned int parent2Idx = selection(mMasterPopulation, randomValues.v[2], randomValues.v[3]);

    mEmigrantsToSend->fitness[emigrantIdx    ] = mMasterPopulation->fitness[parent1Idx];
    mEmigrantsToSend->fitness[emigrantIdx + 1] = mMasterPopulation->fitness[parent2Idx];

    memcpy(&mEmigrantsToSend->population [emigrantIdx * mParams.getChromosomeSize()],
           &mMasterPopulation->population[parent1Idx * mParams.getChromosomeSize()],
           mParams.getChromosomeSize() * sizeof(Gene));

    memcpy(&mEmigrantsToSend->population [(emigrantIdx + 1) * mParams.getChromosomeSize()],
           &mMasterPopulation->population[parent2Idx * mParams.getChromosomeSize()],
           mParams.getChromosomeSize() * sizeof(Gene));
  }
}// end of packEmigrantsToSend
//----------------------------------------------------------------------------------------------------------------------

/**
 * Include Received population.
 */
void Evolution::unpackReceivedEmigrants()
{
  //-- other emigrants --//
  RNG_4x32  rng_4x32;
  RNG_4x32::key_type key     = {{uint32_t(mParams.getIslandIdx()), 0xaacc8844}};
  RNG_4x32::ctr_type counter = {{0, getSeed().seed1, getSeed().seed1 ,0xfeeff00d}};
  RNG_4x32::ctr_type randomValues;

  // unpack emigrants
  for (int emigrantIdx = 0; emigrantIdx < mParams.getEmigrantCount(); emigrantIdx++)
  {

    counter.incr();
    randomValues = rng_4x32(counter, key);
    unsigned int parentIdx = selection(mMasterPopulation, randomValues.v[0], randomValues.v[1]);

    if (mEmigrantsToReceive->fitness[emigrantIdx] > mMasterPopulation->fitness[parentIdx])
    {
      mMasterPopulation->fitness[parentIdx] = mEmigrantsToReceive->fitness[emigrantIdx];
      memcpy(&mMasterPopulation  ->population[parentIdx   * mParams.getChromosomeSize()],
             &mEmigrantsToReceive->population[emigrantIdx * mParams.getChromosomeSize()],
             mParams.getChromosomeSize() * sizeof(Gene));
    }
  }
}// end of unpackReceivedEmigrants
//----------------------------------------------------------------------------------------------------------------------
