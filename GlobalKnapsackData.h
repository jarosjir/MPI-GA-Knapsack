/**
 * @file        GlobalKnapsackData.h
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
 * @brief       Header file of the knapsack global data class.
 *              This class maintains the benchmark data.
 *
 * @date        06 June      2012, 00:00 (created)
 *              15 February  2022, 12:14 (revised)
 *
 * @copyright   Copyright (C) 2012 - 2022 Jiri Jaros.
 *
 * This source code is distribute under OpenSouce GNU GPL license.
 * If using this code, please consider citation of related papers
 * at http://www.fit.vutbr.cz/~jarosjir/pubs.php
 *
 */

#ifndef GLOBAL_KNAPSACK_DATA_H
#define GLOBAL_KNAPSACK_DATA_H


/// Data type for item prices.
using PriceType = int ;
/// Data type for item weights.
using WeightType = int;

/**
 * @class GlobalKnapsackData
 * @brief Global data for Knapsack benchmark.
 */
class GlobalKnapsackData
{
  public:
    /// Constructor.
    GlobalKnapsackData();
    /// Disable copy constructor.
    GlobalKnapsackData(const GlobalKnapsackData&) = delete;
    /// Destructor.
    virtual ~GlobalKnapsackData();

    /// Disable assignment operator.
    GlobalKnapsackData& operator=(const GlobalKnapsackData&) = delete;

    /// Load data from file (uses TParameters to get the FileName)
    void loadFromFile();

    /// Number of items in knapsack padded to nearest multiply of 32.
    int         numberOfItems;
    /// Original size without padding to multiple of 32.
    int         originalNumberOfItems;
    /// Total knapsack capacity
    int         knapsackCapacity;
    /// What is the best Price-Weight ratio (to penalization)
    float       maxPriceWightRatio;

    /// Array listing all item prices.
    PriceType*  itemPrice;
    /// Array listing all item weights.
    WeightType* itemWeight;

  protected:
    /// Allocate memory for arrays.
    void allocateMemory(const int NumberOfItems);
    /// Free memory of arrays.
    void freeMemory();
};// end of GlobalKnapsackData
//----------------------------------------------------------------------------------------------------------------------

#endif	/* GLOBAL_KNAPSACK_DATA_H */
