/* 
 * File:        GlobalKnapsackData.h
 * Author:      Jiri Jaros
 * Affiliation: Brno University of Technology
 *              Faculty of Information Technology
 *              
 *              and
 * 
 *              The Australian National University
 *              ANU College of Engineering & Computer Science
 *
 * Email:       jarosjir@fit.vutbr.cz
 * Web:         www.fit.vutbr.cz/~jarosjir
 * 
 * Comments:    Header file of the knapsack global data class. 
 *              This class maintains the benchmark data
 *
 * 
 * License:     This source code is distribute under OpenSource GNU GPL license
 *                
 *              If using this code, please consider citation of related papers
 *              at http://www.fit.vutbr.cz/~jarosjir/pubs.php        
 *      
 *
 * 
 * Created on 06 June 2012, 00:00 PM
 */


#ifndef GLOBALKNAPSACKDATA_H
#define	GLOBALKNAPSACKDATA_H




typedef int TPriceType;
typedef int TWeightType;

/*
 * GLobal data for Knapsack Benchmark
 */
class TGlobalKnapsackData{        
public:        
    int           NumberOfItems;                // Number of items in knapsack
    int           OriginalNumberOfItems;        // Original size without padding to multiple of 32
    int           KnapsackCapacity;             // Total knapsack capacity
    float         MaxPriceWightRatio;           // What is the best Price-Weight ration (to penalization)
    
    TPriceType  * ItemPrice;                    // An array listing all item prices
    TWeightType * ItemWeight;                   // An array listing all item weights 

    TGlobalKnapsackData();
    
    virtual ~TGlobalKnapsackData();
    
    void LoadFromFile();                        // Load data from file (uses TParameters to get the FileName)
    
    
protected:
        
    void AllocateMemory(int NumberOfItems);     // Allocate memory for arrays
    void FreeMemory();                          // Free memory of arrays

        
    
    
};// end of TGlobalKnapsackData
//------------------------------------------------------------------------------







#endif	/* GLOBALKNAPSACKDATA_H */

