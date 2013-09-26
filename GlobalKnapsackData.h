/** 
 * @file        GlobalKnapsackData.h
 * @author      Jiri Jaros
 *   	 	Brno University of Technology \n
 *              Faculty of Information Technology \n
 *              
 *              and			\n
 * 
 *              The Australian National University	\n
 *              ANU College of Engineering & Computer Science	\n
 *
 * 		jarosjir@fit.vutbr.cz
 * 	        www.fit.vutbr.cz/~jarosjir
 * 
 *  
 * @brief       Header file of the knapsack global data class. 
 *              
 * @section     This source code is distribute under OpenSource GNU GPL license
 *                
 *              If using this code, please consider citation of related papers
 *              at http://www.fit.vutbr.cz/~jarosjir/pubs.php        
 *      
 *
 * @version	1.0
 * @date	06 June      2012, 00:00 (created)
 *		26 September 2013, 09:50 (revised)
 */


#ifndef GLOBALKNAPSACKDATA_H
#define	GLOBALKNAPSACKDATA_H



/**
 * @typedef  TPriceType
 * @brief Datatype for Knapsack item price
*/
typedef int TPriceType;

/**
 * @typedef  TWeightType
 * @brief    Datatype for Knapsack item weight
*/
typedef int TWeightType;


/**
 * @class TGlobalKnapsackData
 * @brief Global data for Knapsack Benchmark
 */
class TGlobalKnapsackData{        
public:        
    /// Number of items in knapsack
    int           NumberOfItems;        
    /// Original size without padding to multiple of 32        
    int           OriginalNumberOfItems;
    /// Total knapsack capacity
    int           KnapsackCapacity;     
    /// What is the best Price-Weight ration (to penalization)
    float         MaxPriceWightRatio;           
    
    /// An array listing all item prices
    TPriceType  * ItemPrice;                    
    /// An array listing all item weights 
    TWeightType * ItemWeight;                   

    /// Constructor
    TGlobalKnapsackData();
    /// Destructor
    virtual ~TGlobalKnapsackData();
    
    /// Load data from file (uses TParameters to get the FileName)
    void LoadFromFile();                        
    
    
protected:
    /// Allocate memory for arrays        
    void AllocateMemory(const int NumberOfItems);     
    /// Free memory of arrays
    void FreeMemory();                          

        
    
    
};// end of TGlobalKnapsackData
//------------------------------------------------------------------------------


#endif	/* GLOBALKNAPSACKDATA_H */

