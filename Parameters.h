/**
 * @file:       Parameters.h
 * @author	Jiri Jaros \n
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
 * @brief 	Header file the class that maintains all the parameters of evolution.
 *
 * 
 * @section	License
 *		This source code is distribute under OpenSource GNU GPL license
 *                
 *              If using this code, please consider citation of related papers
 *              at http://www.fit.vutbr.cz/~jarosjir/pubs.php        
 *      
 *
 * 
 * @version	1.0
 * @date	06 June      2012, 00:00 (created)
 *		20 September 2013, 14:58 (revised)
 */



#ifndef PARAMETERS_H
#define	PARAMETERS_H

#include <string>

using namespace std;


/**
 * @struct TEvolutionParameters
 * @brief  Structure with all parameters for the evolutionary algorithm
 */
struct TEvolutionParameters{    
    /// Population size
    int   PopulationSize;                
    /// Offspring population size
    int   OffspringPopulationSize;        
    /// Length of binary chromosome in int chunks
    int   ChromosomeSize;
    /// Total number of generations to evolve                 
    int   NumOfGenerations;               
    
    /// Crossover rate (as flaot)
    float CrossoverPst;                 
    /// Mutaion rate   (as float      
    float MutationPst;                  
    /// Crossover rate as uint 
    unsigned int CrossoverUINTBoundary; 
    /// Mutation rate as uint 
    unsigned int MutationUINTBoundary;  
    
    /// Number of migrating individuals    
    int EmigrantCount;
    /// Number of CPU threads                  
    int MigrationInterval;
    /// Index of Island              
    int IslandIdx;
    /// Number of independent islands                      
    int IslandCount;
    /// How often to print statistics                    
    int StatisticsInterval;             
    
    /// size of int block (32 bin genes)
    int IntBlockSize;                   
};// end of TEvolutionParameters
//------------------------------------------------------------------------------



/**
 * @class TParameters
 * @brief Singleton class with Parameters maintaining them in CPU and GPU constant memory
 */
class TParameters {
public:
    
    /// Get instance of the singleton class
    static TParameters* GetInstance();
    
    
    /// Destructor
    virtual ~TParameters() {
        pTParametersInstanceFlag = false;
    };
      
    /// Parse command line and populate the class
    void  LoadParametersFromCommandLine(int argc, char **argv);
      
    /// Get population size
    int   PopulationSize() 		const {return EvolutionParameters.PopulationSize; };
    /// Get chromosome length
    int   ChromosomeSize() 		const {return EvolutionParameters.ChromosomeSize; };        
    /// Set chromosome size
    void  SetChromosomeSize(const int Value) {EvolutionParameters.ChromosomeSize = Value; };    
    /// Get total number of generations
    int   NumOfGenerations() 		const {return EvolutionParameters.NumOfGenerations; };
    
    /// Get Crossover pst
    float CrossoverPst() 		const {return EvolutionParameters.CrossoverPst; };
    /// Get Mutation pst
    float MutationPst() 		const {return EvolutionParameters.MutationPst; };

    /// Get decision boundary for crossover as an uint
    unsigned int CrossoverUINTBoundary() const {return EvolutionParameters.CrossoverUINTBoundary; };    
    /// Get decision boundary formutation as an uint
    unsigned int MutationUINTBoundary() const {return EvolutionParameters.MutationUINTBoundary; };
    
    /// Get offspring population size
    int OffspringPopulationSize() 	const {return EvolutionParameters.OffspringPopulationSize; };
    /// Get number of emigants
    int EmigrantCount() 		const {return EvolutionParameters.EmigrantCount; };
    /// Get migration interval
    int MigrationInterval()    		const {return EvolutionParameters.MigrationInterval; };
    /// Get island index
    int IslandIdx()            		const {return EvolutionParameters.IslandIdx;};
    /// Get total number of islands
    int IslandCount()          		const {return EvolutionParameters.IslandCount; };    
    /// Get interval in which to collect statistics 
    int StatisticsInterval()   		const {return EvolutionParameters.StatisticsInterval;};    
        
    /// Get IntBlock size
    int IntBlockSize()         		const {return EvolutionParameters.IntBlockSize;};  
            
    /// Get filename with global data
    string BenchmarkFileName()    	const {return GlobalDataFileName;};
        
    /// Print best solution?
    bool GetPrintBest()         	const {return FPrintBest;};
    
    /// Print usage end exit
    void PrintUsageAndExit();
    /// print parameters to stdout
    void PrintAllParameters();
    
private:        
    /// Evolution  parameters
    TEvolutionParameters EvolutionParameters;   
    ///Global data filename
    string               GlobalDataFileName;
        
    /// Singleton instance
    static bool        	 pTParametersInstanceFlag;
    /// Singleton instance
    static TParameters *pTParametersSingle;
    /// Shall it print the best solution?    
    bool                FPrintBest;      
        
    /// prevent default constructor
    TParameters();

    ///Prevent copy-construction
    TParameters(const TParameters&);

    ///Prevent assignment
    TParameters& operator=(const TParameters&);
    
};// end of TParameters
//------------------------------------------------------------------------------


#endif	/* PARAMETERS_H */

