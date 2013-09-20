/* 
 * File:        CPU_Evolution.cpp
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
 * Comments:    Implementation file of the GA evolution
 *              This class controls the evolution process on multicore CPU
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


#include <iostream>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <mpi.h>
#include <sys/time.h>

#include "Random123/philox.h"

#include "CPU_Evolution.h"
#include "CPU_Statistics.h"
#include "Parameters.h"

using namespace std;
using namespace r123;


typedef r123::Philox4x32 RNG_4x32;
typedef r123::Philox2x32 RNG_2x32;



//----------------------------------------------------------------------------//
//                              Definitions                                   //
//----------------------------------------------------------------------------//


const int MPI_TAG_DATA_LEFT  = 100;
const int MPI_TAG_DATA_RIGHT = 101;

//----------------------------------------------------------------------------//
//                              Implementation                                //
//                              public methods                                //
//----------------------------------------------------------------------------//

/*
 * Constructor of the class
 * 
 * @param argc
 * @param argv
 */
TCPU_Evolution::TCPU_Evolution(int argc, char **argv){
        
    // Create parameters    
    Params  = TParameters::GetInstance();

    // Read parameters from command line
    Params->LoadParametersFromCommandLine(argc,argv);  
    
    // Load data from input file
    GlobalData.LoadFromFile();    
    Params->PrintAllParameters();
    
    // Create populations
    MasterPopulation    = new TCPU_Population(Params->PopulationSize(), Params->ChromosomeSize());    
    OffspringPopulation = new TCPU_Population(Params->OffspringPopulationSize(), Params->ChromosomeSize());
    
    EmigrantsToSend     = new TCPU_Population(Params->EmigrantCount(), Params->ChromosomeSize());
    EmigrantsToReceive  = new TCPU_Population(Params->EmigrantCount(), Params->ChromosomeSize());
    
    // create Statostics
    CPUStatistics       = new TCPU_Statistics();
    
    ActGeneration = 0;
  
    
}// end of TGPU_Evolution
//------------------------------------------------------------------------------
  

/*
 * Destructor of the class
 */
TCPU_Evolution::~TCPU_Evolution(){
    
    delete MasterPopulation;    
    delete OffspringPopulation;
            
    delete EmigrantsToSend;
    delete EmigrantsToReceive;
    
    delete CPUStatistics;
    
}// end of Destructor
//------------------------------------------------------------------------------

/*
 * Run Evolution
 */
void TCPU_Evolution::Run(){
    
    Initialize();
    
    RunEvolutionCycle();
    
}// end of Run
//------------------------------------------------------------------------------


//----------------------------------------------------------------------------//
//                              Implementation                                //
//                              protected methods                             //
//----------------------------------------------------------------------------//

/*
 * Initialization evolution
 */
void TCPU_Evolution::Initialize(){
            
   ActGeneration = 0;
      
   GenerateFirstPopulation(); 
       
   MasterPopulation->CalculateFitness_SSE(&GlobalData);
   
}// end of TGPU_Evolution
//------------------------------------------------------------------------------


/*
 * Get seed for random generator 
 */
r123_seed TCPU_Evolution::GetSeed() {                
      
  struct r123_seed seed;
  struct timeval tp1;  
  
  gettimeofday(&tp1, NULL);
  
  // use time in seconds divided by IslandID and 
  seed.seed1 = (unsigned long) (tp1.tv_sec) ;
  seed.seed2 = (unsigned long) tp1.tv_usec;
  
        
  return seed;
}// end of GetSeed
//------------------------------------------------------------------------------



/*
 * Generate first population
 */
void TCPU_Evolution::GenerateFirstPopulation(){
  
   const int PopulationDim = MasterPopulation->ChromosomeSize * MasterPopulation->PopulationSize;

   // Init Random generator
   RNG_4x32  rng_4x32;    
   RNG_4x32::key_type key    ={{Params->IslandIdx(),0xcaffe123}};
   RNG_4x32::ctr_type counter={{0, GetSeed().seed1 ,GetSeed().seed2 ,0xbeeff000}};
   RNG_4x32::ctr_type RandomValues;


    // Randomly init  genes, 4x unroll 
   for (int i = 0; i < (PopulationDim >> 2) << 2; i+=4){       
       counter.incr();
       RandomValues = rng_4x32(counter, key);    

       MasterPopulation->Population[i]   = RandomValues.v[0];
       MasterPopulation->Population[i+1] = RandomValues.v[1];
       MasterPopulation->Population[i+2] = RandomValues.v[2];
       MasterPopulation->Population[i+3] = RandomValues.v[3];              
   }


   
   // Init the rest 
    counter.incr();
    RandomValues = rng_4x32(counter, key);    

    int IDX = 0;
    for (int i = (PopulationDim >> 2) <<2; i < PopulationDim; i++){            
         MasterPopulation->Population[i]   = RandomValues.v[IDX];       
         IDX++;   
    }    



   // init fitness values
   for (int i = 0; i < MasterPopulation->PopulationSize; i++){       
       MasterPopulation->Fitness[i] = TFitness(0);
   } 

}// end of GenerateFirstPopulation
//------------------------------------------------------------------------------





/*
 * Run evolutionary cycle for defined number of generations
 * 
 */
void TCPU_Evolution::RunEvolutionCycle(){
    
    // Execute N generations
    for (ActGeneration = 1; ActGeneration < Params->NumOfGenerations(); ActGeneration++) {

        // If it's time for migration then migrate
        if (ActGeneration % Params->MigrationInterval() == 0) {            
            Migrate();
        }
        
        // Create new population
        GeneticManipulation();
        
        // Evaluate fitness
        OffspringPopulation->CalculateFitness_SSE(&GlobalData);       
        
        
        // Merge populations together
        Replacement();

        // if necessary, calculate and print statistics   
        if (ActGeneration % Params->StatisticsInterval() == 0){
              CPUStatistics->Calculate(MasterPopulation);
             
              if (IsMaster()) {
                printf("Generation %6d, BestIsland %d, MaxFitness %6f, MinFitness %6f, AvgFitness %6f, Diver %6f \n", 
                        ActGeneration, CPUStatistics->GetBestIslandIdx(),
                        CPUStatistics->GetMaxFitness(), CPUStatistics->GetMinFitness(),
                        CPUStatistics->GetAvgFitness(), CPUStatistics->GetDivergence());
              
                if (Params->GetPrintBest())  printf("%s\n", CPUStatistics->GetBestIndividualStr(&GlobalData).c_str());
              }  
          }
                  
    }// for ActGeneration

    
    //-- After evolution has finished --//                  
    CPUStatistics->Calculate(MasterPopulation);

    if (IsMaster()){
        printf("---------------------------------------------------------------------------------------------\n");
        printf("BestIsland %d, FinalMaxFitness %6f, FinalMinFitness %6f, FinalAvgFitness %6f, FinalDiver %6f \n", 
              CPUStatistics->GetBestIslandIdx(),
              CPUStatistics->GetMaxFitness(), CPUStatistics->GetMinFitness(),
              CPUStatistics->GetAvgFitness(), CPUStatistics->GetDivergence());
        
           
       printf("Best solution: \n");
       printf("%s\n", CPUStatistics->GetBestIndividualStr(&GlobalData).c_str());
    }

    
}// end of RunEvolutionCycle
//------------------------------------------------------------------------------
    
/*
 * Genetic manipulation
 * 
 */    
void TCPU_Evolution::GeneticManipulation(){
    
    
    const TParameters * Params = TParameters::GetInstance();
    
        
    //-- Init Random --//
    RNG_4x32  rng_4x32;    
    RNG_4x32::key_type key    ={{Params->IslandIdx(), 0}};
    RNG_4x32::ctr_type counter={{0, 0xabcd1234,GetSeed().seed1, GetSeed().seed2 }};
    RNG_4x32::ctr_type RandomValues;

    // go over master population
    for (int ChromosomeIdx = 0; ChromosomeIdx < Params->OffspringPopulationSize(); ChromosomeIdx+=2){

        //--------------------------------------------------------------------//
        //---------------------------- Selection -----------------------------//
        //--------------------------------------------------------------------//
        counter.incr();
        RandomValues = rng_4x32(counter, key);

        unsigned int Parent1_Idx = Selection(MasterPopulation, RandomValues.v[0], RandomValues.v[1]);                        
        unsigned int Parent2_Idx = Selection(MasterPopulation, RandomValues.v[2], RandomValues.v[3]);

        
        counter.incr();
        RandomValues = rng_4x32(counter, key);

        if ( RandomValues.v[0] < Params->CrossoverUINTBoundary()) {
        //--------------------------------------------------------------------//
        //------------------------ Crossover ---------------------------------//
        //--------------------------------------------------------------------//
            int GeneMod4Idx = 0;

            for ( int GeneIdx = 0; GeneIdx < Params->ChromosomeSize(); GeneIdx++ ){
                TGene GeneParent1 = MasterPopulation->Population[Parent1_Idx * Params->ChromosomeSize() + GeneIdx];
                TGene GeneParent2 = MasterPopulation->Population[Parent2_Idx * Params->ChromosomeSize() + GeneIdx];

                TGene GeneOffspring1 = 0;
                TGene GeneOffspring2 = 0;


                RNG_4x32::ctr_type RandomForMutation1;
                RNG_4x32::ctr_type RandomForMutation2;

                if (GeneMod4Idx == 4) {
                    GeneMod4Idx = 0;
                    counter.incr();
                    RandomValues = rng_4x32(counter, key);
                }
                 // crossover
                CrossoverUniformFlip(GeneOffspring1, GeneOffspring2, GeneParent1, GeneParent2, RandomValues.v[GeneMod4Idx]);
                GeneMod4Idx++;


                for (int BitID = 0; BitID < Params->IntBlockSize(); BitID+=4){                    

                    counter.incr();            
                    RandomForMutation1 = rng_4x32(counter, key);

                    counter.incr();
                    RandomForMutation2 = rng_4x32(counter, key);                        


                    for (int i = 0; i < 4 ; i++){                    
                        // mutation
                        MutationBitFlip(GeneOffspring1, GeneOffspring2, RandomForMutation1.v[i],RandomForMutation2.v[i], BitID+i);
                    }
                 }// mutation   

                 OffspringPopulation->Population[ChromosomeIdx     * Params->ChromosomeSize() + GeneIdx] = GeneOffspring1;
                 OffspringPopulation->Population[(ChromosomeIdx+1) * Params->ChromosomeSize() + GeneIdx] = GeneOffspring2;

              }
       //------------------------------------------------------------------------//
       //-------------------------- Cloning  ------------------------------------//
       //------------------------------------------------------------------------//
        }else {  //-- Go through two chromosomes and do uniform crossover --//
                for (int GeneIdx = 0; GeneIdx < Params->ChromosomeSize(); GeneIdx++ ){
                TGene GeneOffspring1 = MasterPopulation->Population[Parent1_Idx * Params->ChromosomeSize() + GeneIdx];
                TGene GeneOffspring2 = MasterPopulation->Population[Parent2_Idx * Params->ChromosomeSize() + GeneIdx];

                RNG_4x32::ctr_type RandomForMutation1;
                RNG_4x32::ctr_type RandomForMutation2;


                for (int BitID = 0; BitID < Params->IntBlockSize(); BitID+=4){

                    counter.incr();            
                    RandomForMutation1 = rng_4x32(counter, key);

                    counter.incr();
                    RandomForMutation2 = rng_4x32(counter, key);                        

                    for (int i = 0; i < 4 ; i++){
                        // mutation
                        MutationBitFlip(GeneOffspring1, GeneOffspring2, RandomForMutation1.v[i],RandomForMutation2.v[i], BitID+i);
                    }
                 }// one cloning
                OffspringPopulation->Population[ChromosomeIdx     * Params->ChromosomeSize() + GeneIdx] = GeneOffspring1;
                OffspringPopulation->Population[(ChromosomeIdx+1) * Params->ChromosomeSize() + GeneIdx] = GeneOffspring2;

            }// for 1 id            
        }      


        OffspringPopulation->Fitness[ChromosomeIdx]   = 0.0f;
        OffspringPopulation->Fitness[ChromosomeIdx+1] = 0.0f;

    }// for individuals        
   
}// end of GeneticManipulation
//------------------------------------------------------------------------------



/*
 * Binary tournament selection
 * 
 * @param ParentData - Parent Population
 * @param Random 1   - Frist Random value
 * @param Random 2   - Second Random value 
 * 
 * @return Selected chrromosome idx
 */
unsigned int TCPU_Evolution::Selection(TCPU_Population * ParentsData, unsigned int Random1, unsigned int Random2){
            
    unsigned int Idx1 = Random1 % (ParentsData->PopulationSize);
    unsigned int Idx2 = Random2 % (ParentsData->PopulationSize);
    
    return (ParentsData->Fitness[Idx1] > ParentsData->Fitness[Idx2]) ? Idx1 : Idx2;
}// Selection
//------------------------------------------------------------------------------


/*
 * Flip bites of parents to produce parents
 * 
 * @param  Offspring1
 * @param  Offspring2
 * @param  Parent1
 * @param  Parent2
 * @param  BitMask
 */
void TCPU_Evolution::CrossoverUniformFlip(TGene& GeneOffspring1, TGene& GeneOffspring2,
                                          TGene& GeneParent1   , TGene& GeneParent2,
                                          unsigned int RandomValue){

    GeneOffspring1 =  (~RandomValue  & GeneParent1) | ( RandomValue  & GeneParent2);
    GeneOffspring2 =  ( RandomValue  & GeneParent1) | (~RandomValue  & GeneParent2);
    
    
//    if (RandomValue > UINT_MAX>>1) { // flip bits           
//        GeneOffspring1 |=  (GeneParent2 & (1 << BitID));
//        GeneOffspring2 |=  (GeneParent1 & (1 << BitID));                            
//    } else {              
//        GeneOffspring1 |=  (GeneParent1 & (1 << BitID));
//        GeneOffspring2 |=  (GeneParent2 & (1 << BitID));                            
//    } 

}// end of CrossoverUniformFlip
//------------------------------------------------------------------------------



/*
 * Bit Flip Mutation
 * 
 * @param  Offspring1
 * @param  Offspring2
 * @param  Random Value 1
 * @param  Random Value 2
 * @param  Bit to mutate
 * 
 */
 void TCPU_Evolution::MutationBitFlip(TGene& GeneOffspring1, TGene& GeneOffspring2,
                                      unsigned int RandomValue1,unsigned int RandomValue2, int BitID){

 static TParameters * Params = TParameters::GetInstance();    
  
 
 if (RandomValue1 < Params->MutationUINTBoundary()) GeneOffspring1 ^= (1 << BitID); 
 if (RandomValue2 < Params->MutationUINTBoundary()) GeneOffspring2 ^= (1 << BitID); 
 
 //GeneOffspring1 ^= (((unsigned int) RandomValue1 < Params->MutationUINTBoundary()) << BitID);
 //GeneOffspring2 ^= (((unsigned int) RandomValue2 < Params->MutationUINTBoundary() ) << BitID);          
         
    
}// end of MutationBitFlip
//------------------------------------------------------------------------------
            
/*
 * Replacement kernel
 * 
 */
void TCPU_Evolution::Replacement(){
   
    TParameters * Params = TParameters::GetInstance();
    

    // Init Random Number Generator 
    RNG_2x32  rng_4x32;    
    RNG_2x32::key_type key    ={{Params->IslandIdx()}};
    RNG_2x32::ctr_type counter={{GetSeed().seed1 ,GetSeed().seed2}};
    RNG_2x32::ctr_type RandomValues;
    
    for (int ChromosomeIdx = 0; ChromosomeIdx < Params->PopulationSize(); ChromosomeIdx+=2){
        counter.incr();
        RandomValues = rng_4x32(counter, key);

        // Inline selection
        unsigned int Idx1 = RandomValues.v[0] % OffspringPopulation->PopulationSize;
        unsigned int Idx2 = RandomValues.v[1] % OffspringPopulation->PopulationSize;

        //  Repalace individuals if necessary
        if (OffspringPopulation->Fitness[Idx1] > MasterPopulation->Fitness[ChromosomeIdx]){

            memcpy(&MasterPopulation->Population   [ChromosomeIdx * Params->ChromosomeSize()], 
                   &OffspringPopulation->Population[Idx1          * Params->ChromosomeSize()],
                    Params->ChromosomeSize() * sizeof(TGene));

            MasterPopulation->Fitness[ChromosomeIdx] = OffspringPopulation->Fitness[Idx1];
        }

        if (OffspringPopulation->Fitness[Idx2] > MasterPopulation->Fitness[ChromosomeIdx+1]){

            memcpy(&MasterPopulation->Population   [(ChromosomeIdx+1) * Params->ChromosomeSize()], 
                   &OffspringPopulation->Population[Idx2              * Params->ChromosomeSize()],
                    Params->ChromosomeSize() * sizeof(TGene));
            MasterPopulation->Fitness[ChromosomeIdx+1] = OffspringPopulation->Fitness[Idx2];

        }
    }// ChromosomeIDx
  
}// end of Replacement()
//------------------------------------------------------------------------------

/*
 * Unidirectional migration
 */
void TCPU_Evolution::Migrate(){
  
  // Four messages running in parallel  
  static const int NumOfMessages = 4;  
  
  // MPI status and request
  MPI_Status  status [NumOfMessages];
  MPI_Request request[NumOfMessages];

  int Target;
  int Source;
 
  Target = (Params->IslandIdx() + 1) % Params->IslandCount(); // Send to the right
  if (Params->IslandIdx() == 0) Source = Params->IslandCount()-1;
  else Source = Params->IslandIdx() - 1;              // Send to the left
  
  // receive new fitnesses and chromosomes  
  MPI_Irecv(EmigrantsToReceive->Fitness   ,Params->EmigrantCount()                           , MPI_FLOAT   , Source, MPI_TAG_DATA_RIGHT, MPI_COMM_WORLD, &request[2]); 
  MPI_Irecv(EmigrantsToReceive->Population,Params->EmigrantCount() * Params->ChromosomeSize(), MPI_UNSIGNED, Source, MPI_TAG_DATA_RIGHT, MPI_COMM_WORLD, &request[3]); 
  

  PackEmigrantsToSend();
  
  // dispatch new emigrants
  MPI_Isend(EmigrantsToSend->Fitness      ,Params->EmigrantCount()                           , MPI_FLOAT   , Target, MPI_TAG_DATA_RIGHT, MPI_COMM_WORLD, &request[0]);         
  MPI_Isend(EmigrantsToSend->Population   ,Params->EmigrantCount() * Params->ChromosomeSize(), MPI_UNSIGNED, Target, MPI_TAG_DATA_RIGHT, MPI_COMM_WORLD, &request[1]);         

  MPI_Waitall(NumOfMessages, request, status);
  
  UnpackReceivedEmigrants();
    
  
}// end of Migrate
//------------------------------------------------------------------------------


/*
 * Create Population to send
 * 
 */
void TCPU_Evolution::PackEmigrantsToSend(){

    CPUStatistics->SetBestLocalIndividualAndMaxFintess(MasterPopulation);
    
    //-- Best individual --//
    EmigrantsToSend->Fitness[0] = CPUStatistics->GetMaxFitness();
    memcpy(&EmigrantsToSend->Population[0], 
           &MasterPopulation->Population[CPUStatistics->GetBestLocalSolutionIdx() * Params->ChromosomeSize()], 
            Params->ChromosomeSize() * sizeof(TGene));
    
    
    //-- other emigrants --//
     RNG_4x32  rng_4x32;    
     RNG_4x32::key_type key    ={{Params->IslandIdx(),0xcaffe123}};
     RNG_4x32::ctr_type counter={{0, GetSeed().seed1, GetSeed().seed2 ,0xbeeff00d}};
     RNG_4x32::ctr_type RandomValues;
    
    // Pack all emigrants into a buffer
    for (int EmigrantIdx = 1; EmigrantIdx < Params->EmigrantCount(); EmigrantIdx+=2){
        
       counter.incr();
       RandomValues = rng_4x32(counter, key);

       unsigned int Parent1_Idx = Selection(MasterPopulation, RandomValues.v[0], RandomValues.v[1]);                        
       unsigned int Parent2_Idx = Selection(MasterPopulation, RandomValues.v[2], RandomValues.v[3]); 
       
       EmigrantsToSend->Fitness[EmigrantIdx    ] = MasterPopulation->Fitness[Parent1_Idx];
       EmigrantsToSend->Fitness[EmigrantIdx + 1] = MasterPopulation->Fitness[Parent2_Idx];
       
       memcpy(&EmigrantsToSend->Population [EmigrantIdx * Params->ChromosomeSize()], 
              &MasterPopulation->Population[Parent1_Idx * Params->ChromosomeSize()], 
               Params->ChromosomeSize() * sizeof(TGene));
       
       memcpy(&EmigrantsToSend->Population [(EmigrantIdx + 1) * Params->ChromosomeSize()], 
              &MasterPopulation->Population[Parent2_Idx * Params->ChromosomeSize()], 
              Params->ChromosomeSize() * sizeof(TGene));
       
    }
    
}// end of PackEmigrantsToSend
//------------------------------------------------------------------------------
    

/*
 * Include Received population 
 * 
 */
void TCPU_Evolution::UnpackReceivedEmigrants(){
    
    
    
    //-- other emigrants --//
    RNG_4x32  rng_4x32;    
    RNG_4x32::key_type key    ={{Params->IslandIdx(),0xaacc8844}};
    RNG_4x32::ctr_type counter={{0, GetSeed().seed1, GetSeed().seed1 ,0xfeeff00d}};
    RNG_4x32::ctr_type RandomValues;
    
    // unpack emigrants    
    for (int EmigrantIdx = 0; EmigrantIdx < Params->EmigrantCount(); EmigrantIdx++){
        
      counter.incr();
      RandomValues = rng_4x32(counter, key);
      unsigned int ParentIdx = Selection(MasterPopulation, RandomValues.v[0], RandomValues.v[1]); 
       
      if (EmigrantsToReceive->Fitness[EmigrantIdx] > MasterPopulation->Fitness[ParentIdx]){
           MasterPopulation->Fitness[ParentIdx] = EmigrantsToReceive->Fitness[EmigrantIdx];
           memcpy(&MasterPopulation  ->Population[ParentIdx   * Params->ChromosomeSize()],
                  &EmigrantsToReceive->Population[EmigrantIdx * Params->ChromosomeSize()],                   
                   Params->ChromosomeSize() * sizeof(TGene));           
       }
    }  
       
     
}// end of UnpackReceivedEmigrants
//------------------------------------------------------------------------------