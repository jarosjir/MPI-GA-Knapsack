# 
# File:        Makefile
# Author:      Jiri Jaros
# Affiliation: Brno University of Technology
#              Faculty of Information Technology
#              
#              and
# 
#              The Australian National University
#              ANU College of Engineering & Computer Science
#
# Email:       jarosjir@fit.vutbr.cz
# Web:         www.fit.vutbr.cz/~jarosjir
# 
# Comments:    Efficient MPI implementation of the island based Genetic Algorithm, 
#              solving the Knapsack problem.
#
# 
# License:     This source code is distribute under OpenSource GNU GPL license
#                
#              If using this code, please consider citation of related papers
#              at http://www.fit.vutbr.cz/~jarosjir/pubs.php        
#      
#
# 
# Created on 06 June 2012, 00:00 PM
#


# Environment
CC=mpicc
CXX=mpic++
CXXFLAGS=-m64 -O3 -msse4.1 -ffast-math -Wall
TARGET= mpi_ga_knapsack

all:		$(TARGET)	

$(TARGET):	main.o CPU_Statistics.o Parameters.o CPU_Population.o CPU_Evolution.o GlobalKnapsackData.o
	$(CXX) $(LDFLAGS) main.o CPU_Statistics.o Parameters.o CPU_Population.o CPU_Evolution.o GlobalKnapsackData.o -lm -o $@ $(LIBS) 



%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $<


# Clean Targets
clean: 
	/bin/rm -f *.o *.~ $(TARGET)

run:
#you are likely to be required to modify MPI running configuration 
	mpirun -np 4 ./mpi_ga_knapsack -f ./Data/knap_40.txt -p 100 -g 50 -s 10 -m 5

	


