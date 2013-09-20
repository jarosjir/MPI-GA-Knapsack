MPI-GA-knapsack
===============

Distributed genetic algorithm solving the 0/1 knapsack problem (MPI version of the island model)


Outline

This package contains an efficient CPU implementation of the Island based GA running the
knapsack benchmark. The implementation is done in C++ and MPI. 

You can compile the code by typing:	make
If you want to run a simple demo, type:	make run

In order to specify the number of islands and their distribution over MPI nodes, you may have to
modify the makefile

For more information visit: http://www.fit.vutbr.cz/~jarosjir/pubs.php?id=9860&shortname=1
and read the content of 
Jaros, J.: Multi-GPU Island-Based Genetic Algorithm Solving the Knapsack Problem, 
In: 2012 IEEE World Congress on Computational Intelligence, CA, US, IEEE, 2012, p. 217-224, 
ISBN 978-1-4673-1508-1

CPU Requirements:
Intel Core i7 (SSE4.1) for SSE 
Intel Core 2  (NO SSE)



Software Requirements:
Compiler: g++-4.4 or newer
	  mpic++ 1.4.3 (1.5.4)