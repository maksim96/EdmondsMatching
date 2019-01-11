# Edmonds' Matching Algorithm
This is an efficient implementation of the (unweighted) Edmonds' matching algorithm as described in [Korte, Vygen Combinatorial Optimization] for the course "Combinatorial Optimization" held in the winter semester 2018/19 at the Univerity of Bonn.

## Compile and Run
Compile with: compile.sh

Run with: ./EdmondsMatching file.dmx or ./all_files.sh

## Data Format
The algorithm expects the input graph to be given as an .dmx file (check bin/graphs) and outputs the matching also in this format.

bin/graphs also contains an optima.md file with the correct matching sizes for all the graphs.

## Implementation Details
We first perform a greedy initialization and then compute the matching as described in [Korte, Vygen Combinatorial Optimization].

Our rho is only implicitly there. It doesn't contain the real rho values, but is only a pointer to the "next" contracted blossom root. So, if you want to have the real rho value you have to traverse all the way to the top. Because of this, we could get rid of the loop, traversing all nodes in V(G) in the end of step (6), and just traverse over the paths, which heavily speedup the algorithm.  
