## dynamic_source_localization_toolbox

The dynamic_source_localization_toolbox is a toolbox written in C++ to compute
the metric dimension (MD), the dynamic metric dimension (DynMD) of certain random graph models.
For details, please see the 

Odor, Gergely, and Patrick Thiran. "Sequential metric dimension for random graphs."

companion paper (there the DynMD is referred to as SMD).

### Installation:

Requirements:
- *igraph-0.7.1* or higher (C package)

In order to install the dynamic_source_localization_toolbox you install the igraph library
(potentially change the inlcude path in the header files), and then run the *make* command.
### Description:

```
NAME: 
MD_main -- compute the MD and/or the DynMD

SYNOPSIS: 
MD_main [parameters]  

DESCRIPTION:
The program generates a file in which each line has the basic properties of the generated
random graph and also the property of interest (MD, DynMD). The first line of each file
is a header which explains which column corresponds to which property.


PARAMETERS:
-net    network
        - Gnp: Erdos-Renyi network with N,p
        - Grid: 2D (NxN) grid, with random edges rewired with probability p
        - RGG: Random geometric graph, with random edges rewired with probability p
-N	graph size (number of nodes)
-Ns	graph size if a range of values is desired
        - pow2: powers of two from 2 to 1024
        - pow2_large: powers of 2 from 2048 to 32768
        - range10: multiples of 100 from 100 to 1000
-p	probability of edges (Gnp) or probability of rewiring (Grid, RGG)
        - 0: The single value 0
        - 1/2: The single value 0.5
        - 1: The single value 1
        - range20: Multiples of 0.05 from 0 to 1
        - range30: Multiples of 0.05 from 0 to 1 and multiples of 0.005 from 0 to 0.05
        - N^-1/4: The single value of N^(-1/4)
        - N^-1/2: The single value of N^(-1/2)
        - N^-2/3: The single value of N^(-2/3)
        - N^-3/4: The single value of N^(-3/4)
        - N^-4/5: The single value of N^(-4/5)
        - N^-5/6: The single value of N^(-5/6)
-prop   Property of interest
        - MD: Metric dimension (with greedy entropy approximation)
        - DynMD: Dynamic metric dimension (with greedy algroithm computed for each source)
        - MD+DynMD: Compute both the MD and the DynMD
-cd	connection distance for grid
        - default value is 1 (which is the usual grid)
-rad    radius of RGG
-i      number of iterations
        - default value is 1 
-o	output file
        - default value is stdout which will print on the screen
        - note: a suffix _iter%d.txt will be added to the output file name

OUTPUT:
The program generates a file specified by the -o parameter (or if no -o parameter is given on the screen).
The format of the file will be specified by the first line (header).

EXAMPLE:

# Generate an Erdos-Renyi graph with N=2,4,...,1024 and with p=N^(-1/2),
# and compute its Metrix dimension. Do this for 5 iterations and print it
# on the screen.

./MD_main -Ns pow2 -i 5 -p N^-2/3 -net Gnp -prop MD
```
