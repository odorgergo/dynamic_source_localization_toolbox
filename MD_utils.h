#include <igraph/igraph.h>
#include <map>
#include <vector>
#include <set>
#include <cmath>
#include "stdio.h"

using namespace std;

/*bool testMD(int N, vector<int> sensors);
bool testDMD(int N, vector<int> sensors);
void printD(int N);*/
void free_complist(igraph_vector_ptr_t *complist);
void largest_connected(const igraph_t *graph, igraph_t *giant);