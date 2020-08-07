#include <igraph/igraph.h>
#include <map>
#include <vector>
#include <set>
#include <cmath>
#include "stdio.h"

using namespace std;

pair<float, int> H(igraph_matrix_t *D,vector<int> **part, int N);
int H_Dyn(igraph_matrix_t *D, set<int>* candidates, int N, int source);
int cut_approx(igraph_matrix_t* D, int** confusion, int N, int theta);
int compute_MD(igraph_matrix_t *D, int gN);          
int compute_DynMD(igraph_matrix_t *D,int gN);
int compute_RMD(igraph_matrix_t *D, int gN, int rho, int theta);
void dist_sets(igraph_matrix_t *D, int N, FILE* ki,char* network,int ecount,double p,int cd, igraph_vector_t *deg);
