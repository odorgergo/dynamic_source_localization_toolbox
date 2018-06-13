#include <igraph/igraph.h>
#include <map>
#include <vector>
#include <set>
#include <cmath>
#include "stdio.h"

#include "MD_compute.h"
#include "MD_utils.h"

using namespace std;

int read_arguments(char network[], char property[], char outF[],char header[], vector<int>* Ns, vector<float>* ps, int*a, int* iter_max, int* cd, double* rad, int* rho_max, int* theta_max,int argc,char * const  argv[]);
