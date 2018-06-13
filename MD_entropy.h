#include <igraph/igraph.h>
#include <map>
#include <vector>
#include <cmath>
#include "stdio.h"

using namespace std;

vector<int> **part;
igraph_matrix_t D;
int** confusion;

bool testMD(int N, vector<int> sensors) {
    for(int v=0; v<N; v++) {
      for(int w=v+1; w<N; w++) {
        bool foundS=false;
         for(vector<int>::iterator it = sensors.begin(); it != sensors.end(); ++it) {
           if (MATRIX(D,v,*it)!=MATRIX(D,w,*it)) {
             foundS=true;
             break;
           }
         }
       if (!foundS) {
         printf("%d %d %d %d %d\n",v,w,(int)MATRIX(D,v,sensors[0]),(int)MATRIX(D,w,sensors[0]),sensors[0]);
         return false;
       }
     }
    }
    return true;
}

bool testDMD(int N, vector<int> sensors) {
    for(int v=0; v<N; v++) {
      for(int w=v+1; w<N; w++) {
        bool foundS=false;
         for(vector<int>::iterator it = sensors.begin(); it != sensors.end(); ++it) {
           for(vector<int>::iterator it2 = sensors.begin(); it2 != sensors.end(); ++it2) {
             if ((MATRIX(D,v,*it)-MATRIX(D,w,*it)) !=  (MATRIX(D,v,*it2)-MATRIX(D,w,*it2))) {
               foundS=true;
               break;
             }
            }
         }
       if (!foundS) {
         printf("%d %d %d %d\n",v,w,(int)MATRIX(D,v,sensors[0]),(int)MATRIX(D,w,sensors[0]));
         return false;
       }
     }
    }
    return true;
}


pair<float, int> H(int N) {
    float best_res=100000000;
    int best_guess=-1;
    
    for(int v=0; v<N; v++) {
      //TODO Here we could check not to chose a v that's already added. 
      //Compute the size of the refined part with node v
      map<pair<int,vector<int> >,int> buckets;
	  for(int i=0; i<N; i++) {
	     pair<int,vector<int> > guess_pair = make_pair(MATRIX(D,v,i), *part[i]);
	     if (buckets.find(guess_pair)==buckets.end() )
	       buckets.insert(make_pair(guess_pair,1));
	     else
	       buckets.at(guess_pair)++;
	  }
	  
	  // Compute entropy
	  float res=0;
	  for(map<pair<int,vector<int> >,int>::iterator it=buckets.begin(); it!=buckets.end(); ++it) {
	    for(int i = 0; i<it->second; i++) {
	      res+=log(i+1);
	    }
	  }
	  
	  //Save best
	  if (res<best_res) {
	    best_res=res;
	    best_guess=v;
	  }
	  
	}
	
	// Add the best partition to the partition data structure
	for(int i=0; i<N; i++) {
	  part[i]->push_back(MATRIX(D,best_guess,i));
	}
	
	return make_pair(best_res,best_guess);
}


pair<float, int> cut_approx(int N, int theta) {
    float best_count=-1;
    int best_guess=-1;
    
    for(int v=0; v<N; v++) {
      //TODO Here we could check not to chose a v that's already added. 
      //Compute the number of edges cut by node v
      int cut_count=0;
      
	  for(int i=0; i<N; i++) {
	    for(int j=i+1; j<N; j++) {
	      if (confusion[i][j]==1) 
	        if (abs(MATRIX(D,v,i)-MATRIX(D,v,j))>theta)
	          cut_count++;
	    }
	  }
	  
	  
	  //Save best
	  if (cut_count>best_count) {
	    best_count=cut_count;
	    best_guess=v;
	  }
	  
	}
	
	int still_left=0;
	// Cut the edges
	for(int i=0; i<N; i++) {
	  for(int j=i+1; j<N; j++) {
	     if (confusion[i][j]==1) {
	       if (abs(MATRIX(D,best_guess,i)-MATRIX(D,best_guess,j))>theta) {
	          confusion[i][j]=0;
	       }
	       else {
	          still_left++;
	       }
	     }
	  }
	}
	
	return make_pair(still_left,best_guess);
}

void printD(int N) {
   for(int i=0; i<N;i++) {
     printf("i=%d: ", i);
     for(int j=0; j<N;j++) {
       printf("%d ", (int) MATRIX(D,i,j));
     }
     printf("\n");
   }
}


void free_complist(igraph_vector_ptr_t *complist) {
  long int i;
  for (i=0; i<igraph_vector_ptr_size(complist); i++) {
    igraph_t *comp=(igraph_t *) VECTOR(*complist)[i]; 
    igraph_destroy(comp);
    free(VECTOR(*complist)[i]);
  }
}

void largest_connected(const igraph_t *graph, igraph_t *giant) {
     igraph_vector_ptr_t complist;
     igraph_vector_ptr_init(&complist, 0);
     igraph_decompose(graph, &complist, IGRAPH_WEAK, -1, 0);  
     int lcomp_size=-1;
     int lcomp_i=-1;
     for (int lci=0; lci< igraph_vector_ptr_size(&complist); lci++) {
       int temp=igraph_vcount((igraph_t *) VECTOR(complist)[lci]);
       if (temp>lcomp_size) {
          lcomp_size= temp;
          lcomp_i=lci;
       }
     }
     igraph_copy(giant,(igraph_t *) VECTOR(complist)[lcomp_i]); 
     free_complist(&complist);
}


int compute_MD(const igraph_t *giant) {     
     int gN=(int )igraph_vcount(&giant);

     // Computing all pairs of shortest paths
     igraph_matrix_init(&D, 0, 0);
     igraph_shortest_paths(&giant, &D, igraph_vss_all(), igraph_vss_all(), IGRAPH_ALL);
    // printD(gN);
     
     // Init partition data structure
     part=new vector<int>*[gN];
     for(int i =0 ; i < gN ; i++ ) part[i] = new vector<int>;

   /* // int myints[] = {0,1};
     //vector<int> sensors (myints, myints + sizeof(myints) / sizeof(int) );
     vector<int> sensors;
     for(int i=0; i<cd+1; i++) sensors.push_back(i);
   //  for(int i=0; i<cd; i++) sensors.push_back(i+500);

     printf("%d: %d %lu\n",cd,testMD(gN,sensors), sensors.size());*/
     igraph_integer_t cl_num;
     igraph_integer_t num_mcl;
         
     // Approximating MD
     bool cont=true;
     int MD=0;
     while(cont) {
       pair<float, int> ret = H(gN);
      // printf("%f %d\n",ret.first, ret.second);
     //  fprintf(ki_graph,"s %d\n",ret.second);
       if (ret.first==0) cont=false;
       MD++;
     }
     
     for(int i =0 ; i < gN ; i++ ) delete part[i];
     delete[] part;
          
     return MD;
}
