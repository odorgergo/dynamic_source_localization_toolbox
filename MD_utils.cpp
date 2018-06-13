#include "MD_utils.h"

bool testMD(igraph_matrix_t *D, int N, vector<int> sensors) {
    for(int v=0; v<N; v++) {
      for(int w=v+1; w<N; w++) {
        bool foundS=false;
         for(vector<int>::iterator it = sensors.begin(); it != sensors.end(); ++it) {
           if (MATRIX(*D,v,*it)!=MATRIX(*D,w,*it)) {
             foundS=true;
             break;
           }
         }
       if (!foundS) {
         printf("%d %d %d %d %d\n",v,w,(int)MATRIX(*D,v,sensors[0]),(int)MATRIX(*D,w,sensors[0]),sensors[0]);
         return false;
       }
     }
    }
    return true;
}

bool testDMD(igraph_matrix_t *D, int N, vector<int> sensors) {
    for(int v=0; v<N; v++) {
      for(int w=v+1; w<N; w++) {
        bool foundS=false;
         for(vector<int>::iterator it = sensors.begin(); it != sensors.end(); ++it) {
           for(vector<int>::iterator it2 = sensors.begin(); it2 != sensors.end(); ++it2) {
             if ((MATRIX(*D,v,*it)-MATRIX(*D,w,*it)) !=  (MATRIX(*D,v,*it2)-MATRIX(*D,w,*it2))) {
               foundS=true;
               break;
             }
            }
         }
       if (!foundS) {
         printf("%d %d %d %d\n",v,w,(int)MATRIX(*D,v,sensors[0]),(int)MATRIX(*D,w,sensors[0]));
         return false;
       }
     }
    }
    return true;
}


void printD(igraph_matrix_t *D, int N) {
   for(int i=0; i<N;i++) {
     printf("i=%d: ", i);
     for(int j=0; j<N;j++) {
       printf("%d ", (int) MATRIX(*D,i,j));
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
