#include "MD_compute.h"


pair<float, int> H(igraph_matrix_t *D,vector<int> **part, int N) {
    float best_res=100000000;
    int best_guess=-1;
    
    for(int v=0; v<N; v++) {
      //TODO Here we could check not to chose a v that's already added. 
      //Compute the size of the refined part with node v
      map<pair<int,vector<int> >,int> buckets;
	  for(int i=0; i<N; i++) {
	     pair<int,vector<int> > guess_pair = make_pair(MATRIX(*D,v,i), *part[i]);
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
		
	return make_pair(best_res,best_guess);
}


int H_Dyn(igraph_matrix_t *D, set<int>* candidates, int N, int source) {
    int min_max_count=100000000;
    int min_max_guess=-1;

    for(int v=0; v<N; v++) {
      //TODO Here we could check not to chose a v that's already added. 
      //Compute the size of the candidate set with node v
      int dist_counter[N];
      int max_dist=0;
      for(int i=0; i<N; i++) dist_counter[i]=0;

      for(set<int>::iterator it=candidates->begin(); it!=candidates->end(); ++it) {
        int temp=MATRIX(*D,v,*it);
        dist_counter[temp]++;
        if (temp>max_dist) max_dist=temp;
      }
      int max_count=0;
      for(int i=0; i<=max_dist; i++)
        if(dist_counter[i]>max_count) {
          max_count=dist_counter[i];
        }
      
      if (max_count<min_max_count) {
        min_max_count=max_count;
        min_max_guess=v;
      }	  
	}
	
	return min_max_guess;
}


int cut_approx(igraph_matrix_t* D, int** confusion, int N, int theta) {
    float best_count=-1;
    int best_guess=-1;
    
    for(int v=0; v<N; v++) {
      //TODO Here we could check not to chose a v that's already added. 
      //Compute the number of edges cut by node v
      int cut_count=0;
      
	  for(int i=0; i<N; i++) {
	    for(int j=i+1; j<N; j++) {
	      if (confusion[i][j]==1) 
	        if (abs(MATRIX(*D,v,i)-MATRIX(*D,v,j))>theta)
	          cut_count++;
	    }
	  }
	  
	  //Save best
	  if (cut_count>best_count) {
	    best_count=cut_count;
	    best_guess=v;
	  }
	  
	}
		
	return best_guess;
}

int compute_MD(igraph_matrix_t *D, int gN) {          
     // Init partition data structure
     vector<int> **part=new vector<int>*[gN];
     for(int i =0 ; i < gN ; i++ ) part[i] = new vector<int>;
         
     // Approximating MD
     bool cont=true;
     int MD=0;
     while(cont) {
       pair<float, int> ret = H(D,part,gN);
            
       // Add the best partition to the partition data structure
	   for(int i=0; i<gN; i++) {
	  		part[i]->push_back(MATRIX(*D,ret.second,i));
	   }

       if (ret.first==0) cont=false;
       MD++;
     }
     
     for(int i =0 ; i < gN ; i++ ) delete part[i];
     delete[] part;
          
     return MD;
}

int compute_DynMD(igraph_matrix_t *D,int gN) {
     set<int>* candidates=new set<int>;

     // Approximating DynMD
     int DynMD=0;

     //Find the source for each possible source
     for(int source =0; source<gN; source++) {
       bool cont=true;
       int steps_needed_for_source=0;
       
       // Init candidates->set 
       for (int i=0; i<gN; i++) candidates->insert(i);

       while(cont) {
         int best = H_Dyn(D,candidates, gN,source);
         
         // make measurement
		 int measurement = MATRIX(*D,best,source);
	
		 // Update the candidate set
	     for(int i=0; i<gN; i++) {
	    	if (MATRIX(*D,best,i)!=measurement)
	      		candidates->erase(i);
         }
		      
         //printf("%d %f %d\n",source, ret.first, ret.second);
         if (candidates->size()==1) cont=false;
           steps_needed_for_source++;
         }
         
         if(steps_needed_for_source>DynMD) DynMD=steps_needed_for_source;
         candidates->clear();
    }
    return DynMD;
}

int compute_RMD(igraph_matrix_t *D, int gN, int rho, int theta) {

     // Init cut_edge data structure
     int** confusion=new int*[gN];
     int count_conf=0;
     for(int i =0 ; i < gN ; i++ ) {
      confusion[i] = new int[gN];
      for(int j =i+1 ; j < gN ; j++ ) {
        if (MATRIX(*D,i,j)<rho+1) confusion[i][j]=0;
        else { 
         confusion[i][j]=1;
         count_conf++;
        }
      }
    }
          
     // Approximating MD
     bool cont=true;
     int MD=0;
     while(cont) {
       int best_guess = cut_approx(D, confusion, gN,theta);
     
       int still_left=0;
	   // Cut the edges
	   for(int i=0; i<gN; i++) {
	     for(int j=i+1; j<gN; j++) {
	        if (confusion[i][j]==1) {
	          if (abs(MATRIX(*D,best_guess,i)-MATRIX(*D,best_guess,j))>theta) {
	             confusion[i][j]=0;
	          }
	          else {
	             still_left++;
	          }
	        }
	     }
  	   }
     
       if (still_left==0) cont=false;
       MD++;
     }
     
     // Freeing memory
     for(int i =0 ; i < gN ; i++ ) delete[] confusion[i];
     delete[] confusion;

     return MD;
}


void dist_sets(igraph_matrix_t *D, int N, FILE* ki,char* network,int ecount,double p,int cd) {
	fprintf(ki,"\n"); 
    for(int v=0; v<N; v++) {
	  for(int i=v+1; i<N; i++) {
	     int dist_set=0;
	     //for(int j=0; j<N; j++)  {if(MATRIX(*D,v,j)!=MATRIX(*D,i,j)) dist_set++; else fprintf(ki,"%d %d %d\n",j,(int)MATRIX(*D,v,j),(int)MATRIX(*D,i,j)); }
	     for(int j=0; j<N; j++)  if(MATRIX(*D,v,j)!=MATRIX(*D,i,j)) dist_set++;

	     fprintf(ki,"%s %d %d %f %d ",network,N,ecount,p,cd);
	     fprintf(ki,"%d %d %d\n",v,i,dist_set);
	  }
	}
}
