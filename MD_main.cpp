#include "MD_main.h"

int main(int argc, char * const argv[]){
	igraph_rng_seed(igraph_rng_default(), (int) time(NULL));
 /***********************************************************************
            Read arguments and initialize random generator	
 ************************************************************************/ 
	
	int a, iter_max, cd, rho_max, theta_max;
	double rad;
	char network [200],property [200],header[1200],outF [1200];
	vector<int> Ns;
	vector<float> ps;
	
	read_arguments(network, property, outF,header, &Ns, &ps, &a, &iter_max, &cd, &rad, &rho_max, &theta_max, argc, argv);
	
	for(int iter=0; iter<iter_max; iter++) {
    	printf("Iter: %d\n",iter);
		FILE* ki;
    	if(!strcmp(outF,"stdin")) {
    		ki=stdout;
    	}
    	else {
			char outF_final[1200];
			sprintf(outF_final ,"%s_iter%d.txt",outF,iter);
    		ki = fopen(outF_final,"w");
    	}
    	//Header
    	fprintf(ki,"%s\n",header);

        for(vector<int>::iterator it = Ns.begin(); it != Ns.end(); ++it) {
          int N = *it;
	      for(vector<float>::iterator it2 = ps.begin(); it2 != ps.end(); ++it2) {
	    	double p = *it2;
	    	if (p<0) p=pow(N,p); //p is given in the form N^(i/i+1)
	    	igraph_t graph;
	    	igraph_t giant;
	    	
	    	// Generating graph
	    	if(!strcmp(network,"Gnp")) igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNP, N, p, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
	    	else if (!strcmp(network,"Grid")) {
	    	     igraph_vector_t dimvector;      
	    	     igraph_vector_init(&dimvector, 2);
	    	     VECTOR(dimvector)[0]=a;
	    	     VECTOR(dimvector)[1]=a;
	    	     igraph_lattice(&graph,&dimvector,cd, false, false, false);
	    	     igraph_rewire_edges(&graph, p, false, false);
	    	}
	    	else if (!strcmp(network,"RGG")) {
	    		igraph_grg_game(&graph,N,rad, 1,0,0);
	    		igraph_rewire_edges(&graph, p, false, false);
	    	}
	    	
	    	// Taking the largest connected component
	    	largest_connected(&graph,&giant);
	    	
	    	//Output network properties
	    	int gN = (int )igraph_vcount(&giant);
	    	fprintf(ki,"%s %d %d %f ",network,gN,(int )igraph_ecount(&giant),p);
	    	if(!strcmp(network,"Grid")) fprintf(ki,"%d ",cd);
	    	if(!strcmp(network,"RGG")) fprintf(ki,"%f ",rad);	    	
	    	
	    	// Computing all pairs of shortest paths
	    	igraph_matrix_t D;
     		igraph_matrix_init(&D, 0, 0);
     		igraph_shortest_paths(&giant, &D, igraph_vss_all(), igraph_vss_all(), IGRAPH_ALL);
	    	
	    	//Computing property of interest
	    	if (!strcmp(property,"MD")) fprintf(ki,"%d\n",compute_MD(&D,gN));
	    	else if (!strcmp(property,"DynMD")) fprintf(ki,"%d\n",compute_DynMD(&D,gN));
	    	else if (!strcmp(property,"RMD")) 
	    		for(int rho=0; rho<rho_max; rho++) 
	    			for(int theta=0; theta<rho+1; theta++) 
	    				fprintf(ki,"%d\n",compute_RMD(&D,gN,rho,theta));
	    	else if (!strcmp(property,"MD+DynMD")) fprintf(ki,"%d %d\n",compute_MD(&D,gN),compute_DynMD(&D,gN));
	    	else if (!strcmp(property,"dist_sets")) dist_sets(&D,gN,ki,network,(int )igraph_ecount(&giant),p,cd); 
	    	
	    	// Freeing memory
	    	igraph_matrix_destroy(&D);   
            igraph_destroy(&giant);
            igraph_destroy(&graph);
          }
	    }
	    
	    if(strcmp(outF,"stdin")) fclose(ki);
	}
	
	
	return 1;
 }


/// a function that read the arguments of the program
int read_arguments(char network[], char property[], char outF[],char header[], vector<int>* Ns, vector<float>* ps, int*a, int* iter_max, int* cd, double* rad, int* rho_max, int* theta_max,int argc,char * const  argv[]){
	char p_input[200];
	char N_input[200];
	int N=0;
  /******  Default values *******/
	
    *a 			= 0;
    *cd			= 1;
    *rho_max    = 1;
    *theta_max  = 1;
    *rad		= 0;
    *iter_max	= 1;
    sprintf(network ,"none");
    sprintf(property,"none");
    sprintf(outF,"stdin");
    sprintf(N_input,"none");
    sprintf(p_input,"none");
    sprintf(header,"network N E p ");

    	
  /******  Reading arguments *******/

	int i;
	for (i=1; i<argc; ++i){
                    
		if      (!strcmp(argv[i],"-N")) N = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-cd"    )) *cd     = atoi(argv[++i]);
		else if (!strcmp(argv[i],"-rho"   )) *rho_max    = atoi(argv[++i]);
		else if (!strcmp(argv[i],"-the"   )) *theta_max    = atoi(argv[++i]);
		else if (!strcmp(argv[i],"-i"   )) *iter_max    = atoi(argv[++i]);
		else if (!strcmp(argv[i],"-rad"   )) *rad    = atof(argv[++i]);
        else if (!strcmp(argv[i],"-net"   )) sprintf(network   ,"%s",argv[++i]);
        else if (!strcmp(argv[i],"-prop"   )) sprintf(property   ,"%s",argv[++i]);
        else if (!strcmp(argv[i],"-o"    )) sprintf(outF    ,"%s",argv[++i]);
        else if (!strcmp(argv[i],"-p"    )) sprintf(p_input    ,"%s",argv[++i]);
        else if (!strcmp(argv[i],"-Ns"    )) sprintf(N_input    ,"%s",argv[++i]);
		
		else {
		  fprintf(stderr,"\nError:    '%s' not recognized as an argument type.\n",argv[i]);
		  fprintf(stderr,"          Before any argument you should first write the argument -key and then the value\n"); 
          fprintf(stderr,"        For more information read the README file \n\n");
		  exit(-1);
		}
       
	}
    
    if(!strcmp(N_input,"none")){        
		if(N==0) {
	    	fprintf(stderr,"N not specified\n\n");
    		exit(-1);
    	}
    	else Ns->push_back(N);
    }
    else if (!strcmp(N_input,"pow2")) for(int ni=1; ni<11; ni++) Ns->push_back(pow(2,ni)); 
    else if (!strcmp(N_input,"pow2_large")) for(int ni=11; ni<16; ni++) Ns->push_back(pow(2,ni)); 
    else if (!strcmp(N_input,"range10")) for(int ni=1; ni<11; ni++) Ns->push_back(ni*100); 
    else {
       	fprintf(stderr,"given N not implemented\n\n");
    	exit(-1);
    }

    
    if(!strcmp(p_input,"none")){        
   	    fprintf(stderr,"p not specified\n\n");
    	exit(-1);
    }
    else if (!strcmp(p_input,"0")) ps->push_back(0);
    else if (!strcmp(p_input,"1")) ps->push_back(1);
    else if (!strcmp(p_input,"1/2")) ps->push_back(0.5);
    else if (!strcmp(p_input,"range20")) for(int pi=0; pi<21; pi++) ps->push_back(pi*0.05);  
    else if (!strcmp(p_input,"range30")) {
      for(int pi=0; pi<10; pi++) ps->push_back(pi*0.005);
      for(int pi=1; pi<21; pi++) ps->push_back(pi*0.05);  
    }
    else if (!strcmp(p_input,"N^-1/4")) ps->push_back(-float(1)/4);
    else if (!strcmp(p_input,"N^-1/2")) ps->push_back(-float(1)/2);
    else if (!strcmp(p_input,"N^-2/3")) ps->push_back(-float(2)/3);
    else if (!strcmp(p_input,"N^-3/4")) ps->push_back(-float(3)/4);
    else if (!strcmp(p_input,"N^-4/5")) ps->push_back(-float(4)/5);
    else if (!strcmp(p_input,"N^-5/6")) ps->push_back(-float(5)/6);
    else {
       	fprintf(stderr,"given p not implemented\n\n");
    	exit(-1);
    }
    
    
    printf("ARGUMENTS:\n");
    printf("   Network : %s\n",network);
    if(!strcmp(network,"none")){        
        fprintf(stderr,"You did not specify network model\n\n");
        fprintf(stderr,"        For more information read the README file \n\n");
        exit(-1);
    }
    else if(!strcmp(network,"Gnp")) {
    	printf("   N= ");
    	for(vector<int>::iterator it = Ns->begin(); it != Ns->end(); ++it) printf("%d ",*it);
    	printf("\n   p=");
    	for(vector<float>::iterator it = ps->begin(); it != ps->end(); ++it) printf("%f ",*it);
	   	printf(" (%s)\n",p_input);
    }
	else if(!strcmp(network,"Grid")) {
	    *a=N;
	    N=N*N;
	    printf("   a= %d, N=a^2= %d\n",*a,N);
    	printf("   connection distance cd=%d\n",*cd);
    	printf("   rewiring probability p=");
    	for(vector<float>::iterator it = ps->begin(); it != ps->end(); ++it) printf("%f ",*it);
	   	printf(" (%s)\n",p_input);
	   	
	   	sprintf(header,"%scd, ",header);
    }
	else if(!strcmp(network,"RGG")) {
    	printf("   N= %d\n",N);
		if(*rad==0) {
	    	fprintf(stderr,"RGG radius not specified\n\n");
        	exit(-1);
        }
        else  printf("   RGG radius r=%f\n",*rad);
        printf("   rewiring probability p=");
    	for(vector<float>::iterator it = ps->begin(); it != ps->end(); ++it) printf("%f ",*it);
	   	printf(" (%s)\n",p_input);
	   	
	   	sprintf(header,"%srad, ",header);
    }
    else {
        fprintf(stderr,"given network model not implemented\n\n");
    	exit(-1);
    }
    
    printf("   Property to compute : %s\n",property);
    if(!strcmp(property,"none")){        
        fprintf(stderr,"You did not specify a property to compute\n\n");
        fprintf(stderr,"        For more information read the README file \n\n");
        exit(-1);
    }
    else if(!strcmp(property,"MD")) {
        sprintf(header,"%sMD ",header);
    }
    else if(!strcmp(property,"DynMD")) {
        sprintf(header,"%sDynMD ",header);
    }
    else if(!strcmp(property,"MD+DynMD")) {
        sprintf(header,"%sMD DynMD ",header);
    }
    else if(!strcmp(property,"dist_sets")) {
        sprintf(header,"%si j dist_set",header);
    }
    else if(!strcmp(property,"RMD")) {
    	printf("   rho_max= %d\n",*rho_max);
    	printf("   theta_max= %d\n",*theta_max);
        sprintf(header,"%srho theta RMD ",header);
    }
    else {
        fprintf(stderr,"given property computation not implemented\n\n");
    	exit(-1);
    }
    
    printf("   Iterations : %d\n",*iter_max);
    printf("   Output: %s\n",outF);


	return 1;
 }

