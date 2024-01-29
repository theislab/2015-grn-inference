
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <string>

#include "paramestFuns.h"

using namespace std;

int main(int argc, char** argv) {

	if(argc < 3)
	{
		std::cerr << "input <temperature scale> <burning>\n";
       	return EXIT_FAILURE;
    }    
  	const double temperatureN = atof(argv[1]); // temperature discretisation (e.g. 100)
    const int    burn         = atof(argv[2]);  	
  	const int    geneInd      = 2;
    const int    emu1Ind      = 1;
    const int    emu2Ind      = 3;
    const int    Ob1    	  = 100;         // number observations
  	
    // -------------------------------------------------------------------------
    
    double    	dt     = 1;
    int       	N      = 100;        // vector length

    double 		sigmaC = 10;         // noise standard deviation for likelihood
    double      sigma1[Ob1];
    
    const int   ngenes = 3;
    const int 	nPar   = 6;          // number parameters
    const int   nReps  = 3; 		 // n. Replicates
    
    const int   totalT = 6;
    int temperature[totalT] = {1,10,20,30,40,50,60,70,80,90,100};

    double                              params[nPar];
    std::vector<double>  				observed1(Ob1);
    
    // -------------------------------------------------------------------------
    
    std::ifstream inFile;
    int k = 0;

    double data1[Ob1][nReps];       // data (noisy observations)
    
    inFile.open("DataBranch_1_2.txt");
    for(int i = 0; i < Ob1; i++){
        for(int j = 0; j < nReps; j++){
            inFile >> data1[i][j];
        }
    }
    inFile.close();
    
    // @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @
    
    int times1[Ob1];
    
    inFile.open("TimesBranch_1_2.txt");
    for(int i = 0; i < Ob1; i++){
        inFile >> times1[i];
    }
    for(int i = 0; i < Ob1; i++){
        times1[i] = times1[i]+1;
    }
    inFile.close();
    
    // @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @
    
    std::vector< std::vector<double> >  emulTotA(N, std::vector<double>(ngenes));

	// In order to pass these arrays to a function, the second dimension must be defined
	// inside the function [see solveODE()]
    double emul1A[N][1];
	double emul2A[N][1];
    
    inFile.open("Interpolators1.txt");
    for(int i = 0; i < N; i++){
        for(int j = 0; j < ngenes; j++){
            inFile >> emulTotA[i][j];
        }
    }
    inFile.close();
    
    for(int i = 0; i < N; i++){
    	emul1A[i][0] = emulTotA[i][emu1Ind-1];
	}
    for(int i = 0; i < N; i++){
    	emul2A[i][0] = emulTotA[i][emu2Ind-1];
	}
	
    // -------------------------------------------------------------------------
    
	for(int i=0; i<Ob1; i++){
    	sigma1[i] = sigmaC;
    }
        
    // -------------------------------------------------------------------------      

    double  init1[nReps];
    
    for(int i=0; i<nReps; i++){
        init1[i] = emulTotA[0][geneInd-1];
    }
    
    // -------------------------------------------------------------------------
 	
 	double value;
 	double logEvidence = 0;
 	std::vector<double> likeSum(totalT);
 		
    std::string fileName = "Risulta_";
    
	int 	counter = 0;
	int 	denomin = 0; 
	double 	like    = 0;
	
		logEvidence = 0;
		counter     = 0;
		denomin     = 0;
		like        = 0;
	
		for(int t=0; t<totalT; t++){	
			
			fileName += to_string((int)(temperature[t]));
	    	inFile.open(fileName.c_str());	    	
	
			counter    = 0;
			denomin    = 0;
			likeSum[t] = 0;		
		
			while(!inFile.eof( ))
			{

				inFile >> value;
				counter++;
		
				if(counter > (nPar+1)*burn)
				{
        			for(int i=0; i<nPar; i++){        	
        				params[i] = value;
        				inFile >> value;
		        	}		        	

		        	observed1 = solveODE(params, nReps, init1, N, Ob1, emul1A, emul2A, dt, times1);
	    		
	    			like = 0;
	    			
					for(int i = 0; i < Ob1; i++){
			    	    for(int j = 0; j < nReps; j++){
    	    			    like = like + normPDF(data1[i][j],observed1[i],sigma1[i],1) + log(pow(2*PI*pow(sigma1[i],2),-0.5));
			        	}
				    }
				    like = like/nReps;
	    			    			
    				denomin++;
    				if(like>-double(RAND_MAX))
	    				likeSum[t] = likeSum[t] + like;	    				
    	    	
		        }		        	        
	                
			}
		
			likeSum[t] = likeSum[t]/denomin;
		
			std::cout << likeSum[t] << "\t";
    
	    	inFile.close( );
	    	fileName.erase(8);
	
		}
		
		// We consider the value at temperature 0.01 as if it is at temperature 0
		temperature[0] = 0;
		
		for(int t=1; t<totalT; t++){
			logEvidence = logEvidence + ((temperature[t]/temperatureN)-(temperature[t-1]/temperatureN))*0.5*(likeSum[t]+likeSum[t-1]);	
		}
	
		std::cout << logEvidence << "\n";


    return 0;

}


