
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

	if(argc < 2)
	{
		std::cerr << "input <burning>\n";
       	return EXIT_FAILURE;
    }    
  	const double temperatureN = 100;
    const int    burn         = atof(argv[1]);  	
  	const int    geneInd      = 3;
    const int    emu1Ind      = 1;
    const int    emu2Ind      = 4;
    const int    Ob1    	  = 106;         // number observations
    const int    Ob2    	  = 118;         // number observations
  	
    // -------------------------------------------------------------------------
    
    double    	dt     = 1;
    int       	N      = 120;        // vector length

    double 		sigmaC = 20.0;        // noise standard deviation for likelihood
    double      sigma1[Ob1];
    double      sigma2[Ob2];
    
    const int   ngenes = 6;
    const int 	nPar   = 6;          // number parameters
    const int   nReps  = 3; 		 // n. Replicates
    
    const int   totalT = 11;
    int temperature[totalT] = {1,10,20,30,40,50,60,70,80,90,100};

    double                              params[nPar];
    std::vector<double>  				observed1(Ob1);
    std::vector<double>  				observed2(Ob2);
    
    // -------------------------------------------------------------------------
    
    std::ifstream inFile;
    int k = 0;

    double data1[Ob1][nReps];       // data (noisy observations)
    double data2[Ob2][nReps];       // data (noisy observations)
    
    inFile.open("DataBranch_1_3.txt");
    for(int i = 0; i < Ob1; i++){
        for(int j = 0; j < nReps; j++){
            inFile >> data1[i][j];
        }
    }
    inFile.close();
    
    inFile.open("DataBranch_2_3.txt");
    for(int i = 0; i < Ob2; i++){
        for(int j = 0; j < nReps; j++){
            inFile >> data2[i][j];
        }
    }
    inFile.close();
    
    // @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @
    
    int times1[Ob1];
    int times2[Ob2];
    
    inFile.open("TimesBranch_1_3.txt");
    for(int i = 0; i < Ob1; i++){
        inFile >> times1[i];
    }
    for(int i = 0; i < Ob1; i++){
        times1[i] = times1[i]+1;
    }
    inFile.close();
    
    inFile.open("TimesBranch_2_3.txt");
    for(int i = 0; i < Ob2; i++){
        inFile >> times2[i];
    }
    for(int i = 0; i < Ob2; i++){
        times2[i] = times2[i]+1;
    }
    inFile.close();
    
    // @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @
    
    std::vector< std::vector<double> >  emulTotA(N, std::vector<double>(ngenes));
    std::vector< std::vector<double> >  emulTotB(N, std::vector<double>(ngenes));

	// In order to pass these arrays to a function, the second dimension must be defined
	// inside the function [see solveODE()]
    double emul1A[N][1];
	double emul2A[N][1];
    double emul1B[N][1];
	double emul2B[N][1];
    
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
	
    inFile.open("Interpolators2.txt");
    for(int i = 0; i < N; i++){
        for(int j = 0; j < ngenes; j++){
            inFile >> emulTotB[i][j];
        }
    }
    inFile.close();
    
    for(int i = 0; i < N; i++){
    	emul1B[i][0] = emulTotB[i][emu1Ind-1];
	}
    for(int i = 0; i < N; i++){
    	emul2B[i][0] = emulTotB[i][emu2Ind-1];
	}    
    
    // -------------------------------------------------------------------------
    
	for(int i=0; i<Ob1; i++){
    	sigma1[i] = sigmaC;
    }
    for(int i=0; i<Ob2; i++){
    	sigma2[i] = sigmaC;
    }
    
    // -------------------------------------------------------------------------      

    double  init1[nReps];
    double  init2[nReps];
    
    for(int i=0; i<nReps; i++){
        init1[i] = emulTotA[0][geneInd-1];
        init2[i] = emulTotB[0][geneInd-1];
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
			    	observed2 = solveODE(params, nReps, init2, N, Ob2, emul1B, emul2B, dt, times2);
	    		
	    			like = 0;
	    			
					for(int i = 0; i < Ob1; i++){
			    	    for(int j = 0; j < nReps; j++){
    	    			    like = like + normPDF(data1[i][j],observed1[i],sigma1[i],1) + log(pow(2*PI*pow(sigma1[i],2),-0.5));
			        	}
				    }
					for(int i = 0; i < Ob2; i++){
    	    			for(int j = 0; j < nReps; j++){
			    	        like = like + normPDF(data2[i][j],observed2[i],sigma2[i],1) + log(pow(2*PI*pow(sigma2[i],2),-0.5));
        				}
				    }
				    like = like/(2*nReps);
	    			    			
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


