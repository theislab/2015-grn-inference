
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
    
	if(argc < 4)
	{
		std::cerr << "input <temperature> <number of temperatures> <gene index> <emulator-1> <Ob>\n";
       	return EXIT_FAILURE;
    }
    const double temperature  = atof(argv[1]);
    const double temperatureN = atof(argv[2]);
    const int    geneInd      = atof(argv[3]);
    const int    emu1Ind      = atof(argv[4]);
    const int    Ob1    	  = atof(argv[5]);         // number observations
    
    srand(time(NULL));  // NB: doesn't work for multiple time a second        

    double    	dt     = 1;
    int       	T      = 1e7;        // number iterations    
    int       	N      = 120;        // vector length
    
    double 		sigmaC = 1;          // noise standard deviation for likelihood
    double      sigma1[Ob1];

    const double ngenes = 3;
    const int 	 nPar   = 4;          // number parameters
    const int    nReps  = 1; 		 // n. Replicates
    
    // Standard deviation for Gaussian proposal
    double tau[nPar];
    
    double alpha  = 0;		// temporary variable
    double u      = 0;		// temporary variable
    
    // -------------------------------------------------------------------------
    
    std::ifstream inFile;
    int k = 0;

    double data1[Ob1][nReps];       // data (noisy observations)
    
    inFile.open("/home/icb/andrea.ocone/Real/inference/DataBranch_3.txt");
    for(int i = 0; i < Ob1; i++){
        for(int j = 0; j < nReps; j++){
            inFile >> data1[i][j];
        }
    }
    inFile.close();
        
    // @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @
    
    int times1[Ob1];
    
    inFile.open("/home/icb/andrea.ocone/Real/inference/TimesBranch_3.txt");
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
    
    inFile.open("/home/icb/andrea.ocone/Real/inference/Interpolators.txt");
    for(int i = 0; i < N; i++){
        for(int j = 0; j < ngenes; j++){
            inFile >> emulTotA[i][j];
        }
    }
    inFile.close();
    
    for(int i = 0; i < N; i++){
    	emul1A[i][0] = emulTotA[i][emu1Ind-1];
	}
	    
    // -------------------------------------------------------------------------
    
	for(int i=0; i<Ob1; i++){
    	sigma1[i] = sigmaC;
    }
    
    // -------------------------------------------------------------------------  
    
    std::vector<double>  observed1(Ob1);
    
    double                              params[nPar];
    double                              paramsNew[nPar];
    double                              init1[nReps];
            
    // initial conditions
    for(int i=0; i<nReps; i++){
        init1[i] = emulTotA[0][geneInd-1];
    }
        
    std::vector<double>                accept(T, 0);                            // vector of acceptance
    std::vector<double>                Likel(T, 0);                             // vector of likelihood
    std::vector< std::vector<double> > state(T, std::vector<double>(nPar));	    // vector of states
    
    // -------------------------------------------------------------------------
    
    // uniform prior
    double priorMax[nPar] = {log(2000), log(50), log(2000), log(30)};
        
    double like    = 0;
    double likeNew = 0;
    
    // -------------------------------------------------------------------------
    
    double sumAccept         = 0;
    int    elem              = 0;    
    double averageSumAccept  = 0;
    
    std::vector< std::vector<double> > S(T, std::vector<double>(nPar));
	std::vector<double>                sLik(T,0);    
 	std::string                        fileName = "/home/icb/andrea.ocone/Real/inference/3/1p/Risulta_";
 	
 	int StepTau = 50000;
 	int negative = 0;
    
    tau[0] = {0.1};
	tau[1] = {0.01};
	tau[2] = {0.1};
	tau[3] = {0.01};
	
	// initial parameters sampled from a Gaussian with given mean and std
    params[0] = log(randn(150,50));
	params[1] = log(randn(1,0.1));
    params[2] = log(randn(150,50));
	params[3] = log(randn(1,0.1));
    
    double paramsFinal[nPar];
    
    for(int i=0; i<nPar; i++){
    	paramsFinal[i] = params[i];
    }
    
	for(int iter=1; iter<31; iter++){
	
		for(int i=0; i<T; i++){
			accept[i] = 0;
		}
	
	    for(int i=0; i<nPar; i++){
    		params[i] = paramsFinal[i];
    	}

	    observed1 = solveODE(params, nReps, init1, N, Ob1, emul1A, dt, times1);
    	
    	like = 0;
    	
		for(int i = 0; i < Ob1; i++){
    	    for(int j = 0; j < nReps; j++){
    	        like = like + normPDF(data1[i][j],observed1[i],sigma1[i],temperature/temperatureN);
        	}
	    }
	
	    for(int ii=1; ii<T; ii++){
	        
	        for (int i=0; i<nPar; i++) {
	            paramsNew[i] = randn(params[i],tau[i]);
    	    }
        
	        observed1 = solveODE(paramsNew, nReps, init1, N, Ob1, emul1A, dt, times1);
        
        	likeNew = 0;	// condition of parameters [0 < parameter < limPrior]
        	
        	negative = 0;
        	for(int i=0; i<Ob1; i++){
        		if(observed1[i]<0){
        			negative++;
        		}
			}
        
	        if((paramsNew[0]<priorMax[0])&(paramsNew[1]<priorMax[1])&(paramsNew[2]<priorMax[2])&(paramsNew[3]<priorMax[3])&(negative<6))
    	    {
            
        	    for(int i=0; i<Ob1; i++){
            	    for(int j=0; j<nReps; j++){
                	        likeNew = likeNew + normPDF(data1[i][j],observed1[i],sigma1[i],temperature/temperatureN);
	                }
    	        }
            
        	    alpha = min(0,likeNew-like);
            
	            u = ((double) rand() / (RAND_MAX));
            
    	        if(log(u)<alpha){
                
        	        for(int i=0; i<nPar; i++){
	                    params[i] = paramsNew[i];
    	            }
                
        	        like       = likeNew;                
            	    accept[ii] = 1;
            
	            }
        
    	    }
        
	        for(int i=0; i<nPar; i++){
    	        state[ii][i] = params[i];
        	}
        
	        Likel[ii] = likeNew;
        
	        if((ii%StepTau)==0){
        
		        sumAccept = 0;
		        if(elem>nPar-1)
	    	    	elem = 0;
        	
	    		for(int m=ii-StepTau; m<ii; m++){
    	    		sumAccept = sumAccept + accept[m];
			    }
		    
	        	std::cout << "Iteration: " << ii/StepTau + 1 << " of " << T/StepTau << "(acceptance rate = " << sumAccept*100/StepTau << "\%)\n";
        	
    	    	if(sumAccept/StepTau > 0.10)
        		{
					tau[elem] = tau[elem]*1.2;
					elem++;
	        	}
    	    	if(sumAccept/StepTau < 0.10)
        		{
					tau[elem] = tau[elem]/1.2;
					elem++;
	        	}
        	
    	    	averageSumAccept = averageSumAccept + sumAccept;
        	
	        }      	
    
    	}
    	
    	for(int i=0; i<nPar; i++){
    		paramsFinal[i] = params[i];
    	}
    
	    std::cout << "[Chain " << iter << " of 30]  " << "Acceptance ratio (temperature " << temperature/temperatureN << "): " << averageSumAccept/T << "\n" ;
    
    	//----------------------- Select accepted samples --------------------------
        
	    int jj = 0;
    
	    for(int i=0; i<T; i++){
        
    	    if(accept[i]>0){
            
        	    for(int k=0; k<nPar; k++){
            	    S[jj][k] = state[i][k];
	            }
            
    	        sLik[jj] = Likel[i];
            
        	    jj = jj+1;
	        }

    	}
    
		//------------------------------- Write files ------------------------------

    	fileName += to_string((int)(temperature));
    	fileName += "_";
    	fileName += to_string((int)(iter));

	    std::ofstream inFileOut;
    	inFileOut.open(fileName.c_str());

    	//-------------------------------------------------------------------------
    	
    	if(jj<50000)
		{
			double V = round(jj/5);
		
	    	for(int j=(int)V; j<jj-1; j++){
    	    	for(int i=0; i<nPar; i++){
        	    	inFileOut << S[j][i] << " ";
		        }
    		    inFileOut << sLik[j] << "\n";
	    	}
	    	for(int j=jj-1; j<jj; j++){
		        for(int i=0; i<nPar; i++){
    		        inFileOut << S[j][i] << " ";
        		}
		        inFileOut << sLik[j];
    		}
    	}
    	
    	if(jj>50000)
		{
			double V = round(jj/2);
			double W = round(V/50000);
		
	    	for(int j=(int)V; j<jj-1; j=j+(int)W){
    	    	for(int i=0; i<nPar; i++){
        	    	inFileOut << S[j][i] << " ";
		        }
    		    inFileOut << sLik[j] << "\n";
	    	}
	    	for(int j=jj-1; j<jj; j++){
		        for(int i=0; i<nPar; i++){
    		        inFileOut << S[j][i] << " ";
        		}
		        inFileOut << sLik[j];
    		}
    	}
    	
    	//-------------------------------------------------------------------------

	    inFileOut.close();
    	fileName.erase(51);
    	
	}

	return 0;	

}
