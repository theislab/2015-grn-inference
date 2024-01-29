
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
		std::cerr << "input <temperature>\n";
       	return EXIT_FAILURE;
    }
    const double temperature  = atof(argv[1]);
    const double temperatureN = 100;
    
    const int    geneIndA      = 3;  // index of target gene
    const int    emu1IndA      = 1;  // index of 1st input gene
    const int    emu2IndA      = 4;  // index of 2nd input gene
    
    const int    Ob1    	  = 106; // number observations of target gene (1st branch)
    const int    Ob2    	  = 118; // number observations of target gene (2nd branch)
    
    srand(time(NULL));

    double    	dt     = 1;
    int       	T      = 1e7;        // number MCMC iterations (for each for loop)
    int       	N      = 120;        // number of points of interpolators
    
    double 		sigmaC = 0.5;        // noise standard deviation for likelihood
    double      sigma1[Ob1];
    double      sigma2[Ob2];

    const double ngenes = 6;         // total number of GRN genes
    const int 	 nPar   = 6;         // number ODE parameters
    const int    nReps  = 3; 		 // number of replicates
    
    double tau[nPar];				 // standard deviation for Gaussian proposal
    
    double alpha  = 0;
    double u      = 0;
    
    // -------------------------------------------------------------------------
    
    std::ifstream inFile;
    int k = 0;

    double data1A[Ob1][nReps];       // data (noisy observations 1st branch)
    double data2A[Ob2][nReps];       // data (noisy observations 2nd branch)
    
    inFile.open("DataBranch_1_3.txt");
    for(int i = 0; i < Ob1; i++){
        for(int j = 0; j < nReps; j++){
            inFile >> data1A[i][j];
        }
    }
    inFile.close();
    
    inFile.open("DataBranch_2_3.txt");
    for(int i = 0; i < Ob2; i++){
        for(int j = 0; j < nReps; j++){
            inFile >> data2A[i][j];
        }
    }
    inFile.close();
        
    // @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @
    
    int times1[Ob1];				 // pseudo-times for noisy data (1st branch)
    int times2[Ob2];				 // pseudo-times for noisy data (2nd branch)
    
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
    
    double emul1Aa[N][1];			 // interpolator 1st input gene (1st branch)
	double emul2Aa[N][1];			 // interpolator 2nd input gene (1st branch)
    double emul1Ba[N][1];			 // interpolator 1st input gene (2st branch)
	double emul2Ba[N][1];			 // interpolator 2nd input gene (2st branch)
    
    inFile.open("Interpolators1.txt");
    for(int i = 0; i < N; i++){
        for(int j = 0; j < ngenes; j++){
            inFile >> emulTotA[i][j];
        }
    }
    inFile.close();
    
    for(int i = 0; i < N; i++){
    	emul1Aa[i][0] = emulTotA[i][emu1IndA-1];
	}
    for(int i = 0; i < N; i++){
    	emul2Aa[i][0] = emulTotA[i][emu2IndA-1];
	}
	
    inFile.open("Interpolators2.txt");
    for(int i = 0; i < N; i++){
        for(int j = 0; j < ngenes; j++){
            inFile >> emulTotB[i][j];
        }
    }
    inFile.close();
    
    for(int i = 0; i < N; i++){
    	emul1Ba[i][0] = emulTotB[i][emu1IndA-1];
	}
    for(int i = 0; i < N; i++){
    	emul2Ba[i][0] = emulTotB[i][emu2IndA-1];
	}
    
    // -------------------------------------------------------------------------
    
	for(int i=0; i<Ob1; i++){
    	sigma1[i] = sigmaC;
    }
    for(int i=0; i<Ob2; i++){
    	sigma2[i] = sigmaC;
    }
    
    // -------------------------------------------------------------------------  
    
    std::vector<double>  observed1a(Ob1);
    std::vector<double>  observed2a(Ob2);
    
    double params[nPar];
    double paramsNew[nPar];
    double init1a[nReps];
    double init2a[nReps];
            
    // initial conditions
    for(int i=0; i<nReps; i++){
        init1a[i] = emulTotA[0][geneIndA-1];
        init2a[i] = emulTotB[0][geneIndA-1];
    }
        
    std::vector<double>                accept(T, 0);                            // vector of acceptance
    std::vector<double>                Likel(T, 0);                             // vector of likelihood
    std::vector< std::vector<double> > state(T, std::vector<double>(nPar));	    // vector of states
    
    // -------------------------------------------------------------------------
    
    // uniform prior over kinetic (ODE) parameters: the following define max values
    double priorMax[nPar] = {log(2000), log(50), log(2000), log(30), log(2000), log(30)};
        
    double like    = 0;
    double likeNew = 0;
    
    // -------------------------------------------------------------------------
    
    double sumAccept         = 0;
    int    elem              = 0;    
    double averageSumAccept  = 0;
    
    std::vector< std::vector<double> > S(T, std::vector<double>(nPar));
	std::vector<double>                sLik(T,0);    
 	std::string                        fileName = "Risulta_";
 	
 	int StepTau  = 50000;
 	int negative = 0;
    
    tau[0] = {0.5};   // standard deviation for MCMC sampling (parameter 1)
	tau[1] = {0.01};  // standard deviation for MCMC sampling (parameter 2)
	tau[2] = {0.1};   // standard deviation for MCMC sampling (parameter 3)
	tau[3] = {0.01};  // standard deviation for MCMC sampling (parameter 4)
	tau[4] = {0.1};   // standard deviation for MCMC sampling (parameter 5)
	tau[5] = {0.01};  // standard deviation for MCMC sampling (parameter 6)
	
	// initial parameters values are sampled from a Gaussian with given mean and std
    params[0] = log(randn(800,50));
	params[1] = log(randn(1,0.1));
    params[2] = log(randn(150,50));
	params[3] = log(randn(1,0.1));
    params[4] = log(randn(150,50));
    params[5] = log(randn(1,0.1));
        
    double paramsFinal[nPar];
    
    for(int i=0; i<nPar; i++){
    	paramsFinal[i] = params[i];
    }
    
    
   for(int iter=1; iter<21; iter++){
   
	   for(int i=0; i<T; i++){
		   accept[i] = 0;
	   }	
   
	   for(int i=0; i<nPar; i++){
		   params[i] = paramsFinal[i];
	   }

	   observed1a = solveODE(params, nReps, init1a, N, Ob1, emul1Aa, emul2Aa, dt, times1);
	   observed2a = solveODE(params, nReps, init2a, N, Ob2, emul1Ba, emul2Ba, dt, times2);
	   
	   like = 0;
	   
	   for(int i = 0; i < Ob1; i++){
		   for(int j = 0; j < nReps; j++){
			   like = like + normPDF(data1A[i][j],observed1a[i],sigma1[i],temperature/temperatureN);
		   }
	   }
	   for(int i = 0; i < Ob2; i++){
		   for(int j = 0; j < nReps; j++){
			   like = like + normPDF(data2A[i][j],observed2a[i],sigma2[i],temperature/temperatureN);
		   }
	   }
   
	   for(int ii=1; ii<T; ii++){
		   
		   for (int i=0; i<nPar; i++) {
			   paramsNew[i] = randn(params[i],tau[i]);
		   }
	   
		   observed1a = solveODE(paramsNew, nReps, init1a, N, Ob1, emul1Aa, emul2Aa, dt, times1);
		   observed2a = solveODE(paramsNew, nReps, init2a, N, Ob2, emul1Ba, emul2Ba, dt, times2);
	   
		   likeNew = 0;
		   
		   // variable "negative" is a control over negative gene expression values: if ODE solution
		   // has more than a certain number of negative values, than the sample is not accepted
		   
		   negative = 0;
		   for(int i=0; i<Ob1; i++){
			   if(observed1a[i]<0){
				   negative++;
			   }
		   }			
		   for(int i=0; i<Ob2; i++){
			   if(observed2a[i]<0){
				   negative++;
			   }
		   }
	   
		   if((paramsNew[0]<priorMax[0])&(paramsNew[1]<priorMax[1])&(paramsNew[2]<priorMax[2])&(paramsNew[3]<priorMax[3])&(paramsNew[4]<priorMax[4])&(paramsNew[5]<priorMax[5])&(negative<6))
		   {
		   
			   for(int i=0; i<Ob1; i++){
				   for(int j=0; j<nReps; j++){
						   likeNew = likeNew + normPDF(data1A[i][j],observed1a[i],sigma1[i],temperature/temperatureN);
				   }
			   }
			   for(int i=0; i<Ob2; i++){
				   for(int j=0; j<nReps; j++){
						   likeNew = likeNew + normPDF(data2A[i][j],observed2a[i],sigma2[i],temperature/temperatureN);
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
		   
		   		
		   		// change std of different parameters (one per time), to maintain an acceptance rate around 0.23
			   if(sumAccept/StepTau > 0.23)
			   {
				   tau[elem] = tau[elem]*1.2;
				   elem++;
			   }
			   if(sumAccept/StepTau < 0.23)
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
   
	   std::cout << "[Chain " << iter << " of 20]  " << "Acceptance ratio (temperature " << temperature/temperatureN << "): " << averageSumAccept/T << "\n" ;
   
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
	   
	   if(jj<100000)
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
	   
	   if(jj>100000)
	   {
		   double V = round(jj/2);
		   double W = round(V/100000);
	   
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
	   fileName.erase(8);
	   
   }

	return 0;	

}
