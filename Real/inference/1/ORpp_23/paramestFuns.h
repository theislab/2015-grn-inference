
#ifndef TOGGLEFUN_H
#define	TOGGLEFUN_H

#define PI 3.14159265358979323846

double min(double a, double b) {
    return (a<b)?a:b;
}

double max(double a, double b) {
    return (a>b)?a:b;
}

double randn(double mu, double sigma) {
        static bool deviateAvailable=false;        // flag
        static float storedDeviate;                // deviate from previous calculation
        double dist, angle;
       
        // If no deviate has been stored, the standard Box-Muller transformation is performed, producing two independent
	// normally-distributed random deviates. One is stored for the next round, and one is returned.
        if (!deviateAvailable) {
               
                // choose a pair of uniformly distributed deviates, one for the distance and one for the angle, and perform transformations
                dist=sqrt( -2.0 * log(double(rand()) / double(RAND_MAX)) );
                angle=2.0 * PI * (double(rand()) / double(RAND_MAX));
               
                // calculate and store first deviate and set flag
                storedDeviate=dist*cos(angle);
                deviateAvailable=true;
               
                // calculate return second deviate
                return dist * sin(angle) * sigma + mu;
        }
      
        // If a deviate is available from a previous call to this function, it is returned, and the flag is set to false.
        else {
                deviateAvailable=false;
                return storedDeviate*sigma + mu;
        }
}

std::vector<double> solveODE(double para[], int nReps, double init[], int N, int Ob, double emul1[][1], double emul2[][1], double dt, int times[])
{
	
	std::vector<double> x(N+1);
	std::vector<double> xobs(Ob);

	x[0] = init[0];

	double k1 = 0;
	double k2 = 0;
	double k3 = 0;
	double k4 = 0;

	for(int t=1; t<N; t++){
            
            k1 = exp(para[0])*((1/(1+pow(exp(para[2])/emul1[t-1][0],exp(para[3]))))+(1/(1+pow(exp(para[4])/emul2[t-1][0],exp(para[5]))))) - exp(para[1])*x[t-1];
            k2 = exp(para[0])*((1/(1+pow(exp(para[2])/(0.5*(emul1[t-1][0]+emul1[t][0])),exp(para[3]))))+(1/(1+pow(exp(para[4])/(0.5*(emul2[t-1][0]+emul2[t][0])),exp(para[5]))))) - exp(para[1])*(x[t-1]+0.5*k1);
            k3 = exp(para[0])*((1/(1+pow(exp(para[2])/(0.5*(emul1[t-1][0]+emul1[t][0])),exp(para[3]))))+(1/(1+pow(exp(para[4])/(0.5*(emul2[t-1][0]+emul2[t][0])),exp(para[5]))))) - exp(para[1])*(x[t-1]+0.5*k2);
            k4 = exp(para[0])*((1/(1+pow(exp(para[2])/emul1[t][0],exp(para[3]))))+(1/(1+pow(exp(para[4])/emul2[t][0],exp(para[5]))))) - exp(para[1])*(x[t-1]+k3);
            
            x[t] = x[t-1] + (k1+2*k2+2*k3+k4)/6;
		                
	}
	x[N] = x[N-1] + dt*(exp(para[0])*((1/(1+pow(exp(para[2])/emul1[N-1][0],exp(para[3]))))+(1/(1+pow(exp(para[4])/emul2[N-1][0],exp(para[5]))))) - exp(para[1])*x[N-1]);
		
	for (int i=0; i<Ob; i++) {
		xobs[i] = x[times[i]-1];
	}

	return xobs;

}

double normPDF(double y, double x, double sigma, double temperature)
{
	double like;
	like = -0.5 * pow(((y - x)/sigma),2) * temperature;
	return like;
}

#endif	/* TOGGLEFUN_H */