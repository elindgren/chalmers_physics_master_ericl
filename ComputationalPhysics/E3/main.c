#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>
#define PI 3.14


double integrand1D(double x){
    return x*(1.0-x);
}

double integrand1DImportance(double x){
    return 2.0*x*(1.0-x)/(PI * sin(PI*x));
}

double inverseCDF1D(double y){
    return acos(1-2*y)/PI;
}


void monteCarlo(int N, double vals[], double errs[], double samples[], int idx, int importanceSampling){
    /* Monte-Carlo integrates the integrand from 0 to 1 */

    // Initialize vector of N points
    double *f = malloc(N * sizeof(double));
    double xi;
    double yi; 
    

    /* Initialize random number generator */
	double u;
	const gsl_rng_type *T;
	gsl_rng *q;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	q = gsl_rng_alloc(T);
	gsl_rng_set(q,time(NULL));

    // Populate with N random uniform numbers between 0 and 1
    for(int i=0; i<N; i++){
        // generate uniform random number
        u = gsl_rng_uniform(q);
        f[i] = u;
	    
    }
    
    // Evaluate integral for each such number and calculate mean value
    for(int i=0; i<N; i++){
        if(importanceSampling){
            yi = f[i];
            xi = inverseCDF1D(yi);
            f[i] = integrand1DImportance(xi);
            samples[i] = xi;
        }else{
            xi = f[i];
            f[i] = integrand1D(xi);
            samples[i] = xi;
        }
    }

    double mean = gsl_stats_mean(f, 1, N);  // Function, stride, length
    double var_f = gsl_stats_variance_m(f, 1, N, mean); // Function, stride, length, mean - Variance in function

    // We estimate the error according to the central limit theorem. Then the error of the MC-integration sigma_N = sigma_f/sqrt(N)
    errs[idx] = sqrt(var_f)/sqrt(N);
    vals[idx] = mean;


    // Free variables
    free(f); f=NULL;
    // free memory
	gsl_rng_free(q);
}

void task1() {
    /* Number of Ns to evaluate the integral for */
    int Ns = 4;
	int N[] = {10, 100, 1000, 10000};
    double *vals = malloc(Ns * sizeof(double));
    double *errs = malloc(Ns * sizeof(double));

    FILE *f;  // Output file

    // Iterate over all Ns, MC integrate function
    for(int i = 0; i<Ns; i++){
        double *samples = malloc(N[i] *sizeof(double));  // sampled points
        monteCarlo(N[i], vals, errs, samples, i, 0);
        // Write sampled points to file
        free(samples); samples=NULL;
    }

    //  Write to file 
    
    f = fopen("task1/task1.dat", "w");
    for (int i = 0; i < Ns; i++)
    {
        fprintf(f, "%d \t %.16f \t %.16f \n", N[i], vals[i], errs[i]);
    }
    fclose(f);

    // Free variables
    free(vals); vals=NULL; free(errs); errs=NULL;
}
 
void task2() {
    /* Number of Ns to evaluate the integral for */
    int Ns = 4;
	int N[] = {10, 100, 1000, 10000};
    double *vals = malloc(Ns * sizeof(double));
    double *errs = malloc(Ns * sizeof(double));

    FILE *f;  // Output file
    

    // Iterate over all Ns, MC integrate function
    for(int i = 0; i<Ns; i++){
        double *samples = malloc(N[i] *sizeof(double));  // sampled points
        monteCarlo(N[i], vals, errs, samples, i, 1);
        
        // Write sampled points to file
        char sample_filename[100] = "task2/N=";
        char size[50];
        gcvt(N[i], 5, size);
        strcat(sample_filename, size);
        strcat(sample_filename, "_samples.dat");
        f = fopen(sample_filename, "w");
        for (int j = 0; j < N[i]; j++)
        {
            fprintf(f, "%.16f \n", samples[j]);
        }
        fclose(f);
    }

    //  Write to file 
    
    f = fopen("task2/task2.dat", "w");
    for (int i = 0; i < Ns; i++)
    {
        fprintf(f, "%d \t %.16f \t %.16f \n", N[i], vals[i], errs[i]);
    }
    fclose(f);

    // Free variables
    free(vals); vals=NULL; free(errs); errs=NULL;
}

int main(){
    // Run tasks 
    task1();
    task2();
}
