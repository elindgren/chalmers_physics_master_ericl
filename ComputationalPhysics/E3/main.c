#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>
#define PI 3.141592653589793238


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

double weigth_function(double x, double y, double z){
    return pow(PI, -1.5) * exp(-( x*x + y*y + z*z));
}

double integrand3D(double x, double y, double z){
    return (x*x + x*x*y*y + x*x*y*y*z*z) * weigth_function(x,y,z);
}

void task3(){
    /* Perform Metropolis sampled Monte Carlo Integration of 3D integrand */
    /* Sampling parameters */
    int nbr_burn = 1000;  // Burn-in steps
    int nbr_prod = 1000;   // Production steps
    double delta = 2;
    double accepted_steps = 0;

    /* Define data structures */
    double *r = malloc(3 * sizeof(double));  // Current position
    double *r_proposal = malloc(3 * sizeof(double));  // Proposal position
    double *f_trace = malloc(nbr_prod * sizeof(double));  // Trace matrix containing all function evaluations

    /* RNG */
    double u;
	const gsl_rng_type *T;
	gsl_rng *q;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	q = gsl_rng_alloc(T);
	gsl_rng_set(q,time(NULL));
    
    /* Probabilities*/
    double P_current;
    double P_proposal;
    double accept;
    
    /* File handling */
    FILE *f;  // Output file

    /* Set starting position to 0, since weight function is symmetric */
	r[0] = 0; r[1] = 0; r[2] = 0;

    /* Perform Metropolis sampling */
    int steps;
    for(steps = 0; steps < nbr_burn + nbr_prod; steps++){
        /* Calculate current probability */
        P_current = weigth_function(r[0], r[1], r[2]);

        /* Suggest a new position */
        for(int i=0; i<3; i++){
            u = gsl_rng_uniform(q);
            r_proposal[i] = r[i] + delta*(u-0.5);
        }

        /* Calculate proposal position */
        P_proposal = weigth_function(r_proposal[0], r_proposal[1], r_proposal[2]);

        /* Calculate quotient */
        accept = P_proposal/P_current;

        /* Check if region of higher probability */
        /* If higher probability step there, or step there with probability accept */
        if(accept > 1 || accept > gsl_rng_uniform(q)){
            for(int i=0; i<3; i++){
                r[i] = r_proposal[i];
            }
            accepted_steps++;
        }

        /* Calculate integrand value for current position */
        if(steps >= nbr_burn){
            f_trace[steps-nbr_burn] = integrand3D(r[0], r[1], r[2]); 
        }
    }
    printf("Steps: %d \n", steps);
    printf("Accepted steps: %.2f \n", accepted_steps);
    printf("Acceptance ratio: %.2f \n", accepted_steps/(double)steps);

    /* Write trace to file */
    f = fopen("task3/task3.dat", "w");
    for (int i = 0; i < nbr_prod; i++)
    {
        fprintf(f, "%.16f \n", f_trace[i]);
    }
    fclose(f);

    /* Free variables */
    free(r); r = NULL; free(f_trace); f_trace = NULL;
}

void correlation_function(int N, double f[N], double phi[N]){
    /* Calculate sampled function average and average squared */
    double avg_f = 0;
    for(int i=0; i<N; i++){
        avg_f += f[i];
    }
    avg_f /= N; 
    printf("Average function value: %.2f \n", avg_f);
    // Remove the mean from the signal, and just have the deviations from the mean
    double *d = malloc(N * sizeof(double));
    for(int i=0; i<N; i++){
        d[i] = f[i] - avg_f;
    }
    // Now the average value of the deviations is 0, simplifying or correlation function

    double avg_d_sq = 0;
    for(int i=0; i<N; i++){
        avg_d_sq += d[i]*d[i];
    }
    avg_d_sq /= N;
    printf("Average deviation squard: %.2f \n", avg_d_sq);

    /* Calculate autocorrelation function */
    printf("Calculating correlation function \n");
    double d_kd_i = 0;
    for(int k=0; k<N; k++){
        if(k%100000 == 0){
            printf("k=%d \n", k);
        }
        d_kd_i = 0;
        for(int i=0; i<N-k; i++){
            d_kd_i += d[i+k]*d[i];
        }
        d_kd_i /= N;
        phi[k] = d_kd_i / avg_d_sq;
    }
}

double block_average(int N, int B, int M,double f[N], double F[M]){
    /* Calculates the block averaged value for s. M = N/B is the number of blocks. */

    /* Calculate F */
    double Fj;
    for (int j = 0; j<M; j++){
        Fj = 0;
        // Calculate the average value for block j
        for(int i = 0; i<B; i++){
            Fj += f[i+j*B];
        }
        Fj /= B;
        F[j] = Fj;
    }

    /* Calculate the variance of F and f */
    double mean_f = gsl_stats_mean(f, 1, N);  // Function, stride, length
    double var_f = gsl_stats_variance_m(f, 1, N, mean_f); // Function, stride, length, mean - Variance in f
    double mean_F = gsl_stats_mean(F, 1, M);  // Function, stride, length
    double var_F = gsl_stats_variance_m(f, 1, M, mean_F); // Function, stride, length, mean - Variance in F
    
    /* Estimate s */
    return B*var_F/var_f;


}

void task4(int calc_corr){
    /* Define variables */
    int nbr_of_lines = 1e6;
    double *phi = malloc(nbr_of_lines * sizeof(double));
    int B = 1000; // Block size
    int nbr_blocks = nbr_of_lines/B;
    double *F = malloc(nbr_blocks * sizeof(double));

    /* Load data */
    int i;
    FILE *f;
    FILE *in_file;
    double *data = malloc((nbr_of_lines) * sizeof (double));
    
    /* Read data from file. */
    in_file = fopen("MC.txt","r");
    for (i=0; i<nbr_of_lines; i++) {
        fscanf(in_file,"%lf",&data[i]);
    }
    fclose(in_file);

    /* Calculate block averaged approximation of s */
    double s = block_average(nbr_of_lines, B, nbr_blocks, data, F);

    printf("Block averaging: s=%.3f \n",s);

    if(calc_corr){
        /* Calculate correlation function */
        correlation_function(nbr_of_lines, data, phi);

        /* Write trace to file */
        f = fopen("task4/corr_func.dat", "w");
        for (int i = 0; i < nbr_of_lines; i++)
        {
            fprintf(f, "%.16f \n", phi[i]);
        }
        fclose(f);
    }
    

}


int main(){
    // Run tasks 
    // task1();
    // task2();
    // task3();
    task4(0);
}
