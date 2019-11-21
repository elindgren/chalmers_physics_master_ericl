#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>


double magnetization(int N, double lat[][N]){
    // Calculation of average magnetization of lattice goes here
    
}


void metropolis(int N, double T, double J, double kB, double lat[][N]){
    // Implementation of metropolis algorithm goes here
}

void simulate_lattice(int N, double T,  double J, double kB){
    /* RNG */
    srand(time(NULL)); // Set the seed for rand
    double rand_disp;  // random displacement

    /* Variables */
    double (*lat)[N] = malloc(sizeof(double [N][N]));
    double *m = malloc(N * sizeof(double));  // Average magnetizations for each iteration
    double tol = 1;  // Tolerance for convergence
    FILE *f;  // Output file
    char filename[50] = "";
    char dir[100] = "problem3_output/";
    char Temperature[50];
    char output[50];
    gcvt(T, 10, Temperature);
    strcat(dir, Temperature);

    /* Initialize lattice randomly */
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            rand_disp = round((double)rand() / (double)RAND_MAX);  // -1 or +1
            if(rand_disp == 0.0){
                rand_disp = -1.0;
            }
            lat[i][j] = rand_disp;
        }
    }
    
    // Run metropolis to tol

    // Save results
    filename[0] = '\0';  // Empty filename string
    strcat(filename, dir);
    strcat(filename, "_lattice.dat");
    f = fopen(filename, "w");
    for (int i = 0; i < N; i++)
    {
        for (int j=0; j<N; j++){
            fprintf(f, "%.1f \t", lat[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);


    // Free variables
    free(lat); lat=NULL;
}


int main(){
    // Define constants: J, kB
    double J = 1;
    double kB = 1; 

    // Define temperatures
    double T[] = {0, 10, 20};

    // Define sizes N
    int N[] = {2, 4, 6};

    /* Iterate over all temperatures and all sizes */
    simulate_lattice(2, 1, J, kB);




}