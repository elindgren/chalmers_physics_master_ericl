# include "c2_functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>


int main(){
    int len = 0;
    printf("Enter array length: ");
    scanf("%d", &len);  // Unsigned int
    printf("Scanned length: %d. \n", len);
    // Declare arrays
    double * A; 
    double * B;
    // Allocate memory
    A = malloc(len*sizeof(double));
    B = malloc(len*sizeof(double));
    // Pre-allocate arrays
    for (int i=0; i<len; i=i+1){
        A[i] = 1.0;
        B[i] = 2.0;
    }
    // Calculate scalar product
    problem1_scalar_product(A,B, len);
    // Do the same using VLAs
    double A_vla[len];
    double B_vla[len];
    for (int i=0; i<len; i=i+1){
        A_vla[i] = 1.0;
        B_vla[i] = 2.0;
    }
    problem1_scalar_product(A_vla,B_vla, len);

    // Calculate pointwise distances
    double points[len][3];  // Switch directions, transpose matrix
    // Declare positions
    for (int i=0; i<len; i=i+1){
        for (int j=0; j<3; j=j+1){
            double val = (double)i;
            points[i][j] = val;
            //printf("Value: %.2f \n", points[i][j]);
        }
    }
    // Calculate distance between points 0 and 1
    calc_distance(0,2, points);
    // Use gsl to create large array of uniformly drawn random numbers
    gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
    int size = 100000;
    // Allocate vector
    double *rand= malloc(size*sizeof(double));
    for(int i=0; i<size; i++){
        rand[i] = gsl_rng_uniform(r);
    } 
    // Define file to write to
    char *filename = "uniform.txt";
    FILE * fp;
    fp = fopen(filename, "w");
    // Write the vector to a file
    for (int i=0; i<size; i++){
        fprintf(fp, "%.2f\n", rand[i]);
    }
    // Free rng and vector
    free(rand);
    gsl_rng_free(r);
}   