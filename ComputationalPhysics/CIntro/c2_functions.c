#include <stdio.h>
#include <math.h>
#include <gsl/gsl_blas.h>


void problem1_scalar_product(double A[], double B[], int len){
    double res = 0.0;
    for (int i=0; i<len; i=i+1){
        res = res + A[i] * B[i];
    }
    printf("Result is: %.2f. \n", res);
}

void calc_distance(int p1, int p2, double ps[][3]){
    double dist = 0.0;
    for (int i=0; i<3; i++){
        dist = dist + pow(ps[p1][i]-ps[p2][i], 2.0);
    }
    printf("Distance is: %.2f. \n", dist);
}