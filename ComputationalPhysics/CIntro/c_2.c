#include <stdio.h>
#include <math.h>
#include <gsl/gsl_blas.h>

void problem1_scalar_product(gsl_vector A[], gsl_vector B[]){
    double res = 0.0;
    gsl_blas_ddot(A, B, &res);
    printf("Result is: %.2f. \n", res);
}

int main(){
    int len = 0;
    printf("Enter array length: ");
    scanf("%d", &len);  // Unsigned int
    printf("Scanned length: %d. \n", len);
    // Declare arrays
    gsl_vector *A; 
    gsl_vector *B;
    // Allocate memory
    A = gsl_vector_alloc(len);
    B = gsl_vector_alloc(len);
    // Pre-allocate arrays
    gsl_vector_set_all(A, 1.0);
    gsl_vector_set_all(B, 2.0);
    // Calculate scalar product
    problem1_scalar_product(A,B);
}