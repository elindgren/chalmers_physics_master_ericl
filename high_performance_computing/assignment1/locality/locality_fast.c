#include <stdio.h>
#include <stdlib.h>

// faster summation implementation
static inline
void
row_sums(
  double * sums,
  const double **matrix, 
  size_t nrs, 
  size_t ncs
){
  double sum = 0;
  int row = 0;
  for(size_t ix=0; ix<nrs*ncs; ++ix){
    sum += *(matrix[0] + ix); // we know that we have all the memory after the first entry in matrix
    if(ix%nrs == 0 && ix>0){
      sums[row] = sum;
      row += 1;
      sum = 0;
    }
  }
}

static inline
void
col_sums(
  double *sums,
  const double **matrix,
  size_t nrs,
  size_t ncs
){
  //Since the matrix is stored in row major form I don't think we can roll out
  //the loop here as in the row-optimized version. 
  for (size_t jx=0; jx<ncs; jx++){
    double sum = 0;
    for(size_t ix=0; ix<nrs; ix++)
      sum += matrix[ix][jx];
    sums[jx] = sum;
  }
}

int
main(
  int argc,
  char *argv[]
){
  // declare matrix
  int benchmark_iters = 5000;
  size_t size = 1000;
  size_t nrs = size;
  size_t ncs = size; 
  double *mentries = (double*) malloc(sizeof(double)*nrs*ncs); //matrix entries
  const double **matrix = (const double**) malloc(sizeof(double*)*nrs); // saved in row major, i.e. rows after each other 
  double *rsums = (double*) malloc(sizeof(double)*nrs);
  double *csums = (double*) malloc(sizeof(double)*ncs);

  // Initialize matrix entries
  for(int i=0; i<size*size; i++)
    mentries[i] = i%10;

  //Initialize matrix with pointers to mentries (rows)
  for(int i=0; i<size; i++)
    matrix[i] = mentries + i*size;

  // Perform benchmark
  for(int i=0; i<benchmark_iters; i++){
    row_sums(rsums, matrix, nrs, ncs);
    col_sums(csums, matrix, nrs, ncs); 
  }
  srand(0); // Set a fixed seed
  int rand_idx = rand() % size;
  printf("Random value from: csums: %f --- rsums: %f \n", csums[rand_idx], rsums[rand_idx]);

  free(mentries);
  free(matrix);
  free(rsums);
  free(csums);
}
