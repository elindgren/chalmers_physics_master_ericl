#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// naive summation implementations
static inline
void
row_sums(
  double * sums,
  const double **matrix, 
  size_t nrs, 
  size_t ncs
){
  // This is easy to speed up - since the memory is contigous
  double sum = 0;
  int row = 0;
  for(size_t ix=0; ix<nrs*ncs; ++ix){
    sum += *(matrix[0]+ix);
    if(ix%nrs == 0 && ix>0){
      sums[row] = sum;
      row += 1;
      sum = 0;
    }
  }
  for(size_t ix=0; ix<nrs; ++ix){
    double sum = 0;
    for (size_t jx = 0; jx<ncs; ++jx)
      sum += matrix[ix][jx];
    sums[ix] = sum;
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
  // column summation is the slow one - implement it more efficiently
  // TODO: How to speed this up? Still row-major matrix. 
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
  int benchmark_iters = 5000;
  struct timespec bench_start_time_row;
  struct timespec bench_start_time_col;
  struct timespec bench_stop_time_row;
  struct timespec bench_stop_time_col;

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

  
  //********* Perform benchmark - first of row summation, then col summation
  timespec_get(&bench_start_time_row, TIME_UTC);
  for(int i=0; i<benchmark_iters; i++)
    row_sums(rsums, matrix, nrs, ncs);
  timespec_get(&bench_stop_time_row, TIME_UTC);

  timespec_get(&bench_start_time_col, TIME_UTC);
  for(int i=0; i<benchmark_iters; i++)
    col_sums(csums, matrix, nrs, ncs); 
  timespec_get(&bench_stop_time_col, TIME_UTC);

  double bench_row_time = difftime(bench_stop_time_row.tv_sec, bench_start_time_row.tv_sec)*1000 + (bench_stop_time_row.tv_nsec - bench_start_time_row.tv_nsec) /1000000;
  double bench_col_time = difftime(bench_stop_time_col.tv_sec, bench_start_time_col.tv_sec)*1000 + (bench_stop_time_col.tv_nsec - bench_start_time_col.tv_nsec) /1000000;
  
  printf("Timings for summations:: rows: %fmus ---- cols: %fmus \n", bench_row_time/benchmark_iters, bench_col_time/benchmark_iters);
  //******** benchmark over
  //
  srand(0); // Set a fixed seed
  int rand_idx = rand() % size;
  printf("Random value from: csums: %f --- rsums: %f \n", csums[rand_idx], rsums[rand_idx]);

  free(mentries);
  free(matrix);
  free(rsums);
  free(csums);
}
