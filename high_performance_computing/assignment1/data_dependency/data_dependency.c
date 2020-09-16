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
  for(size_t ix=0; ix<nrs; ++ix){
    double sum = 0;
    for (size_t jx = 0; jx<ncs; ++jx)
      sum += matrix[ix][jx];
    sums[ix] = sum;
  }
}

static inline
void
row_sums_unrolled2(
  double * sums,
  const double **matrix, 
  size_t nrs, 
  size_t ncs
){
  for(size_t ix=0; ix<nrs; ++ix){
    double sum0 = 0;
    double sum1 = 0;
    for (size_t jx = 0; jx<ncs; jx += 2){
      sum0 += matrix[ix][jx];
      sum1 += matrix[ix][jx+1];
    }
    sums[ix] = sum0 + sum1;
  }
}

static inline
void
row_sums_unrolled4(
  double * sums,
  const double **matrix, 
  size_t nrs, 
  size_t ncs
){
  for(size_t ix=0; ix<nrs; ++ix){
    double sum0 = 0;
    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;
    for (size_t jx = 0; jx<ncs; jx += 4){
      sum0 += matrix[ix][jx];
      sum1 += matrix[ix][jx+1];
      sum2 += matrix[ix][jx+2];
      sum3 += matrix[ix][jx+3];
    }
    sums[ix] = sum0 + sum1 + sum2 + sum3;
  } 
}


static inline
void
row_sums_unrolled8(
  double * sums,
  const double **matrix, 
  size_t nrs, 
  size_t ncs
){
  for(size_t ix=0; ix<nrs; ++ix){
    double sum0 = 0;
    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;
    double sum4 = 0;
    double sum5 = 0;
    double sum6 = 0;
    double sum7 = 0;
    for (size_t jx = 0; jx<ncs; jx += 8){
      sum0 += matrix[ix][jx];
      sum1 += matrix[ix][jx+1];
      sum2 += matrix[ix][jx+2];
      sum3 += matrix[ix][jx+3];
      sum4 += matrix[ix][jx+4];
      sum5 += matrix[ix][jx+5];
      sum6 += matrix[ix][jx+6];
      sum7 += matrix[ix][jx+7];
    }
    sums[ix] = sum0 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7;
  } 
}

static inline
void
row_sums_unrolled4_array(
  double * sums,
  const double **matrix, 
  size_t nrs, 
  size_t ncs
){
  for(size_t ix=0; ix<nrs; ++ix){
    double sum[4] = {0,0,0,0};
    for (size_t jx = 0; jx<ncs; ++jx){
      sum[0] += matrix[ix][jx];
      sum[1] += matrix[ix][jx+1];
      sum[2] += matrix[ix][jx+2];
      sum[3] += matrix[ix][jx+3];
    }
    sums[ix] = sum[0] + sum[1] + sum[2] + sum[3];
  }
}



int
main(
  int argc,
  char *argv[]
){
  int benchmark_iters = 5000;
  struct timespec bench_start_time;
  struct timespec bench_stop_time;
  srand(0); // Set a fixed seed
  int rand_idx;

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

  //********* Perform benchmark - for each type of unrolled loop
  timespec_get(&bench_start_time, TIME_UTC);
  for(int i=0; i<benchmark_iters; i++)
    row_sums(rsums, matrix, nrs, ncs);
  timespec_get(&bench_stop_time, TIME_UTC);
  double bench_baseline = difftime(bench_stop_time.tv_sec, bench_start_time.tv_sec)*1000 + (bench_stop_time.tv_nsec - bench_start_time.tv_nsec) /1000000;
  rand_idx = rand() % size;
  printf("Random value from rsums: %f \n", rsums[rand_idx]);

  timespec_get(&bench_start_time, TIME_UTC);
  for(int i=0; i<benchmark_iters; i++)
    row_sums_unrolled2(rsums, matrix, nrs, ncs);
  timespec_get(&bench_stop_time, TIME_UTC);
  double bench_unrolled2 = difftime(bench_stop_time.tv_sec, bench_start_time.tv_sec)*1000 + (bench_stop_time.tv_nsec - bench_start_time.tv_nsec) /1000000;
  rand_idx = rand() % size;
  printf("Random value from rsums: %f \n", rsums[rand_idx]);
  
  timespec_get(&bench_start_time, TIME_UTC);
  for(int i=0; i<benchmark_iters; i++)
    row_sums_unrolled4(rsums, matrix, nrs, ncs);
  timespec_get(&bench_stop_time, TIME_UTC);
  double bench_unrolled4 = difftime(bench_stop_time.tv_sec, bench_start_time.tv_sec)*1000 + (bench_stop_time.tv_nsec - bench_start_time.tv_nsec) /1000000;
  rand_idx = rand() % size;
  printf("Random value from rsums: %f \n", rsums[rand_idx]);
  
  timespec_get(&bench_start_time, TIME_UTC);
  for(int i=0; i<benchmark_iters; i++)
    row_sums_unrolled8(rsums, matrix, nrs, ncs);
  timespec_get(&bench_stop_time, TIME_UTC);
  double bench_unrolled8 = difftime(bench_stop_time.tv_sec, bench_start_time.tv_sec)*1000 + (bench_stop_time.tv_nsec - bench_start_time.tv_nsec) /1000000;
  rand_idx = rand() % size;
  printf("Random value from rsums: %f \n", rsums[rand_idx]);

  timespec_get(&bench_start_time, TIME_UTC);
  for(int i=0; i<benchmark_iters; i++)
    row_sums_unrolled4_array(rsums, matrix, nrs, ncs);
  timespec_get(&bench_stop_time, TIME_UTC);
  double bench_unrolled4_array = difftime(bench_stop_time.tv_sec, bench_start_time.tv_sec)*1000 + (bench_stop_time.tv_nsec - bench_start_time.tv_nsec) /1000000;
  rand_idx = rand() % size;
  printf("Random value from rsums: %f \n", rsums[rand_idx]);

  // benchmark over
  printf("\nTimings for summations:: \n \tbaseline: %fmus\n \tUnrolled-2: %fmus\n \tUnrolled-4: %fmus\n \tUnrolled-8: %fmus\n \tUnrolled-4-array: %fmus\n\n", bench_baseline/benchmark_iters, bench_unrolled2/benchmark_iters, bench_unrolled4/benchmark_iters, bench_unrolled8/benchmark_iters, bench_unrolled4_array/benchmark_iters);

  free(mentries);
  free(matrix);
  free(rsums);
  free(csums);
}
