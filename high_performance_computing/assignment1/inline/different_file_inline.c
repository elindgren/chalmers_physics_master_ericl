#include <stdio.h>
#include <stdlib.h>
#include <time.h>

inline void
mul_cpx(
    double *c_re,
    double *c_im,
    double *a_re,
    double *a_im,
    double *b_re,
    double *b_im
    ){
  // Calculate the complex multiplication of b and c, saving the parts in a_re and a_im each
  // a_re = b_re*c_re - b_im*c_im, a_im = b_re*c_im + b_im*c_re
  *a_re = *b_re * *c_re - *b_im * *c_im;
  *a_im = *b_re * *c_im + *b_im * *c_re;
}

int main(
    int argc,
    char *argv[]
)
{
  int benchmark_iters = 200000;
  int size = 30000;
  struct timespec bench_start_time;
  struct timespec bench_stop_time;
  double bench_diff_time;

  double *as_re = (double*) malloc(sizeof(double*)*size);
  double *as_im = (double*) malloc(sizeof(double*)*size);

  double *bs_re = (double*) malloc(sizeof(double*)*size);
  double *bs_im = (double*) malloc(sizeof(double*)*size);

  double *cs_re = (double*) malloc(sizeof(double*)*size);
  double *cs_im = (double*) malloc(sizeof(double*)*size);
  
  // Initialize arrays with values:
  for(int ix=0; ix<size; ix++){
    as_re[ix] = ix%2;
    as_im[ix] = ix%3+1;
    bs_re[ix] = ix%4+2;
    bs_im[ix] = ix%5+3;
    cs_re[ix] = ix%6+4;
    cs_im[ix] = ix%7+5;
  }
  printf("Before --- a_re: %f - a_im: %f \n", as_re[666], as_im[666]); // printing a random value
  
  // Now start the benchmarking. Perform all of the multiplications
  timespec_get(&bench_start_time, TIME_UTC);
  for(int bench=0; bench<benchmark_iters; bench++)
    for(int ix=0; ix<size; ix++)
      mul_cpx(cs_re+ix, cs_im+ix, as_re+ix, as_im+ix, bs_re+ix, bs_im+ix);
  // Print a random value of as afterwards
 printf("After --- a_re: %f - a_im: %f \n", as_re[666], as_im[666]); // printing a random value

  timespec_get(&bench_stop_time, TIME_UTC);
  bench_diff_time =
      difftime(bench_stop_time.tv_sec, bench_start_time.tv_sec)*1000 + (bench_stop_time.tv_nsec - bench_start_time.tv_nsec) /1000000;
 
  printf("Execution with automatic inlining (same file) took: %fms \n", bench_diff_time);

  free(as_re); free(as_im);
  free(bs_re); free(bs_im);
  free(cs_re); free(cs_im);
  return 0;
}