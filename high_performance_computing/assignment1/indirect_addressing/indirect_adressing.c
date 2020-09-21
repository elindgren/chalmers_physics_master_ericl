#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void
write_y_indirect(int *y, int *x, int *p,int a, int size){
  for(size_t kx=0; kx<size; ++kx){
    size_t jx = p[kx];
    y[jx] += a * x[jx];
  }
}

void
write_y_direct(int *y, int *x, int *p,int a, int size){
  for(size_t kx=0; kx<size; ++kx){
    y[kx] += a * x[kx];
  }
}

void
initialize_xy(int *x, int *y, int size){
  for(size_t i=0; i<size; i++){
    x[i] = 1;
    y[i] = 2;
  }
}


int
main(
  int argc,
  char *argv[]
){
  int benchmark_iters = 1000;
  size_t size = 1000000;
  size_t size_jump = 1000;
  int a = 10;
  int rand_idx;
  struct timespec bench_start_time;
  struct timespec bench_stop_time;
  srand(0); // Set a fixed seed

  // Allocate arrays
  int *x = (int*) malloc(sizeof(int)*size);
  int *y = (int*) malloc(sizeof(int)*size);
  int *p = (int*) malloc(sizeof(int)*size);

  initialize_xy(x,y, size);

  //********************* INDIRECT ACCESS OF X AND Y ************
  // Linear intialization
  timespec_get(&bench_start_time, TIME_UTC);
  for(int i=0; i<benchmark_iters; i++){
    for(size_t ix=0; ix<size; ++ix)
      p[ix] = ix;
    write_y_indirect(y, x, p, a, size);
  }
  timespec_get(&bench_stop_time, TIME_UTC);
  double linear_time = difftime(bench_stop_time.tv_sec, bench_start_time.tv_sec)*1000 + (bench_stop_time.tv_nsec - bench_start_time.tv_nsec) /1000000;
  rand_idx = rand() % size;
  printf("Random value from y: %d \n", y[rand_idx]);
  initialize_xy(x,y, size); // Reinitialize

  // Jumping intialization
  timespec_get(&bench_start_time, TIME_UTC);
  for(int i=0; i<benchmark_iters; i++){
    for(size_t jx=0, kx=0; jx<size_jump; ++jx)
      for(size_t ix = jx; ix<size; ix+=size_jump, kx++)
        p[ix] = kx;
    // jx=0 -> ix: 0, 1000, 2000, 3000 ... jx=1 -> ix: 1, 1001, 2001, 3001 ...
    // kx: 0, 1, 2, 3,....
    write_y_indirect(y, x, p, a, size);
  }
  timespec_get(&bench_stop_time, TIME_UTC);
  double jumping_time = difftime(bench_stop_time.tv_sec, bench_start_time.tv_sec)*1000 + (bench_stop_time.tv_nsec - bench_start_time.tv_nsec) /1000000;
  rand_idx = rand() % size;
  printf("Random value from y: %d \n", y[rand_idx]);
  initialize_xy(x,y, size); // Reinitialize


  //********************* ACCESS OF X AND Y ************
  // No initialization necessary
  initialize_xy(x,y, size); // Reinitialize
  timespec_get(&bench_start_time, TIME_UTC);
  for(int i=0; i<benchmark_iters; i++){
    write_y_direct(y, x, p, a, size);
  }
  timespec_get(&bench_stop_time, TIME_UTC);
  double time_direct = difftime(bench_stop_time.tv_sec, bench_start_time.tv_sec)*1000 + (bench_stop_time.tv_nsec - bench_start_time.tv_nsec) /1000000;
  rand_idx = rand() % size;
  printf("Random value from y: %d \n", y[rand_idx]);



  printf("\nTimings::\n  \tLinear initialization: %fmus\n\tJumping initialization: %fmus \n\t Direct access: %fmus\n", linear_time/benchmark_iters, jumping_time/benchmark_iters, time_direct/benchmark_iters);
  //******** benchmark over
  //


  free(x); free(y); free(p);
}
