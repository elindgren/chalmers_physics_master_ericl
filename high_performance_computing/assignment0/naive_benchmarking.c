// This program calculates the sum of the first billion integers and prints the result

#include <stdio.h>
#include <time.h>

int
main(
    int argc,
    char *argv[]
)
{
    int it = 10; 
    long int iters = 1000000000; // no. of iterations
    clock_t t = clock();
    for(int j=0; j<it; j++) {
      unsigned long long sum = 0;
      // Start counting
      for(int ix=1; ix<=iters; ++ix)
          sum += ix;
    }

    // Calculate time taken
    t = clock() - t; 
    printf("t=%f\n", (double)t);
    double time_taken = ((double)t)/CLOCKS_PER_SEC;
   printf("Elapsed time: %fs - Time per iteration: %fs \n", time_taken, time_taken/it);
    return 0;
}
