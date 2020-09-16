#include <stdio.h>
#include <stdlib.h>

// Commenting out the initialization of as: Valgrind complains that I'm trying to access uninitalized memory - no memory leak, but the program encounters a segmentation fault and ends.
// Comment out the freeing of as: Valgrind notices that there is a memory leak of 40 blocks (bytes) due to not freeing. 
// Amend the code with an additional free(as): Valgrind says that we have already freed the memory in question, and that we are not allowed to do so. 

int
main(
  int argc,
  char* argv[]
)
{
  int *as;
 as = (int*) malloc(10*sizeof(int));
  int sum = 0;

  for(int ix=0; ix<10; ix++)
    as[ix] = ix;

  for(int ix=0; ix<0; ix++)
    sum += as[ix];

  printf("%d\n", sum);
 free(as);
  free(as);
  return 0;
}

