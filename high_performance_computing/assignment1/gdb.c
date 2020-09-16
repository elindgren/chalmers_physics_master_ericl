#include <stdio.h>
#include <stdlib.h>

// When printing p on line 13 we get that it is a pointer to adress 0x0 - i.e. a NULL pointer. Hence we get the segmentation fault since we try to access memory that we haven't allocated. 

int
main(
    int argc,
    char* argv[]
){
  int *as = NULL;
  int sum = 0;

  for(int ix=0; ix<10; ix++)
    as[ix] = ix;

  for(int ix=0; ix<10; ix++)
    sum += as[ix];
  free(as);

  return 0;
}
