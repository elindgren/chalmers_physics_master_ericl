// This program allocates memory on the stack, with varying sizes until it
// reaches a segfault. 

#include <stdio.h>
#include <stddef.h>

int 
main
(
 int argc,
 char* argv[]
)
{
    // This works
    //int size = 10;

    //int as[size];
    //for(size_t ix = 0; ix<size; ++ix)
    //    as[ix] = 0;
    //printf("%d\n", as[0]);

    // This, however, does not work since we allocate more memory
    // than what is physically available on the stack, resulting in a segmentation fault.
    int size = 1000000000;

    int as[size];
    for(size_t ix = 0; ix<size; ++ix)
        as[ix] = 0;
    printf("%d\n", as[0]);
    
    return 0;

}
