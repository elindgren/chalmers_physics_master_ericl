// This program allocates memory on the heap, with sizes that would encounter a segmentation fault on the heap.

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

int 
main
(
 int argc,
 char* argv[]
)
{
    // This works now, but didn't work before
    int size = 1000000000;
    int *as = (int*) malloc(sizeof(int)*size);
    for(size_t ix = 0; ix<size; ++ix)
        as[ix] = 0;
    printf("%d\n", as[0]);

    free(as);

    return 0;
}
