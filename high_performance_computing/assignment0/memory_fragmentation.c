// This program plays around with two types of memory allocation;
// contigous and fragmented memory.

#include <stdlib.h>
#include <stdio.h>

int
main(
    int argc, 
    char* argv[]
){
    int size = 10; // size of the array

    int **as = (int**) malloc(sizeof(int*) * size); // allocating an array of pointers, which point to pointers to int (sub-arrays) - i.e. a 2D array.
    for (size_t ix=0; ix < size; ++ix)
        as[ix] = (int*) malloc(sizeof(int) * size); // fragmented memory - can end up anywhere on the heap

    for(size_t ix = 0; ix < size; ++ix)
        for(size_t jx=0; jx < size; ++jx)
            as[ix][jx] = 0; // allocate the array - jumpt to memory location of each sub array. Note that Each sub-array is contigous!
    printf("%d\n", as[0][0]);

    for(size_t ix=0; ix<size; ++ix)
        free(as[ix]); // We could've done this in one go - allocate as[ix], initialize it, print and then free it. It would accomplish the same thing since the memory is fragmented anyway.
    free(as);

    // Now do the same, but avoid memory fragmentation by allocating memory in one go
    int *bsentries = (int*) malloc(sizeof(int) * size*size); // allocating contigous memory for our arrays.  This will hold all our sub-arrays.
    int **bs = (int**) malloc(sizeof(int*) * size); // As before, list of pointers to pointers to int, 2D array.
    for(size_t ix=0, jx=0; ix<size; ++ix, jx+=size)
        bs[ix] = bsentries + jx; // the list of pointers, with pointers to int - fill as with the pointers to the beginning of each sub-array (size apart from each other in memory)

    for(size_t ix=0; ix<size; ++ix)
        for(size_t jx=0; jx <size; ++jx)
            bs[ix][jx] = 0; // Initialize 2D array with zeros

    printf("%d\n", bs[0][0]);
    
    free(bs);
    free(bsentries);

    return 0;
}
