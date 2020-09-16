// This file writes a 2D-array of size 10*10 to a file.

#include <stdio.h>
#include <stdlib.h>

int
main(
        int argc,
        char* argv[]
    )
{
    int size = 10; // Array size
    int *asentries = (int*) malloc(sizeof(int)*size*size); // Allocate contigously on the heap
    int **as = (int**) malloc(sizeof(int*)*size);

    for(int ix=0; ix<size; ix++)
       as[ix] =  asentries + ix*size; // initialize array of pointers with pointers to int

    // Open file for writing
    FILE *fp;
    fp = fopen("./basic_assignment_text.txt", "w");
    FILE *bp;
    bp = fopen("./basic_assignment_bytes.txt", "w");

    for(int ix=0; ix<size; ix++){
       for(int jx=0; jx<size; jx++){
            as[ix][jx] = ix*jx; // Initialize 2D array
            // Writing file in text format - possible both with contigous and non-contigous memory
            fprintf(fp, "%i ", as[ix][jx]); 
       }                
       fprintf(fp, "\n");
    }
    
    // writing file as bytes - this is often faster, but requires the array as a contigous piece of memory (to be one array).
    fwrite(asentries, sizeof(int), size*size, bp);
    
    // Close file
    fclose(fp);
    fclose(bp);
    
    // Free arrays
    free(as);
    free(asentries);

    return 0;   
}
