// This file reads a 2D-array of size 10*10 from a file.

#include <stdio.h>
#include <stdlib.h>

int check_matrix_data(int rowSize, int* arr){
    int ok = 1; // This is zero if not okay, 1 if okay
    for(int ix = 0; ix<rowSize; ix++){
       for(int jx=0; jx<rowSize; jx++){
            int a = ix*jx;
            int b = arr[ix*rowSize+jx];
            ok = ok && a==b;
        
            
        }
    }
    return ok;
}

int
main(
        int argc,
        char* argv[]
    )
{


    long size = 10*10;
    printf("Number of ints in file: %ld for a total of %ld bytes\n", size, size*sizeof(int));

    // First read byte file since that is easier for our purposes
    FILE *fp;
    fp = fopen("./basic_assignment_bytes.txt", "rb");
    fseek(fp, 0, SEEK_END);
    int bytefile_length = ftell(fp);
    printf("File length of byte file: %li bytes \n", bytefile_length);
    fseek(fp, 0, SEEK_SET);

    int *asentries = (int*) malloc(sizeof(int)*size); // read file into one continous block of memory
    fread(asentries, size, sizeof(int), fp);
    

    printf("Reading of byte file passed? %i \n", check_matrix_data(10, asentries));
    fclose(fp);
    free(asentries);
    
    // Now read text file as well
    FILE *tp;
    tp = fopen("./basic_assignment_text.txt", "r");
    fseek(tp, 0, SEEK_END);
    int filelength = ftell(tp);
    printf("File length of text file: %li bytes\n", filelength);
    fseek(tp, 0, SEEK_SET);

    // Read all ints into a contigous block of memory as above
    int *bsentries = (int*) malloc(sizeof(int)*size);
    for(int ix=0; ix<size; ix++){
        int b;
        int r = fscanf(fp, "%i", &b);
        if(r!= 1)
            printf("Error when reading from file\n");
        bsentries[ix] = b;
    }

    printf("Reading of text file passed? %i \n", check_matrix_data(10, bsentries));

    return 0; 
}
