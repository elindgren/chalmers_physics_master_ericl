// This file prints the input parameters

#include <stdio.h>
#include <string.h>

int
main(
    int argc,
    char* argv[]
    )
{   

    char avalue = '0';
    char bvalue = '0';

    for(int ix=0; ix<argc; ix++){
        char c  = argv[ix][1];
        char val = argv[ix][2];
        switch(c){
            case 'a':
                avalue = val;
                break; 
            case 'b': 
                bvalue = val;
                break;
        }
    }

    printf("A is %c and B is %c\n", avalue, bvalue);

    return 0;
}
