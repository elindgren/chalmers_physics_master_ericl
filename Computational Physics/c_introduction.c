# include <stdio.h>
# include <stdlib.h>
# include <math.h>
int main() {
    double y = sqrt(2);
    if (y<1){
        printf("y is greater than 1!");
    }
    else{
        printf("Errorrorror");
    }
    printf("sqrt(y)=%.2f\n", y);
    int *ypointer = &y;
    printf("address for y is %d\n", *ypointer);
    return 0;
}