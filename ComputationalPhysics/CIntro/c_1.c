# include <stdio.h>
# include <stdlib.h>
# include <math.h>

void problem3() {
    int m = 2;
    int n = 5;
    double x = 5.6;
    double y = m;
    int k = x;
    int k1 = x-n;

    printf("m/n=%d\n", m/n);
    printf("x/n=%.2f\n", x/n);
    printf("double y=m=%.2f\n", y);
    printf("int k=x=%d\n", k);
    printf("int k1=x-n=%d\n", k1);
    
}

void problem4(){
    int m = 5;
    int *pm = &m;
    int n = pm;
    printf("m+1=%d\n", m+1);
    printf("*pm+1=%d\n", *pm+1);
    printf("pm=%x\n", pm);
    printf("&m=%x\n", &m+1);
    printf("*(&m)+1=%d\n", *(&m)+1);
    printf("int n=pm=%d\n", n);
}

double problem5(double *x, double *y){
    return pow(*x, 3) + sqrt(*y);
}

void increment(double *x, double val){
    *x += val;
}

void print_array(double arr[], int len){
    printf("[");
    for (int i=0; i<len; i+=1){
        if (i !=len-1){
            printf("%.2f, ", arr[i]);
        }else{
            printf("%.2f", arr[i]);
        }
        
    }
    printf("]");
}

void write_array(char *file_name[], double arr[], int len){
    FILE *outputfile = fopen(*file_name,"w");    
    fprintf(outputfile,"[");
    for(int i=0; i<len; i+=1){
        if(i !=len-1){
            fprintf(outputfile, "%.2f, ", arr[i]);
        }else{
            fprintf(outputfile, "%.2f", arr[i]);
        }
    }
    fprintf(outputfile, "]");
}

void matrix_multiply(double A[3][3], double B[3][3], double C[3][3]){
    for(int i=0;i<3;i+=1){
        for(int j=0; j<3; j+=1){
            double sum = 0.0;
            for(int k=0; k<3; k+=1){
                sum += A[i][k]*B[k][j];
            }
            C[i][j] = sum;
        }
    }
}

int main() {
    printf("%s########## Problem 3 #########\n");
    problem3();
    printf("\n");
    printf("########## Problem 4 #########\n");
    problem4();
    printf("\n");
    printf("########## Problem 5 #########\n");
    double x = 2.0;
    double x_old = x;
    double y = 9.0;
    printf("Problem 5a: %.2f\n", problem5(&x, &y));
    increment(&x, 1.0);
    printf("Problem 5b: x=%.2f, x+1=%.2f\n", x_old, x);
    printf("\n");
    printf("########## Problem 6 #########\n");
    double arr[] = {1.0, 2.5, 3.1, 4.9};
    char* file_name = "c_introduction.txt\0";
    print_array(arr, 4);
    write_array(&file_name, arr, 4);
    printf("\n");
    double A[3][3] = {{1,2,3}, {4,5,6}, {7,8,9}};
    double B[3][3] = {{1,1,1}, {1,1,1}, {1,1,1}};
    double C[3][3]; 
    matrix_multiply(A, B, C);
    print_array(C[0], 3);
}

