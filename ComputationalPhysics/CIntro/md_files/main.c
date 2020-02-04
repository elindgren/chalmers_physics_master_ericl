#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "initfcc.h"
#include "alpotential.h"

int main()
{
    // Define lattice parameter
    double lat_param = 1.0;
    // length of positions
    int N = 100;
    // Declare 2D array to hold positions
    double pos[4*N*N*N][3];
    //pos = malloc(4*N*N*N * 3 *sizeof(double));
    // Pre-allocate memory
    // for (int i=0; i<4*N*N*N; i++){
    //     pos[i] = (double *)malloc(3*sizeof(double));
    // }
    // Send 2D-array to initfcc along with lattice parameter
    init_fcc(pos, N, lat_param);
    double energy = get_energy_AL(pos, lat_param, N);
    printf("Energy is: %.2f", energy);
}
