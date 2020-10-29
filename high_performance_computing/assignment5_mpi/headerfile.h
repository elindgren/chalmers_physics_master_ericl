void parsing_cmdl_args(int argc, char** argv, int* nbr_iteration, float* d);
void read_header(FILE* file, int *nx, int *ny);
void read_data(FILE* file, int nx, int up, float* h_old);
float calc_average(float* h, int nx, int ny,int up);
float calc_distance_from_average(float* h, int nx, int ny, float average,int up);
