.PHONY : all

all : clean cell_distances

cell_distances : cell_distances.c utils.c parsing_cmdl_arg.c
	gcc -o cell_distances cell_distances.c utils.c parsing_cmdl_arg.c -O2 -lm -fopenmp -march=native 

clean : 
	-rm -f *.o cell_distances
