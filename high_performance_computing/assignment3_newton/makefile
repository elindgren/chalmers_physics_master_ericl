.PHONY : all

all : clean newton

newton : newton.c parsing_cmdl_arg.c newton_iteration.c write_to_file.c
	gcc -o newton newton.c parsing_cmdl_arg.c newton_iteration.c write_to_file.c -O2 -lm -march=native -ffast-math -lpthread 

clean : 
	-rm -f *.o newton
