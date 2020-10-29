.PHONY : all

all : clean heat_diffusion

heat_diffusion : heat_diffusion.c parsing_cmdl_args.c read_file.c utils.c
	gcc -o heat_diffusion heat_diffusion.c parsing_cmdl_args.c read_file.c utils.c -O3 -lm -march=native -lOpenCL

clean : 
	-rm -f *.o heat_diffusion
