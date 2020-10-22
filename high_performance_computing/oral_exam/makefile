.PHONY : all

all : clean naive handin optimized

naive : heat_naive.c parsing_cmdl_args.c read_file.c utils.c
	gcc -o naive heat_naive.c parsing_cmdl_args.c read_file.c utils.c -O3 -lm -march=native -lOpenCL

handin : heat_handin.c parsing_cmdl_args.c read_file.c utils.c
	gcc -o handin heat_handin.c parsing_cmdl_args.c read_file.c utils.c -O3 -lm -march=native -lOpenCL

optimized : heat_optimized.c parsing_cmdl_args.c read_file.c utils.c
	gcc -o optimized heat_optimized.c parsing_cmdl_args.c read_file.c utils.c -O3 -lm -march=native -lOpenCL

clean : 
	-rm -f *.o naive handin optimized
