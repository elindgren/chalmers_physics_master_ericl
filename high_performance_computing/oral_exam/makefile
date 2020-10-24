.PHONY : all

all : clean overhead naive handin optimized_ls optimized_mem

overhead : heat_overhead.c parsing_cmdl_args.c read_file.c utils.c
	gcc -o overhead heat_overhead.c parsing_cmdl_args.c read_file.c utils.c -O3 -lm -march=native -lOpenCL

naive : heat_naive.c parsing_cmdl_args.c read_file.c utils.c
	gcc -o naive heat_naive.c parsing_cmdl_args.c read_file.c utils.c -O3 -lm -march=native -lOpenCL

handin : heat_handin.c parsing_cmdl_args.c read_file.c utils.c
	gcc -o handin heat_handin.c parsing_cmdl_args.c read_file.c utils.c -O3 -lm -march=native -lOpenCL

optimized_ls : heat_optimized_ls.c parsing_cmdl_args.c read_file.c utils.c
	gcc -o optimized_ls heat_optimized_ls.c parsing_cmdl_args.c read_file.c utils.c -O3 -lm -march=native -lOpenCL

optimized_mem : heat_optimized_mem.c parsing_cmdl_args.c read_file.c utils.c
	gcc -o optimized_mem heat_optimized_mem.c parsing_cmdl_args.c read_file.c utils.c -O3 -lm -march=native -lOpenCL
clean : 
	-rm -f *.o naive handin optimized_ls optimized_mem
