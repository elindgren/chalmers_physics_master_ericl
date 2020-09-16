# - GCC optimization options are described at
#   https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html
# - The flag -flto enables link time optimization.
# - The flag -fprofile-generate enables a binary that generates profiling
#   code, which can then be consumed with the flag -fprofile-use.
# - GCC options for x86 are available at
#   https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html
# - This is part of list of machine dependent options at
#   https://gcc.gnu.org/onlinedocs/gcc/Submodel-Options.html#Submodel-Options

BINS = O0 O1 O2 O3 O2nat O3nat O2g O3g

.PHONY : all
all : $(BINS)

.PHONY : clean
clean :
	rm -f $(BINS)


O0 : naive_benchmarking.c
	gcc -o $@ $< -O0 -g -march=x86-64

O1 : naive_benchmarking.c
	gcc -o $@ $< -O1 -g -march=x86-64

O2 : naive_benchmarking.c
	gcc -o $@ $< -O2 -g -march=x86-64

O3 : naive_benchmarking.c
	gcc -o $@ $< -O3 -g -march=x86-64

O2nat: naive_benchmarking.c
	gcc -o $@ $< -O2 -g -march=native

O3nat: naive_benchmarking.c
	gcc -o $@ $< -O3 -g -march=native

O2g : naive_benchmarking.c
	gcc -o $@ $< -O2 -g -march=x86-64 -g

O3g : naive_benchmarking.c
	gcc -o $@ $< -O3 -g -march=x86-64 -g

bench : all
	@echo
	@echo -n "O0: "
	@./O0
	@echo -n "O1: "
	@./O1
	@echo -n "O2: "
	@./O2
	@echo -n "O3: "
	@./O3
	@echo -n "O2nat: "
	@./O2nat
	@echo -n "O3nat: "
	@./O3nat
