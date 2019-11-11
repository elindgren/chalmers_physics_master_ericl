################################################################################
# FKA121/FIM540 Computational Physics 2019
# C programs for exercise E1
#
# Martin Gren, martin.gren@chalmers.se 
################################################################################
The following commands will build the executable powerspectrum:

	gcc -c -o fft_func.o fft_func.c -O3 -W
	gcc -c -o main.o main.c -O3 -W
	gcc -o powerspectrum fft_func.o main.o -O3 -W -lm -lgsl -lgslcblas

Already with two .c files the compilation starts to become a bit messy.
One therefore uses makefiles that automatically executes the commands 
needed to compile the program. If the makefile is called Makefile or 
makefile one can type the command 'make' to compile the code. 

To run the program type './powerspectrum' in the terminal. 

The data is dumped into files with suffix '.dat'. 
The program plot_powerspectrum.py plots the powerspectrum in Python.
Plot by typing:

	python3 plot_powerspectrum.py

The plot will be saved in .pdf format and displayed in an interactive window.
To stop the program from displaying the plot you need to comment plt.show() 
in plot_function.py.

To open the pdf from the command line you can type 

	evince powerspectrum.pdf &

