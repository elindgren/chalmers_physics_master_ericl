################################################################################
# FKA121/FIM540 Computational Physics 2019
# C programs for exercise E1
#
# Martin Gren, martin.gren@chalmers.se 
################################################################################
If you have installed the gcc compiler the program can be compiled with the 
following command:

			gcc -o run main.c -O3 -W -lm

which will create an executable file called "run".

To run the program type './run' in the terminal. 

The data is dumped into files with suffix '.dat'. 
The program plot_function.py plots the function in Python.
To plot using python3 you first need to install numpy and matplotlib: 

	python3 -m pip install --user numpy matplotlib 

This only need to be done once! 
Now you can plot by typing:

	python3 plot_function.py

The plot will be saved in .pdf format and displayed in an interactive window.
To stop the program from displaying the plot you need to comment plt.show() 
in plot_function.py.

To open the pdf from the command line you can type 

	evince h_t.pdf &
