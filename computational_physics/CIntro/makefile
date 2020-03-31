CC = gcc  # c-compiler
CFLAGS = -std=c99 # Compiler flags
LIBS = -lm -lgsl -lgslcblas

OBJECTS = c2_functions.o c_2.o
PROGRAM = c_2
# build .o files
%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)
# compile  program
all: $(PROGRAM)
$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)