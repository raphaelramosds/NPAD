FLAGS=-g -Wall -std=c99

EXES=mpi_hello mpi_hello_1 mpi_hello_2

all: $(EXES)

mpi_hello: 
	mpicc $(FLAGS) -o mpi_hello mpi_hello.c

mpi_hello_1: 
	mpicc $(FLAGS) -o mpi_hello_1 mpi_hello_1.c

mpi_hello_2: 
	mpicc $(FLAGS) -o mpi_hello_2 mpi_hello_2.c

clean:
	rm $(EXES)