FLAGS=-g -Wall -std=c99

EXES=mpi_hello mpi_trap1 mpi_trap1_as

all: $(EXES)

mpi_hello: 
	mpicc $(FLAGS) -o mpi_hello mpi_hello.c

mpi_trap1:
	mpicc $(FLAGS) -o mpi_trap1 mpi_trap1.c

mpi_trap1_as:
	mpicc $(FLAGS) -o mpi_trap1_as mpi_trap1_as.c

clean:
	rm $(EXES) *.out