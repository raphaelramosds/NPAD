FLAGS=-g -Wall -std=c99

EXES=mpi_hello mpi_trap1

all: $(EXES)

mpi_hello: 
	mpicc $(FLAGS) -o mpi_hello mpi_hello.c

mpi_trap1:
	mpicc $(FLAGS) -o mpi_trap1 mpi_trap1.c

clean:
	rm $(EXES) *.out