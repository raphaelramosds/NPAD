FLAGS=-g -Wall -std=c99

EXES=mpi_hello mpi_hello1 mpi_hello2 mpi_vector mpi_vector_ndiv mpi_prefix_sum

all: $(EXES)

mpi_hello: 
	mpicc $(FLAGS) -o mpi_hello mpi_hello.c

mpi_hello1: 
	mpicc $(FLAGS) -o mpi_hello1 mpi_hello1.c

mpi_hello2: 
	mpicc $(FLAGS) -o mpi_hello2 mpi_hello2.c

mpi_vector:
	mpicc $(FLAGS) -o mpi_vector mpi_vector.c

mpi_vector_ndiv:
	mpicc $(FLAGS) -o mpi_vector_ndiv mpi_vector_ndiv.c

mpi_prefix_sum:
	mpicc $(FLAGS) -o mpi_prefix_sum mpi_prefix_sum.c

clean:
	rm $(EXES)