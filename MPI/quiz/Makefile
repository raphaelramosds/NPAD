FLAGS=-g -Wall -std=c99

EXES=mpi_hello mpi_hello_max_string mpi_hello_without_eof mpi_vector mpi_vector_ndiv mpi_prefix_sum mpi_odd_even_time

all: $(EXES)

mpi_hello: 
	mpicc $(FLAGS) -o mpi_hello mpi_hello.c

mpi_hello_max_string: 
	mpicc $(FLAGS) -o mpi_hello_max_string mpi_hello_max_string.c

mpi_hello_without_eof: 
	mpicc $(FLAGS) -o mpi_hello_without_eof mpi_hello_without_eof.c

mpi_vector:
	mpicc $(FLAGS) -o mpi_vector mpi_vector.c

mpi_vector_ndiv:
	mpicc $(FLAGS) -o mpi_vector_ndiv mpi_vector_ndiv.c

mpi_prefix_sum:
	mpicc $(FLAGS) -o mpi_prefix_sum mpi_prefix_sum.c

mpi_odd_even_time:
	mpicc $(FLAGS) -o mpi_odd_even_time mpi_odd_even_time.c

clean:
	rm $(EXES) *.out *.err
