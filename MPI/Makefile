FLAGS=-g -Wall -std=c99

EXES=mpi_hello mpi_trap1 mpi_trap1_any_source mpi_type_vect mpi_odd_even mpi_mat_triangle

all: $(EXES)

mpi_hello: 
	mpicc $(FLAGS) -o mpi_hello mpi_hello.c

mpi_trap1:
	mpicc $(FLAGS) -o mpi_trap1 mpi_trap1.c

mpi_trap1_any_source:
	mpicc $(FLAGS) -o mpi_trap1_any_source mpi_trap1_any_source.c

mpi_type_vect:
	mpicc $(FLAGS) -o mpi_type_vect mpi_type_vect.c

mpi_odd_even:
	mpicc $(FLAGS) -o mpi_odd_even mpi_odd_even.c

mpi_mat_triangle:
	mpicc $(FLAGS) -o mpi_mat_triangle mpi_mat_triangle.c

clean:
	rm $(EXES) *.out