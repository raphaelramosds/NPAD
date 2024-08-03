FLAGS=-g -Wall -std=c99

EXES=mpi_hello mpi_hello_max_string mpi_hello_without_eof mpi_vector mpi_vector_ndiv mpi_prefix_sum mpi_odd_even_time mpi_type_indexed mpi_trap_pack mpi_pp_clock mpi_pp_wtime mpi_trap_time mpi_mat_vect_cols mpi_mat_vect_lines

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

mpi_type_indexed:
	mpicc $(FLAGS) -o mpi_type_indexed mpi_type_indexed.c

mpi_trap_pack:
	mpicc $(FLAGS) -o mpi_trap_pack mpi_trap_pack.c

mpi_ping_pong:
	mpicc $(FLAGS) -o mpi_pp_clock mpi_ping_pong.c
	mpicc -DWTIME $(FLAGS) -o mpi_pp_wtime mpi_ping_pong.c

mpi_trap_time:
	mpicc $(FLAGS) -o mpi_trap_time mpi_trap_time.c

mpi_mat_vect_cols:
	mpicc $(FLAGS) -o mpi_mat_vect_cols mpi_mat_vect_cols.c

mpi_mat_vect_lines:
	mpicc $(FLAGS) -o mpi_mat_vect_lines mpi_mat_vect_lines.c

clean:
	rm $(EXES) *.out *.err
