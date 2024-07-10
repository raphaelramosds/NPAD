all: LU

LU: LU_rand_omp.cpp
	g++ -O3 -ffast-math -Wall -fopenmp -o LU LU_rand_omp.cpp -lm -lmpascalops
clean:
	rm -f LU *.o

