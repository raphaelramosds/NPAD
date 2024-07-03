FLAGS=-fopenmp -O3

EXES=vadd vadd_p vadd_p_gpu heat heat_gpu

all: $(EXES)

vadd:
	gcc $(FLAGS) vadd.c -o vadd

vadd_p:
	gcc $(FLAGS) vadd_p.c -o vadd_p

vadd_p_gpu:
	gcc $(FLAGS) vadd_p_gpu.c -o vadd_p_gpu

heat:
	gcc $(FLAGS) heat.c -o heat -lm

heat_gpu:
	gcc $(FLAGS) heat_gpu.c -o heat_gpu -lm

clean:
	rm $(EXES) *.out
