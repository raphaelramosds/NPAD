FLAGS=-fopenmp -O3

EXES=vadd vadd_p vadd_p_gpu heat

all: $(EXES)

vadd:
	gcc $(FLAGS) vadd.c -o vadd

vadd_p:
	gcc $(FLAGS) vadd_p.c -o vadd_p

vadd_p_gpu:
	gcc $(FLAGS) vadd_p_gpu.c -o vadd_p_gpu

heat:
	gcc $(FLAGS) heat.c -o heat -lm

clean:
	rm $(EXES) *.out
