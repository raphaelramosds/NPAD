FLAGS=-fopenmp -O3

EXES=vadd vadd_p heat

all: $(EXES)

vadd:
	gcc $(FLAGS) vadd.c -o vadd

vadd_p:
	gcc $(FLAGS) vadd_p.c -o vadd_p

heat:
	gcc $(FLAGS) heat.c -o heat -lm

clean:
	rm $(EXES)
