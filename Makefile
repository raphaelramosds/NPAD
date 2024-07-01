FLAGS=-fopenmp -O3

EXES=vadd heat

vadd:
	gcc $(FLAGS) vadd.c -o vadd

heat:
	gcc $(FLAGS) heat.c -o heat -lm

clean:
	rm $(EXES)
