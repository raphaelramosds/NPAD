FLAGS=-fopenmp

EXES=vadd heat

vadd:
	gcc $(FLAGS) vadd.c -o vadd

heat:
	gcc $(FLAGS) heat.c -o heat

clean:
	rm $(EXES)
