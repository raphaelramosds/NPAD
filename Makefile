FLAGS=-fopenmp

EXES=vadd

vadd:
	gcc $(FLAGS) vadd.c -o vadd

clean:
	rm $(EXES)
