#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
    int i, count = 5;
    int comm_sz, my_rank;
    int *x, *partial;
    MPI_Comm comm;

    x = malloc(count * sizeof(int));
    partial = malloc(count * sizeof(int));

    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);

    srand(time(NULL) + my_rank);

    for (i = 0; i < count; i++) x[i] = rand() % 10;

    MPI_Scan(x, partial, count, MPI_INT, MPI_SUM, comm);

    printf("Process %d: \n", my_rank);
    printf("x = ");
    for (i = 0; i < count; i++) {
        printf("%d ", x[i]);
    }
    printf("\n");

    printf("partial = ");
    for (i = 0; i < count; i++) {
        printf("%d ", partial[i]);
    }
    printf("\n");

    free(x);
    free(partial);

    MPI_Finalize();
    return 0;
}