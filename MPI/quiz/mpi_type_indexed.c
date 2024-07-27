#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

void Get_order(int* n_p, int my_rank, int comm_sz, MPI_Comm comm);

void Read_matrix(char prompt[], double* local_A, int n, int my_rank,
                 MPI_Comm comm);

int main(void) {
    double* local_A;
    int n;
    int my_rank, comm_sz;
    MPI_Comm comm;

    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);

    Get_order(&n, my_rank, comm_sz, comm);

    local_A = malloc(n * n * sizeof(double));

    Read_matrix("A", local_A, n, my_rank, comm);

    free(local_A);
    MPI_Finalize();
    return 0;
}

void Get_order(int* n_p, int my_rank, int comm_sz, MPI_Comm comm) {
    if (my_rank == 0) {
        printf("Enter the matrix order:\n");
        scanf("%d", n_p);
    }
    MPI_Bcast(n_p, 1, MPI_INT, 0, comm);
}

void Read_matrix(char prompt[], double* local_A, int n, int my_rank,
                 MPI_Comm comm) {
    int i, j;

    if (my_rank == 0) {

        int* block_lengths = malloc(n * sizeof(int));
        int* displacements = malloc(n * sizeof(int));

        /* Number of elements by block */
        for (i = n; i > 0; i--) block_lengths[n - i] = i;

        /* Distances between the first block and the others: (# lines)*i + i */
        for (i = 0; i < n; i++) displacements[i] = n * i + i;

        printf("Enter the matrix %s\n", prompt);
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) scanf("%lf", local_A + i*n + j);
        
        MPI_Datatype UTriangle_t;

        MPI_Type_indexed(n, block_lengths, displacements, MPI_DOUBLE,
                        &UTriangle_t);
        MPI_Type_commit(&UTriangle_t);

        MPI_Send(local_A, 1, UTriangle_t, 1, 0, comm);

        free(block_lengths);
        free(displacements);

    } else {
        int i;

        /* Number of elements on the upper triangle */
        int count = n * (n + 1) / 2;

        double* received = malloc(count * sizeof(double));

        MPI_Recv(received, count, MPI_DOUBLE, 0, 0, comm, MPI_STATUS_IGNORE);

        printf("Upper triangle matrix (order %d): \n", n);

        for (i = 0; i < count; i++) printf("%lf ", received[i]);

        free(received);
    }
}
