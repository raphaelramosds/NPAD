#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

const int LIMIT = 500;

int main(void) {
    int response = 0;
    int comm_sz;
    int my_rank;
    clock_t start_t, end_t;
    double total_t;
    MPI_Comm comm = MPI_COMM_WORLD;

    /* Start up MPI */
    MPI_Init(NULL, NULL);

    /* Get the number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    /* Program should run with only two processes */
    if (comm_sz != 2) {
        fprintf(stderr, "Run this program with two processes");
        MPI_Abort(comm, 1);
    }

    /* Get my rank among all the processes */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    while (response < LIMIT) {
        if (my_rank == 0) {
            start_t = clock();

            MPI_Send(&response, 1, MPI_INT, 1, 0, comm);
            MPI_Recv(&response, 1, MPI_INT, 1, 1, comm, MPI_STATUS_IGNORE);

            end_t = clock();
            total_t = (double) (end_t - start_t) / CLOCKS_PER_SEC;

            if (total_t > 0) {
                printf("Elapsed (%d ping-pong): %fs\n", response, total_t);
                MPI_Abort(comm, 1);
            }

        } else {
            MPI_Recv(&response, 1, MPI_INT, 0, 0, comm, MPI_STATUS_IGNORE);
            MPI_Send(&response, 1, MPI_INT, 0, 1, comm);
        }
        response++;
    }
    
    MPI_Finalize();

    return 0;
}
