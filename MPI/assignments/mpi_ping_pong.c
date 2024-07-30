#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

const int LIMIT = 1000;

const int MAX_STRING = 100;

int main(void) {
    char msg[MAX_STRING];
    int it = 0;
    int comm_sz;
    int my_rank;
    double total;
#ifdef WTIME
    double start, end;
#else
    clock_t start_t, end_t;
#endif

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

    sprintf(msg, "Greeting from process %d", my_rank);

    while (it < LIMIT) {
        if (my_rank == 0) {
#if WTIME
            start = MPI_Wtime();
#else
            start_t = clock();
#endif

            // printf("[%d] Sending to process 1\n", it);
            // MPI_Send(&it, 1, MPI_INT, 1, 0, comm);
            MPI_Send(msg, strlen(msg) + 1, MPI_CHAR, 1, 0, comm);
            MPI_Recv(msg, MAX_STRING, MPI_CHAR, 1, 0, comm, MPI_STATUS_IGNORE);

            // MPI_Recv(&it, 1, MPI_INT, 1, 1, comm, MPI_STATUS_IGNORE);
            // printf("[%d] Received from process 1\n", it);

#if WTIME
            end = MPI_Wtime();
            total = end - start;
#else
            end_t = clock();
            total = (double)(end_t - start_t) / CLOCKS_PER_SEC;
#endif

            printf("%dth ping-pong = %.3es\n", it, total);

        } else {
            // MPI_Recv(&it, 1, MPI_INT, 0, 0, comm, MPI_STATUS_IGNORE);
            // printf("[%d] Received from process 0\n", it);

            MPI_Recv(msg, MAX_STRING, MPI_CHAR, 0, 0, comm, MPI_STATUS_IGNORE);
            MPI_Send(msg, strlen(msg) + 1, MPI_CHAR, 0, 0, comm);

            // printf("[%d] Sending to process 0\n", it);
            // MPI_Send(&it, 1, MPI_INT, 0, 1, comm);
        }

        it++;
    }

    MPI_Finalize();

    return 0;
}
