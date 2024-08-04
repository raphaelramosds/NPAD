#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

const int limit = 1000;
const int size = 32768;

int main(void) {
    char buffer[size];
    
    int it = 0;
	
    int comm_sz;
    int my_rank;
    double total = 0.0;
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

    while (it < limit) {
        if (my_rank == 0) {
            for (int i=0; i < size - 1; i++) buffer[i] = 'A';

#if WTIME
            start = MPI_Wtime();
#else
            start_t = clock();
#endif

            // printf("[%d] Sending to process 1\n", it);
            // MPI_Send(&it, 1, MPI_INT, 1, 0, comm);
            MPI_Send(buffer, strlen(buffer) + 1, MPI_CHAR, 1, 0, comm);
            MPI_Recv(buffer, size, MPI_CHAR, 1, 0, comm, MPI_STATUS_IGNORE);

            // MPI_Recv(&it, 1, MPI_INT, 1, 1, comm, MPI_STATUS_IGNORE);
            // printf("[%d] Received from process 1\n", it);

#if WTIME
            end = MPI_Wtime();
            total += it > 0 ? end - start : 0.0;
#else
            end_t = clock();
            total += it > 0 ? (double)(end_t - start_t) / CLOCKS_PER_SEC : 0.0;
#endif

            //printf("%dth ping-pong = %.3es\n", it, total);

        } else {
            // MPI_Recv(&it, 1, MPI_INT, 0, 0, comm, MPI_STATUS_IGNORE);
            // printf("[%d] Received from process 0\n", it);

            MPI_Recv(buffer, size, MPI_CHAR, 0, 0, comm, MPI_STATUS_IGNORE);
            MPI_Send(buffer, strlen(buffer) + 1, MPI_CHAR, 0, 0, comm);

            // printf("[%d] Sending to process 0\n", it);
            // MPI_Send(&it, 1, MPI_INT, 0, 1, comm);
        }

        it++;
    }

	if (my_rank==0) printf("Average = %.3es", total/((double)it));

    MPI_Finalize();

    return 0;
}
