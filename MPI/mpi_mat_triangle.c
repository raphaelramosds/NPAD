#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    // Get the number of processes and check only 2 processes are used
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 2) {
        printf("This application is meant to be run with 2 processes.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Get my rank and do the corresponding job
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (my_rank == 0) {
        // Create the datatype
        MPI_Datatype triangle_type;
        int lengths[3] = {1, 2, 3};
        int displacements[3] = {0, 3, 6};
        MPI_Type_indexed(3, lengths, displacements, MPI_INT, &triangle_type);
        MPI_Type_commit(&triangle_type);

        // Send the message
        int buffer[3][3] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        MPI_Send(buffer, 1, triangle_type, 1, 0, MPI_COMM_WORLD);
    } else {
        int received[6];
        MPI_Recv(&received, 6, MPI_INT, 0, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        printf("MPI process %d received values:%d %d %d %d %d %d\n",
               my_rank, received[0], received[1], received[2], received[3],
               received[4], received[5]);
    }

    MPI_Finalize();

    return EXIT_SUCCESS;
}