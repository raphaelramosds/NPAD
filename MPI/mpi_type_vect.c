/* File:     mpi_type_vect.c
 * Purpose:  Illustrates how to create a vector MPI datatype.
 *
 * Compile:  mpicc -g -Wall -o mpi_type_vect.c
 * Run:      mpiexec -n <number of processes> ./mpi_type_vect
 *
 * Algorithm:
 *    1.  Two MPI processes wil exchange a message made of three integers
 *    2.  On the sender, that message is in fact the middle column of an array
 * it holds
 *    3.  This column will be represented by an MPI vector
 *
 * Note:  This program is meant to be run with 2 processes
 *
 * IPP:   Rookie HPC - MPI_Type_vector
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    int comm_sz, my_rank;

    /* Let the system do what it needs to start up MPI */
    MPI_Init(&argc, &argv);

    /* Get my process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Find out how many processes are being used */
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    if (comm_sz != 2) {
        printf("This application is meant to be run with 2 processes.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (my_rank == 0) {
        /* Create the datatype */
        MPI_Datatype column_type;
        MPI_Type_vector(3, 1, 3, MPI_INT, &column_type);
        MPI_Type_commit(&column_type);

        /* Send the message */
        int buffer[3][3] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        MPI_Request request;
        printf("MPI process %d sends values %d, %d and %d.\n", my_rank,
               buffer[0][1], buffer[1][1], buffer[2][1]);
        MPI_Send(&buffer[0][1], 1, column_type, 1, 0, MPI_COMM_WORLD);
    } else {
        /* Receive the message */
        int received[3];
        MPI_Recv(&received, 3, MPI_INT, 0, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        printf("MPI process %d received values: %d, %d and %d.\n", my_rank,
               received[0], received[1], received[2]);
    }

    MPI_Finalize();

    return EXIT_SUCCESS;
}
