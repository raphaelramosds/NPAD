#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

void Check_for_error(int local_ok, char fname[], char message[], MPI_Comm comm);
void Read_scalar(double* s_p, double* local_s_p, int my_rank, int comm_sz,
                 MPI_Comm comm);
void Read_n(int* n_p, int* local_n_p, int* r, int* q, int my_rank, int comm_sz, MPI_Comm comm);
void Read_vector(double local_a[], int local_n, int n, char vec_name[], int r, int q,
                 int my_rank, int comm_sz, MPI_Comm comm);
void Print_vector(double local_b[], int local_n, int n, char title[], int r, int q,
                  int my_rank, int comm_sz, MPI_Comm comm);
void Scale_vector(double local_x[], double local_z[], double local_s,
                  int local_n);
void Dot_product(double local_x[], double local_y[], int local_n,
                 double local_r, double total_r, int my_rank, MPI_Comm comm);

int main(void) {
    double total_r, local_r, s, local_s;
    int n, local_n, comm_sz, my_rank, r, q;
    double *local_x, *local_y, *local_z;
    MPI_Comm comm;

    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);

    Read_scalar(&s, &local_s, my_rank, comm_sz, comm);
    Read_n(&n, &local_n, &r, &q, my_rank, comm_sz, comm);

#ifdef DEBUG
    printf("Proc %d > n = %d, local_n = %d\n", my_rank, n, local_n);
#endif
    Allocate_vectors(&local_x, &local_y, &local_z, local_n, comm);

    Read_vector(local_x, local_n, n, "x", r, q, my_rank, comm_sz, comm);
    Print_vector(local_x, local_n, n, "x is", r, q, my_rank, comm_sz, comm);
    Read_vector(local_y, local_n, n, "y", r, q, my_rank, comm_sz, comm);
    Print_vector(local_y, local_n, n, "y is", r, q, my_rank, comm_sz, comm);

    Scale_vector(local_x, local_z, local_s, local_n);
    Print_vector(local_z, local_n, n, "z is", r, q, my_rank, comm_sz, comm);

    Dot_product(local_x, local_y, local_n, local_r, total_r, my_rank, comm);

    free(local_x);
    free(local_y);
    free(local_z);

    MPI_Finalize();

    return 0;
}

/*-------------------------------------------------------------------
 * Function:  Check_for_error
 * Purpose:   Check whether any process has found an error.  If so,
 *            print message and terminate all processes.  Otherwise,
 *            continue execution.
 * In args:   local_ok:  1 if calling process has found an error, 0
 *               otherwise
 *            fname:     name of function calling Check_for_error
 *            message:   message to print if there's an error
 *            comm:      communicator containing processes calling
 *                       Check_for_error:  should be MPI_COMM_WORLD.
 *
 * Note:
 *    The communicator containing the processes calling Check_for_error
 *    should be MPI_COMM_WORLD.
 */
void Check_for_error(int local_ok /* in */, char fname[] /* in */,
                     char message[] /* in */, MPI_Comm comm /* in */) {
    int ok;

    MPI_Allreduce(&local_ok, &ok, 1, MPI_INT, MPI_MIN, comm);
    if (ok == 0) {
        int my_rank;
        MPI_Comm_rank(comm, &my_rank);
        if (my_rank == 0) {
            fprintf(stderr, "Proc %d > In %s, %s\n", my_rank, fname, message);
            fflush(stderr);
        }
        MPI_Finalize();
        exit(-1);
    }
}

/*-------------------------------------------------------------------
 * Function:  Read_n
 * Purpose:   Get the order of the vectors from stdin on proc 0 and
 *            broadcast to other processes.
 * In args:   my_rank:    process rank in communicator
 *            comm_sz:    number of processes in communicator
 *            comm:       communicator containing all the processes
 *                        calling Read_n
 * Out args:  n_p:        global value of n
 *            local_n_p:  local value of n = n/comm_sz
 *
 * Errors:    n should be positive and evenly divisible by comm_sz
 */
void Read_n(int* n_p /* out */, int* local_n_p /* out */, int* r, int* q, int my_rank /* in  */,
            int comm_sz /* in  */, MPI_Comm comm /* in  */) {

    int local_ok = 1;
    char* fname = "Read_n";

    if (my_rank == 0) {
        printf("What's the order of the vectors?\n");
        scanf("%d", n_p);
    }
    MPI_Bcast(n_p, 1, MPI_INT, 0, comm);
    if (*n_p <= 0 || *n_p % comm_sz != 0) local_ok = 0;
    // Check_for_error(local_ok, fname,
    //                 "n should be > 0 and evenly divisible by comm_sz", comm);

    // *local_n_p = *n_p / comm_sz;
    *r = *n_p % comm_sz;
    *q = *n_p / comm_sz;
    *local_n_p = my_rank < *r ? *q + 1 : *q;
}

/*-------------------------------------------------------------------
 * Function:  Read_scalar
 * Purpose:   Get the scalar from stdin on proc 0 and broadcast to
 *            other processes.
 * In args:   my_rank:    process rank in communicator
 *            comm_sz:    number of processes in communicator
 *            comm:       communicator containing all the processes
 *                        calling Read_n
 * Out args:  s_p:        global value of scalar
 *            local_s_p:  local value of scalar
 *
 * Errors:
 */
void Read_scalar(double* s_p, double* local_s_p, int my_rank, int comm_sz,
                 MPI_Comm comm) {
    int local_ok = 1;
    char* fname = "Read_scalar";

    if (my_rank == 0) {
        printf("What's the scalar?\n");
        scanf("%lf", s_p);
    }
    MPI_Bcast(s_p, 1, MPI_DOUBLE, 0, comm);
    *local_s_p = *s_p;
}

/*-------------------------------------------------------------------
 * Function:   Read_vector
 * Purpose:    Read a vector from stdin on process 0 and distribute
 *             among the processes using a block distribution.
 * In args:    local_n:  size of local vectors
 *             n:        size of global vector
 *             vec_name: name of vector being read (e.g., "x")
 *             my_rank:  calling process' rank in comm
 *             comm:     communicator containing calling processes
 * Out arg:    local_a:  local vector read
 *
 * Errors:     if the malloc on process 0 for temporary storage
 *             fails the program terminates
 *
 * Note:
 *    This function assumes a block distribution and the order
 *   of the vector evenly divisible by comm_sz.
 */
void Read_vector(double local_a[] /* out */, int local_n /* in  */,
                 int n /* in  */, char vec_name[] /* in  */, int r, int q,
                 int my_rank /* in  */, int comm_sz /* in */,  MPI_Comm comm /* in  */) {
    double* a = NULL;
    int i;
    int local_ok = 1;
    char* fname = "Read_vector";

    int* sendcounts = malloc(comm_sz*sizeof(int));

    for (i = 0; i < comm_sz; i++)
        sendcounts[i] = i < r ? q + 1 : q;

    if (my_rank == 0) {
        int curr_idx = 0;
        int* displacements = malloc(comm_sz*sizeof(int));

        for (i = 0; i < comm_sz; i++) {
            displacements[i] = curr_idx;
            curr_idx += sendcounts[i];
        }

        a = malloc(n * sizeof(double));
        if (a == NULL) local_ok = 0;
        Check_for_error(local_ok, fname, "Can't allocate temporary vector",
                        comm);
        printf("Enter the vector %s\n", vec_name);
        for (i = 0; i < n; i++) scanf("%lf", &a[i]);
        // MPI_Scatter(a, local_n, MPI_DOUBLE, local_a, local_n, MPI_DOUBLE, 0,
        //             comm);
        MPI_Scatterv(a, sendcounts, displacements, MPI_DOUBLE, local_a, sendcounts[0] , MPI_DOUBLE, 0, comm);

        free(a);
        free(displacements);
    } else {
        Check_for_error(local_ok, fname, "Can't allocate temporary vector",
                        comm);
        // MPI_Scatter(a, local_n, MPI_DOUBLE, local_a, local_n, MPI_DOUBLE, 0,
        //             comm);
        MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, local_a, sendcounts[my_rank], MPI_DOUBLE, 0, comm);
    }
    free(sendcounts);
} /* Read_vector */

/*-------------------------------------------------------------------
 * Function:  Print_vector
 * Purpose:   Print a vector that has a block distribution to stdout
 * In args:   local_b:  local storage for vector to be printed
 *            local_n:  order of local vectors
 *            n:        order of global vector (local_n*comm_sz)
 *            title:    title to precede print out
 *            comm:     communicator containing processes calling
 *                      Print_vector
 *
 * Error:     if process 0 can't allocate temporary storage for
 *            the full vector, the program terminates.
 *
 * Note:
 *    Assumes order of vector is evenly divisible by the number of
 *    processes
 */
void Print_vector(double local_b[] /* in */, int local_n /* in */,
                  int n /* in */, char title[] /* in */, int r, int q, int my_rank /* in */,
                  int comm_sz, MPI_Comm comm /* in */) {
    double* b = NULL;
    int i;
    int local_ok = 1;
    char* fname = "Print_vector";

    int* sendcounts = malloc(comm_sz*sizeof(int));

    for (i = 0; i < comm_sz; i++)
        sendcounts[i] = i < r ? q + 1 : q;

    if (my_rank == 0) {
        int curr_idx = 0;
        int* displacements = malloc(comm_sz*sizeof(int));

        for (i = 0; i < comm_sz; i++) {
            displacements[i] = curr_idx;
            curr_idx += sendcounts[i];
        }

        b = malloc(n * sizeof(double));
        if (b == NULL) local_ok = 0;
        Check_for_error(local_ok, fname, "Can't allocate temporary vector",
                        comm);
        // MPI_Gather(local_b, local_n, MPI_DOUBLE, b, local_n, MPI_DOUBLE, 0,
        //            comm);
        MPI_Gatherv(local_b, sendcounts[0], MPI_DOUBLE, b, sendcounts, displacements, MPI_DOUBLE, 0, comm);
        printf("%s\n", title);
        for (i = 0; i < n; i++) printf("%f ", b[i]);
        printf("\n");
        free(b);
        free(displacements);
    } else {
        Check_for_error(local_ok, fname, "Can't allocate temporary vector",
                        comm);
        // MPI_Gather(local_b, local_n, MPI_DOUBLE, b, local_n, MPI_DOUBLE, 0,
        //            comm);
        MPI_Gatherv(local_b, sendcounts[my_rank], MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, comm);
    }
    free(sendcounts);
}

/*-------------------------------------------------------------------
 * Function:  Scale_vector
 * Purpose:   Scale a vector by a real factor
 * In args:   local_x:   local storage for vector to be scaled
 *            local_z:   local scaled vector
 *            local_s:   factor
 *            local_n:   order of local vector
 *
 * Error:
 *
 * Note:
 */
void Scale_vector(double local_x[], double local_z[], double local_s,
                  int local_n) {
    int local_i;

    for (local_i = 0; local_i < local_n; local_i++)
        local_z[local_i] = local_x[local_i] * local_s;
}

/*-------------------------------------------------------------------
 * Function:  Dot_product
 * Purpose:   Performs the dot product between two vectors
 *
 * Error:
 *
 * Note:
 */
void Dot_product(double local_x[], double local_y[], int local_n,
                 double local_r, double total_r, int my_rank, MPI_Comm comm) {
    int local_i;

    for (local_i = 0; local_i < local_n; local_i++)
        local_r += local_x[local_i] * local_y[local_i];

    MPI_Reduce(&local_r, &total_r, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (my_rank == 0) {
        printf("Dot product between x and y is: %f", total_r);
    }
}

/*-------------------------------------------------------------------
 * Function:  Allocate_vectors
 * Purpose:   Allocate storage for x, y, and z
 * In args:   local_n:  the size of the local vectors
 *            comm:     the communicator containing the calling processes
 * Out args:  local_x_pp, local_y_pp, local_z_pp:  pointers to memory
 *               blocks to be allocated for local vectors
 *
 * Errors:    One or more of the calls to malloc fails
 */
void Allocate_vectors(double** local_x_pp /* out */,
                      double** local_y_pp /* out */,
                      double** local_z_pp /* out */, int local_n /* in  */,
                      MPI_Comm comm /* in  */) {
    int local_ok = 1;
    char* fname = "Allocate_vectors";

    *local_x_pp = malloc(local_n * sizeof(double));
    *local_y_pp = malloc(local_n * sizeof(double));
    *local_z_pp = malloc(local_n * sizeof(double));

    if (*local_x_pp == NULL || *local_y_pp == NULL || *local_z_pp == NULL)
        local_ok = 0;
    Check_for_error(local_ok, fname, "Can't allocate local vector(s)", comm);
}