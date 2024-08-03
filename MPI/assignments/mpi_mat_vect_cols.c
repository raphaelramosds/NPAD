/* File:     mpi_mat_vect_time.c
 *
 * Purpose:  Implement parallel matrix-vector multiplication using
 *           one-dimensional arrays to store the vectors and the
 *           matrix.  Vectors use block distributions and the
 *           matrix is distributed by block rows.  This version
 *           generates a random matrix A and a random vector x.
 *           It prints out the run-time.
 *
 * Compile:  mpicc -g -Wall -o mpi_mat_vect_time mpi_mat_vect_time.c
 * Run:      mpiexec -n <number of processes> ./mpi_mat_vect_time
 *
 * Input:    Dimensions of the matrix (m = number of rows, n
 *              = number of columns)
 * Output:   Elapsed time for execution of the multiplication
 *
 * Notes:
 *    1. Number of processes should evenly divide both m and n
 *    2. Define DEBUG for verbose output, including the product
 *       vector y
 *
 * IPP:  Section 3.6.2 (pp. 122 and ff.)
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

void Check_for_error(int local_ok, char fname[], char message[], MPI_Comm comm);
void Get_dims(int* m_p, int* local_m_p, int* n_p, int* local_n_p, int my_rank,
              int comm_sz, MPI_Comm comm);
void Allocate_arrays(double** local_A_pp, double** local_x_pp, double** y_pp,
                     double** local_y_pp, int m, int local_m, int n,
                     int local_n, MPI_Comm comm);
void Read_matrix(char prompt[], double local_A[], int m, int local_m, int n,
                 int local_n, int my_rank, MPI_Comm comm);
void Read_vector(char prompt[], double local_vec[], int n, int local_n,
                 int my_rank, MPI_Comm comm);
void Generate_matrix(double local_A[], int local_m, int n);
void Generate_vector(double local_x[], int local_n);
void Print_matrix(char title[], double local_A[], int m, int local_m, int n,
                  int local_n, int my_rank, MPI_Comm comm);
void Print_vector(char title[], double local_vec[], int n, int local_n,
                  int my_rank, MPI_Comm comm);
void Mat_vect_mult(double local_A[], double local_x[], double y[],
                   double local_y[], int m, int local_m, int n, int local_n,
                   MPI_Comm comm);

/*-------------------------------------------------------------------*/
int main(int argc, int** argv) {
    double* y;
    double* local_A;
    double* local_x;
    double* local_y;
    int m, local_m, n, local_n;
    int my_rank, comm_sz;
    MPI_Comm comm;
    double start, finish, loc_elapsed, elapsed;

    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);

    m = atoi(argv[1]); n = atoi(argv[2]);

    Get_dims(&m, &local_m, &n, &local_n, my_rank, comm_sz, comm);
    Allocate_arrays(&local_A, &local_x, &y, &local_y, m, local_m, n, local_n,
                    comm);
    srandom(my_rank);
    Read_matrix("A", local_A, m, local_m, n, local_n, my_rank, comm);
    // Generate_matrix(local_A, local_m, n);

#ifdef DEBUG
    // Print_matrix("A", local_A, m, local_m, n, local_n, my_rank, comm);
#endif

#ifdef DEBUG
    int local_i;
    printf("local_a(%d) =", my_rank);
    for (local_i = 0; local_i < local_n * m; local_i++) {
        printf("%f ", local_A[local_i]);
    }
    printf("\n");
#endif
    // Read_vector("x", local_x, n, local_n, my_rank, comm);
    Generate_vector(local_x, local_n);
#ifdef DEBUG
    Print_vector("x", local_x, n, local_n, my_rank, comm);
#endif

    MPI_Barrier(comm);
    start = MPI_Wtime();
    Mat_vect_mult(local_A, local_x, y, local_y, m, local_m, n, local_n, comm);
    finish = MPI_Wtime();
    loc_elapsed = finish - start;
    MPI_Reduce(&loc_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

#ifdef DEBUG
    // Print_vector("y", local_y, m, local_m, my_rank, comm);
    if (my_rank == 0) {
        int local_i;
        printf("y=");
        for (local_i = 0; local_i < m; local_i++) printf("%f ", y[local_i]);
        printf("\n");
    }
#endif

    if (my_rank == 0) 
        printf("%d,%d,%f\n", comm_sz, m, elapsed);

    free(local_A);
    free(local_x);
    free(local_y);
    MPI_Finalize();
    return 0;
} /* main */

/*-------------------------------------------------------------------*/
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
} /* Check_for_error */

/*-------------------------------------------------------------------*/
void Get_dims(int* m_p /* out */, int* local_m_p /* out */, int* n_p /* out */,
              int* local_n_p /* out */, int my_rank /* in  */,
              int comm_sz /* in  */, MPI_Comm comm /* in  */) {
    int local_ok = 1;

    // if (my_rank == 0) {
    //     printf("Enter the number of rows\n");
    //     scanf("%d", m_p);
    //     printf("Enter the number of columns\n");
    //     scanf("%d", n_p);
    // }
    MPI_Bcast(m_p, 1, MPI_INT, 0, comm);
    MPI_Bcast(n_p, 1, MPI_INT, 0, comm);
    if (*m_p <= 0 || *n_p <= 0 || *m_p % comm_sz != 0 || *n_p % comm_sz != 0)
        local_ok = 0;
    Check_for_error(local_ok, "Get_dims",
                    "m and n must be positive and evenly divisible by comm_sz",
                    comm);

    // *local_m_p = *m_p/comm_sz;
    // *local_n_p = *n_p/comm_sz;

    /* Processes should access all lines of the matrix */
    *local_m_p = *m_p;

    /* Processes should have their own columns */
    *local_n_p = *n_p / comm_sz;

} /* Get_dims */

/*-------------------------------------------------------------------*/
void Allocate_arrays(double** local_A_pp /* out */,
                     double** local_x_pp /* out */, double** y_pp,
                     double** local_y_pp /* out */, int m /* in */,
                     int local_m /* in  */, int n /* in  */,
                     int local_n /* in  */, MPI_Comm comm /* in  */) {
    int local_ok = 1;

    /* (m, local_n) because processes have all lines of A and and part of its
     * columns */
    *local_A_pp = malloc(m * local_n * sizeof(double));

    /* (1, local_n) as this program join all components of x with allgather */
    *local_x_pp = malloc(local_n * sizeof(double));

    /* (1, m) since processes must have a copy of all components of y */
    *y_pp = malloc(m * sizeof(double));
    *local_y_pp = malloc(m * sizeof(double));

    if (*local_A_pp == NULL || local_x_pp == NULL || local_y_pp == NULL)
        local_ok = 0;
    Check_for_error(local_ok, "Allocate_arrays", "Can't allocate local arrays",
                    comm);
} /* Allocate_arrays */

/*-------------------------------------------------------------------*/
void Read_matrix(char prompt[] /* in  */, double local_A[] /* out */,
                 int m /* in  */, int local_m /* in  */, int n /* in  */,
                 int local_n /* in  */, int my_rank /* in  */,
                 MPI_Comm comm /* in  */) {
    double* A = NULL;
    int local_ok = 1;
    int i, j;

    MPI_Datatype col, col_type;
    MPI_Type_vector(m, 1, n, MPI_DOUBLE, &col);
    MPI_Type_commit(&col);

    /* Resize type vector in order to MPI Scatter accesses it contigously */
    MPI_Type_create_resized(col, 0, sizeof(double), &col_type);
    MPI_Type_commit(&col_type);

    if (my_rank == 0) {
        A = malloc(m * n * sizeof(double));
        if (A == NULL) local_ok = 0;
        Check_for_error(local_ok, "Read_matrix",
                        "Can't allocate temporary matrix", comm);
        // printf("Enter the matrix %s\n", prompt);
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
                A[i * n + j] = ((double)random()) / ((double)RAND_MAX);
        //    scanf("%lf", &A[i*n+j]);

        // MPI_Scatter(A, local_m*n, MPI_DOUBLE,
        //       local_A, local_m*n, MPI_DOUBLE, 0, comm);

        /* Send as a column type and receive as an array */
        MPI_Scatter(A, local_n, col_type, local_A, local_n * m, MPI_DOUBLE, 0,
                    comm);

        free(A);
    } else {
        Check_for_error(local_ok, "Read_matrix",
                        "Can't allocate temporary matrix", comm);

        // MPI_Scatter(A, local_m*n, MPI_DOUBLE,
        //       local_A, local_m*n, MPI_DOUBLE, 0, comm);
        MPI_Scatter(NULL, NULL, NULL, local_A, local_n * m, MPI_DOUBLE, 0,
                    comm);
    }
} /* Read_matrix */

/*-------------------------------------------------------------------*/
void Read_vector(char prompt[] /* in  */, double local_vec[] /* out */,
                 int n /* in  */, int local_n /* in  */, int my_rank /* in  */,
                 MPI_Comm comm /* in  */) {
    double* vec = NULL;
    int i, local_ok = 1;

    if (my_rank == 0) {
        vec = malloc(n * sizeof(double));
        if (vec == NULL) local_ok = 0;
        Check_for_error(local_ok, "Read_vector",
                        "Can't allocate temporary vector", comm);
        printf("Enter the vector %s\n", prompt);
        for (i = 0; i < n; i++) scanf("%lf", &vec[i]);
        MPI_Scatter(vec, local_n, MPI_DOUBLE, local_vec, local_n, MPI_DOUBLE, 0,
                    comm);
        free(vec);
    } else {
        Check_for_error(local_ok, "Read_vector",
                        "Can't allocate temporary vector", comm);
        MPI_Scatter(vec, local_n, MPI_DOUBLE, local_vec, local_n, MPI_DOUBLE, 0,
                    comm);
    }
} /* Read_vector */

/*-------------------------------------------------------------------*/
void Generate_matrix(double local_A[] /* out */, int local_m /* in  */,
                     int n /* in  */) {
    int i, j;

    for (i = 0; i < local_m; i++)
        for (j = 0; j < n; j++)
            local_A[i * n + j] = ((double)random()) / ((double)RAND_MAX);
} /* Generate_matrix */

/*-------------------------------------------------------------------*/
void Generate_vector(double local_x[] /* out */, int local_n /* in  */) {
    int i;

    for (i = 0; i < local_n; i++)
        local_x[i] = ((double)random()) / ((double)RAND_MAX);
} /* Generate_vector */

/*-------------------------------------------------------------------*/
void Print_matrix(char title[] /* in */, double local_A[] /* in */,
                  int m /* in */, int local_m /* in */, int n /* in */,
                  int local_n /* in */, int my_rank /* in */,
                  MPI_Comm comm /* in */) {
    double* A = NULL;
    int i, j, local_ok = 1;

    if (my_rank == 0) {
        A = malloc(m * n * sizeof(double));
        if (A == NULL) local_ok = 0;
        Check_for_error(local_ok, "Print_matrix",
                        "Can't allocate temporary matrix", comm);
        MPI_Gather(local_A, local_n * m, MPI_DOUBLE, A, local_n * m, MPI_DOUBLE,
                   0, comm);
        printf("\nThe matrix %s\n", title);
        for (i = 0; i < m; i++) {
            for (j = 0; j < n; j++) printf("%f ", A[i * n + j]);
            printf("\n");
        }
        printf("\n");
        free(A);
    } else {
        Check_for_error(local_ok, "Print_matrix",
                        "Can't allocate temporary matrix", comm);
        MPI_Gather(local_A, local_n * m, MPI_DOUBLE, A, local_n * m, MPI_DOUBLE,
                   0, comm);
    }
} /* Print_matrix */

/*-------------------------------------------------------------------*/
void Print_vector(char title[] /* in */, double local_vec[] /* in */,
                  int n /* in */, int local_n /* in */, int my_rank /* in */,
                  MPI_Comm comm /* in */) {
    double* vec = NULL;
    int i, local_ok = 1;

    if (my_rank == 0) {
        vec = malloc(n * sizeof(double));
        if (vec == NULL) local_ok = 0;
        Check_for_error(local_ok, "Print_vector",
                        "Can't allocate temporary vector", comm);
        MPI_Gather(local_vec, local_n, MPI_DOUBLE, vec, local_n, MPI_DOUBLE, 0,
                   comm);
        printf("\nThe vector %s\n", title);
        for (i = 0; i < n; i++) printf("%f ", vec[i]);
        printf("\n");
        free(vec);
    } else {
        Check_for_error(local_ok, "Print_vector",
                        "Can't allocate temporary vector", comm);
        MPI_Gather(local_vec, local_n, MPI_DOUBLE, vec, local_n, MPI_DOUBLE, 0,
                   comm);
    }
} /* Print_vector */

/*-------------------------------------------------------------------*/
void Mat_vect_mult(double local_A[] /* in  */, double local_x[] /* in  */,
                   double y[], double local_y[] /* out */, int m /* in  */,
                   int local_m /* in  */, int n /* in  */,
                   int local_n /* in  */, MPI_Comm comm /* in  */) {
    //  double* x;
    int i, local_j;
    int local_ok = 1;

    // x = malloc(n*sizeof(double));
    // if (x == NULL) local_ok = 0;
    Check_for_error(local_ok, "Mat_vect_mult",
                    "Can't allocate temporary vector", comm);
    // MPI_Allgather(local_x, local_n, MPI_DOUBLE,
    //       x, local_n, MPI_DOUBLE, comm);

    for (i = 0; i < m; i++) local_y[i] = 0.0;

    /* The product (m elements of a column) and (1 component of x) fills (m*
     * coordinates of y) */
    for (local_j = 0; local_j < local_n; local_j++) {
        for (i = 0; i < m; i++)
            local_y[i] += local_A[local_j * m + i] * local_x[local_j];
    }

    /* Reduce all local_y into y */
    MPI_Allreduce(local_y, y, m, MPI_DOUBLE, MPI_SUM, comm);

} /* Mat_vect_mult */
