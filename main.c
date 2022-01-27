/* 
 * PHYM0004 Assignment 1 / Jay Malhotra / 07/10/2021
 * 
 * main.c features the data collection for the report. It generates a 100x100
 * matrix full of random numbers, and then checks determinant and LU determinant
 * operations on each leading matrix from 2x2 up to 100x100.
 * 
 * Data is gathered by piping the output of the program to a text file; the
 * prints have been designed to be compatible with numpy.loadtxt.
 */

#include "matrix_core.h"
#include "determinant.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

// Comment this definition out if you want to disable brute-force calculation to
// see how LU decomposition performs at larger sizes without having to wait
// several years for brute-force to catch up
// It's a bit of a hack but writing user interfaces in C is tricky
//#define USEBRUTE

long int timespec_diff(struct timespec t1, struct timespec t2) {
    // timeseries.tv_nsec is the nanoseconds within the current second. So
    // (t1.tv_nsec - t2.tv_nsec) will give you a negative number if t1 and t2
    // are in different seconds. This method intends to provide the accurate
    // nanosecond difference between two timespecs.
    return (t1.tv_sec - t2.tv_sec) * 1e9 + (t1.tv_nsec - t2.tv_nsec);
}

int main(void) {
    // Set the size of the master matrix -- i.e. the largest size to do
    // the benchmarking for
    int size = 1000;
    // The numbers in the matrix range between min -> max
    double min = -10;
    double max = 10;

    // Generate random matrix
    double* matrix_data = malloc((size*size) * sizeof(double));
    if (matrix_data == NULL) {
        printf("Could not allocate memory for %d required values", size*size);
        return 1;
    }

    srand(time(0));
    for (int i = 0; i < size*size; i++) {
        *(matrix_data + i) = min + rand() / (RAND_MAX / (max - min + 1) + 1);
    }

    struct DoubleMatrix* huge_matrix = matrix_factory(size, size, matrix_data);
    free(matrix_data);

    struct timespec start_time;
    struct timespec brute_time;
    struct timespec lu_time;

    printf("# Columns:\n# Matrix size / LU time (ns) / Brute-force time (ns) (0=disabled)\n");

    // Iterate through leading matrices
    for (int i = 2; i <= size; i++) {
        struct DoubleMatrix* matrix = matrix_leading(huge_matrix, i);
        printf("%d\t", i);

        /* I have elected not to perform checks on leading matrices to decide
         * whether or not to do row swaps because the use of random decimal
         * numbers makes it incredibly unlikely that any leading matrix will
         * have a determinant of zero. If I was using random integers, this may
         * be a concern, but so far I have not encountered a non-decomposable
         * matrix. */

        clock_gettime(CLOCK_REALTIME, &start_time);
        
        double lu = calculate_determinant_lu(matrix);
        clock_gettime(CLOCK_REALTIME, &lu_time);
        long int lu_nsec = timespec_diff(lu_time, start_time);
        printf("%ld\t", lu_nsec);

        // Print 0 for brute time unless enabled
        long int brute_nsec = 0;

        #ifdef USEBRUTE

        double brute = calculate_determinant(matrix);
        clock_gettime(CLOCK_REALTIME, &brute_time);
        brute_nsec = timespec_diff(brute_time, lu_time);
        // Check that the two methods give the same result
        assert(abs(1 - brute/lu) < 0.01);

        #endif

        printf("%ld\t\n", brute_nsec);
        matrix_free(matrix);
    }

    matrix_free(huge_matrix);

    return 0;
}