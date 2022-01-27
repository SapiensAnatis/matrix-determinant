/* 
 * PHYM0004 Assignment 1 / Jay Malhotra / 18/10/2021
 *
 * determinant.c implements the methods defined in determinant.h. For more info
 * and explanations about what each method does, please see the aforementioned
 * header file.
 */

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "matrix_core.h"
#include "determinant.h"

double calculate_determinant(struct DoubleMatrix* matrix_ptr) {
    // Check square
    assert(matrix_ptr->numRows == matrix_ptr->numCols);

    // Backstop for recursive method. For 1x1 matrix return inner value
    if (matrix_ptr->numVals == 1) {
        return *(matrix_ptr->values);
    }

    double result = 0;

    // Determinant of A for NxN matrix:
    // a_11*C_11 + a_12*C_12 + a_13*C_13 + ... + a_1m*C_1m
    // where C_1m is the cofactor of element a_1m.
    for (int col = 0; col < matrix_ptr->numCols; col++) {
        result += *matrix_get_element(matrix_ptr, 0, col) * calculate_cofactor(matrix_ptr, 0, col);
    }

    return result;
}

double calculate_cofactor(struct DoubleMatrix* matrix_ptr, int n, int m) {
    // Assert that these are valid indices
    assert(n >= 0 && n <= matrix_ptr->numRows);
    assert(m >= 0 && m <= matrix_ptr->numCols);

    double result = 0;

    struct DoubleMatrix* minor_matrix = matrix_minor(matrix_ptr, n, m);

    // If n+m is even, sign=1. Else sign=-1. Equivalent to (-1)^(n+m)
    int sign = 0;
    sign = ((n+m) % 2 == 0) ? 1 : -1;

    result = sign * calculate_determinant(minor_matrix);
    matrix_free(minor_matrix);

    return result;
}

void calculate_lu_matrices(struct DoubleMatrix* matrix_ptr, 
                           struct DoubleMatrix* L_buffer, struct DoubleMatrix* U_buffer) {
    // This is a bit of a mess, but hopefully shouldn't be too hard to follow
    // when comparing against slide 15 of Lecture 5

    // Assert square
    assert(matrix_ptr->numRows == matrix_ptr->numCols);

    // Doesn't matter which, already asserted square
    int n = matrix_ptr->numRows;

    // Create empty matrices. This method works by setting values rather than
    // outputting to array and making matrices later -- easier to preserve
    // structure that way
    struct DoubleMatrix* L = matrix_factory(n, n, NULL);
    struct DoubleMatrix* U = matrix_factory(n, n, NULL);

    // Set diagonals of L to 1
    for (int diag = 0; diag < n; diag++) {
        matrix_set_element(L, diag, diag, 1);
    }

    // Get preliminary U values by copying top row of input matrix
    // 'Step 1' on lecture slides
    for (int col = 0; col < n; col++) {
        double m_element = *matrix_get_element(matrix_ptr, 0, col);
        matrix_set_element(U, 0, col, m_element);
    }

    // Get first column of L values by dividing U elements
    // 'Step 2' on lecture slides
    // Note that the loop starts at 1 instead of 0 due to setting diagonal Ls
    // to 1, so L_00 for instance is not needed
    double a_11 = *matrix_get_element(matrix_ptr, 0, 0);
    for (int row = 1; row < n; row++) {
        // a_21/a_11, a_31/a_11, etc.
        double l_value = *matrix_get_element(matrix_ptr, row, 0) / a_11;
        matrix_set_element(L, row, 0, l_value);
    }

    // Main iteration of method
    for (int i = 1; i < n; i++) {
        // Calculate u_ii
        double u_ii = *matrix_get_element(matrix_ptr, i, i);
        for (int t = 0; t < i; t++) {
            double l_it = *matrix_get_element(L, i, t);
            double u_ti = *matrix_get_element(U, t, i);
            u_ii -= (l_it * u_ti);
        }
        matrix_set_element(U, i, i, u_ii);
        
        for (int j = i+1; j < n; j++) {
            // Calculate u_ij
            double u_ij = *matrix_get_element(matrix_ptr, i, j);
            for (int t = 0; t < i; t++) {
                double l_it = *matrix_get_element(L, i, t);
                double u_tj = *matrix_get_element(U, t, j);
                u_ij -= l_it * u_tj;
            }
            matrix_set_element(U, i, j, u_ij);

            // Calculate l_ji
            double l_ji = *matrix_get_element(matrix_ptr, j, i);
            for (int t = 0; t < i; t++) {
                double l_jt = *matrix_get_element(L, j, t);
                double u_ti = *matrix_get_element(U, t, i);
                l_ji -= l_jt * u_ti;
            }
            // Divide by u_ii
            l_ji /= *matrix_get_element(U, i, i);
            matrix_set_element(L, j, i, l_ji);
        }
    }
    // Finally calculate u_nn. Note that n is the length and not the final
    // index, so we do have to do n-1 here.
    double u_nn = *matrix_get_element(matrix_ptr, n-1, n-1);
    for (int t = 0; t < n-1; t++) {
        double l_jt = *matrix_get_element(L, n-1, t);
        double u_ti = *matrix_get_element(U, t, n-1);
        u_nn -= l_jt * u_ti;
    }
    matrix_set_element(U, n-1, n-1, u_nn);

    // Finished! Copy over result into buffer matrices
    for (int i = 0; i < L->numVals; i++) {
        *(L_buffer->values + i) = *(L->values + i);
        *(U_buffer->values + i) = *(U->values + i);
    }

    matrix_free(L);
    matrix_free(U);
}

double calculate_determinant_lu(struct DoubleMatrix* matrix_ptr) {
    // Check square
    assert(matrix_ptr->numRows == matrix_ptr->numCols);

    int n = matrix_ptr->numRows;

    struct DoubleMatrix* L_matrix = matrix_factory(n, n, NULL);
    struct DoubleMatrix* U_matrix = matrix_factory(n, n, NULL);

    calculate_lu_matrices(matrix_ptr, L_matrix, U_matrix);
    
    double result = 1;
    for (int i = 0; i < n; i++) {
        result *= *matrix_get_element(L_matrix, i, i);
        result *= *matrix_get_element(U_matrix, i, i);
    }

    matrix_free(L_matrix);
    matrix_free(U_matrix);

    return result;
}