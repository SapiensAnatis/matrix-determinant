/* 
 * PHYM0004 Assignment 1 / Jay Malhotra / 07/10/2021
 *
 * matrix_core.c contains the implementations of the methods defined in the
 * header file of the same name (matrix_core.h). The header file contains more
 * detail about these definitions with comments explaining the purpose of each
 * function and definition.
 */

#include "matrix_core.h"
#include <stdio.h>
#include <assert.h>

struct DoubleMatrix* matrix_factory(int numRows, int numCols, double* data) {
    assert(numRows > 0);
    assert(numCols > 0);

    struct DoubleMatrix* result = (struct DoubleMatrix*)malloc(sizeof(struct DoubleMatrix));
    int numVals = numRows * numCols;

    result->numRows = numRows;
    result->numCols = numCols;
    result->numVals = numVals;
    
    // Want to initialize all values to zero rather than having it filled with
    // garbage values by default
    result->values = calloc(numVals, sizeof(double));

    if (data) {
        for (int i = 0; i < numVals; i++) {
            *(result->values + i) = *(data + i);
        }
    }

    return result;
}

void matrix_free(struct DoubleMatrix* matrix_ptr) {
    free(matrix_ptr->values);
    free(matrix_ptr);
}

double* matrix_get_element(struct DoubleMatrix* matrix_ptr, int n, int m) {
    // Assert that these are valid indices
    assert(n >= 0 && n <= matrix_ptr->numRows);
    assert(m >= 0 && m <= matrix_ptr->numCols);

    double* value_ptr = matrix_ptr->values;
    // Go down to the correct row first; row width = numCols, skip forward
    value_ptr += (matrix_ptr->numCols) * n;
    // Then go along that row to the right column
    value_ptr += m;

    return value_ptr;
}

void matrix_set_element(struct DoubleMatrix* matrix_ptr, int n, int m, double value) {
    double* value_ptr = matrix_get_element(matrix_ptr, n, m);
    *value_ptr = value;
}

bool matrix_check_equal(struct DoubleMatrix* a, struct DoubleMatrix* b, double reltol) {
    if (a->numRows != b->numRows || a->numCols != b->numCols) { return false; }
    
    double reldiff;
    for (int i = 0; i < a->numVals; i++) {
        reldiff = 1 - (*(a->values + i) / *(b->values + i));    
        if (reldiff > reltol) {
            return false;
        }
    }

    return true;
}

void matrix_print(struct DoubleMatrix* matrix_ptr) {
    // TODO maybe if I care enough:
    // Add padding intelligently so the columns aren't misaligned when numbers
    // are of differing magnitudes
    
    // e.g. stop it from looking like this
    // [ 8.000 280.000 561.000 ]
    // [ 72.000 2.000 48.000 ]
    // [ 24.000 80.000 16.000 ]

    // Add some flavour text, otherwise printing multiple matrix in a row can
    // make them look like they're all one matrix
    printf("%dx%d matrix:\n", matrix_ptr->numRows, matrix_ptr->numCols);

    for (int i = 0; i < matrix_ptr->numRows; i++) {
        printf("\t[ ");
        for (int j = 0; j < matrix_ptr->numCols; j++) {
            printf("%.3f ", *matrix_get_element(matrix_ptr, i, j));
        }
        printf("]\n");
    }
}

struct DoubleMatrix* matrix_multiply(struct DoubleMatrix* a, struct DoubleMatrix* b) {
    // Check multiplication is possible
    assert(a->numCols == b->numRows);

    double new_data[a->numRows * b->numCols];
    int buffer_offset = 0;
    struct DoubleMatrix* matrix_result;

    double tmp_result;
    double operand_a;
    double operand_b;

    // First populate the double array with the results. Sorry if this is a
    // little opaque, but this is just how I personally do it in my head...pick
    // a row in A, and go sideways along it and down a column in B. Repeat for
    // the next column of B, then go down a row in A and do the same.
    for (int row_a = 0; row_a < a->numRows; row_a++) {
        for (int col_b = 0; col_b < b->numCols; col_b++) {
            // This loop represents the calculation of one element in the result
            tmp_result = 0;
            for (int i = 0; i < a->numCols; i++) {
                // i is the column of a and the row of b
                operand_a = *matrix_get_element(a, row_a, i);
                operand_b = *matrix_get_element(b, i, col_b);

                tmp_result += operand_a * operand_b;
            }

            new_data[buffer_offset] = tmp_result;
            buffer_offset++;
        }
    }

    // Then create the new matrix
    matrix_result = matrix_factory(a->numRows, b->numCols, new_data);

    return matrix_result;
}

struct DoubleMatrix* matrix_multiply_scalar(double s, struct DoubleMatrix* a) {
    // It is probably preferable not to directly access the data pointer of the
    // matrix, however for simple use-cases involving iteration over all values
    // without any regard to their position, it's often far easier and more
    // readable to do it this way rather than get_element in a nested for loop.
    // The same logic applies for the addition and subtraction methods below.
    
    double new_data[a->numVals];

    for(int i = 0; i < a->numVals; i++) {
        *(new_data + i) = s * *(a->values + i);
    }

    return matrix_factory(a->numRows, a->numCols, new_data);
}

struct DoubleMatrix* matrix_add(struct DoubleMatrix* a, struct DoubleMatrix* b) {
    // Cannot add matrices if they have different dimensions
    assert(a->numRows == b->numRows && a->numCols == b->numCols);

    double new_data[a->numVals];

    for (int i = 0; i < a->numVals; i++) {
        *(new_data + i) = *(a->values + i) + *(b->values + i);
    }

    return matrix_factory(a->numRows, a->numCols, new_data);
}

struct DoubleMatrix* matrix_subtract(struct DoubleMatrix* a, struct DoubleMatrix* b) {
    // Very much the same as above...no use in being clever and doing a scalar
    // multiplication by -1 and adding as that'll probably cost more CPU cycles
    // so I'll just copy and paste the code

    assert(a->numRows == b->numRows && a->numCols == b->numCols);

    double new_data[a->numVals];

    for (int i = 0; i < a->numVals; i++) {
        // Witness the one changed line
        *(new_data + i) = *(a->values + i) - *(b->values + i);
    }

    return matrix_factory(a->numRows, a->numCols, new_data);
}

struct DoubleMatrix* matrix_minor(struct DoubleMatrix* matrix_ptr, int n, int m) {
    // 'Black out' the row and column
    // Assert that these are valid indices
    assert(n >= 0 && n <= matrix_ptr->numRows);
    assert(m >= 0 && m <= matrix_ptr->numCols);

    int new_rows = matrix_ptr->numRows - 1;
    int new_cols = matrix_ptr->numCols - 1;
    int buffer_offset = 0;

    double new_data[new_rows * new_cols];

    for (int row = 0; row < matrix_ptr->numRows; row++) {
        // Skip over the row, or column, that matches the index of the element
        if (row == n) { continue; }
        for (int col = 0; col < matrix_ptr->numCols; col++) {
            if (col == m) { continue; }
            new_data[buffer_offset] = *matrix_get_element(matrix_ptr, row, col);
            buffer_offset++;
        }
    }

    return matrix_factory(new_rows, new_cols, new_data);
}

struct DoubleMatrix* matrix_leading(struct DoubleMatrix* matrix_ptr, int n) {
    // Assert square
    assert(matrix_ptr->numRows == matrix_ptr->numCols);
    // Assert that n is valid
    assert(n >= 2 && n <= matrix_ptr->numRows);

    double new_data[n*n];
    int buffer_offset = 0;

    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            new_data[buffer_offset] = *matrix_get_element(matrix_ptr, row, col);
            buffer_offset++;
        }
    }

    return matrix_factory(n, n, new_data);
}