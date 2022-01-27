/* 
 * PHYM0004 Assignment 1 / Jay Malhotra / 07/10/2021
 * 
 * matrix_core.h sets out the key matrix functions that define a matrix datatype
 * and enable basic operations with it, such as getting/setting values, output
 * to stdout, and addition/multiplication/etc.
 * 
 * Note that not all of the methods are used by the program's intended purpose
 * of calculating determinants -- in particular, some arithmetic operations such
 * as multiplication and addition are never used outside of unit testing; I
 * programmed the matrix implementation before really considering what the
 * methods would actually require, and those felt like methods that 'should'
 * be there.
 */

// Include guard
#ifndef MATRIX_CORE_H
#define MATRIX_CORE_H

#include <stdlib.h>
#include <stdbool.h>

/* 
 * =============================================================================
 * CORE METHODS
 * =============================================================================
 * Define the DoubleMatrix object and the function to create and free instances
 * of it. The struct should never be directly created and the factory should
 * always be used, so that memory can be allocated and numVals is correct.
 * 
 * It is always necessary to free a matrix once it is no longer needed, because
 * there is allocated memory that will stick around if you don't, and may cause
 * a memory leak. This will not be done automatically and so it is always the
 * caller's responsibility to do so.
 */

/* DoubleMatrix is the main type definition used to work with matrices. The
 * first two values are the dimension numbers. The fourth is numRows * numCols,
 * a.k.a. how many values are in the matrix. THe fourth is a pointer to the
 * start of an array of the matrix's values, which are 'unwrapped' in row order,
 * e.g. [a_11, a_12, a_21, a_22] for a 2x2 matrix. */
struct DoubleMatrix {
    int numRows;
    int numCols;
    int numVals;
    double* values;
};

/* matrix_factory serves as the way to create a DoubleMatrix struct -- given
 * dimensional information and optionally a pointer to an array of data to fill
 * the matrix with, it will allocate the memory and handle creation of the value
 * array pointer. 
 *
 * It will return a pointer to the created DoubleMatrix. If the data pointer is
 * null, the default allocation will be used and all the values of the matrix
 * will be zero. 
 *
 * Note that if the data pointed to by data_ptr is not the same size as the
 * matrix itself, undefined behaviour will occur. In particular, if data_ptr
 * contains less elements than the matrix, memory beyond the region of data_ptr
 * will be read and inserted into the matrix. */
struct DoubleMatrix* matrix_factory(int numRows, int numCols, double* data);

/* Given a pointer to a matrix, call free() appropriately to dispose of it. */
void matrix_free(struct DoubleMatrix* matrix_ptr);

/* 
 * =============================================================================
 * MISCELLANEOUS METHODS
 * =============================================================================
 */

/* Get the pointer to the value stored in the matrix at the intersection of row 
 * n and column m. Note that matrices are zero-indexed because that makes
 * the most sense when working with pointers and loops, in my opinion. */
double* matrix_get_element(struct DoubleMatrix* matrix_ptr, int n, int m);

/* Write to the specific point in a matrix at the intersection of row n and
 * column m. Note that matrices are zero-indexed because that makes
 * the most sense when working with pointers. */
void matrix_set_element(struct DoubleMatrix* matrix_ptr, int n, int m, double value);

/* Check if two matrices are equal, i.e. have the same dimensions and identical
 * elements to within a given relative tolerance. The tolerance is defined due
 * to the fact that floating point numbers are rarely ever exactly equal. */
bool matrix_check_equal(struct DoubleMatrix* a, struct DoubleMatrix* b, double reltol);

/* Output a matrix to stdout, showing each value to 3 decimal places. */
void matrix_print(struct DoubleMatrix* matrix_ptr);

/* 
 * =============================================================================
 * ARITHMETIC METHODS
 * =============================================================================
 * Each of these methods will create a new matrix rather than modifying the
 * existing one(s). Thus it will be necessary to free all matrices, i.e. if you
 * perform a calculation A * B = C, once you are finished you must still free
 * A, B, and C -- not just C.
 */

/* Multiply two matrices in the order a*b, and get a pointer to a new matrix
 * representing the result. Will fail an assert if the matrices cannot be
 * multiplied (i.e. A is (P by Q) and B is (R by S) and Q != R) */
struct DoubleMatrix* matrix_multiply(struct DoubleMatrix* a, struct DoubleMatrix* b);

/* Multiply all elements of the matrix by the given scalar, and get a pointer
 * to a new matrix representing the result. */
struct DoubleMatrix* matrix_multiply_scalar(double s, struct DoubleMatrix* a);

/* Add two matrices in the order a+b, and get a pointer to a new matrix
 * representing the result. Will fail an assert if the matrices are of different
 * sizes. */
struct DoubleMatrix* matrix_add(struct DoubleMatrix* a, struct DoubleMatrix* b);

/* Subtract two matrices in the order a-b, and get a pointer to a new matrix
 * representing the result. Will fail an assert if the matrices are of different
 * sizes. */
struct DoubleMatrix* matrix_subtract(struct DoubleMatrix* a, struct DoubleMatrix* b);

/* 
 * =============================================================================
 * DETERMINANT METHODS
 * =============================================================================
 * These are not traditional arithmetic operations, but rather operations
 * that are useful when calculating the determinant, either via brute-force or
 * LU decomposition. Similarly to the above methods, they produce a new matrix
 * rather than writing to the existing one.
 */

/* Get the minor matrix of the element at row n and column m. Used when
 * calculating the determinant 'traditionally'. */
struct DoubleMatrix* matrix_minor(struct DoubleMatrix* matrix_ptr, int n, int m);

/* Get the nth leading matrix of a square matrix given by matrix_ptr.
 * Constraints in n: 2 <= n < matrix size. */
struct DoubleMatrix* matrix_leading(struct DoubleMatrix* matrix_ptr, int n);

#endif // include guard end