/* 
 * PHYM0004 Assignment 1 / Jay Malhotra / 18/10/2021
 * 
 * determinant.h is a header that defines the two methods to calculate the
 * determinant of a DoubleMatrix type. One is a brute-force method using the
 * full cofactor formula, and one is a PLU decomposition method.
 */

/* Given a matrix m, calculate its determinant using the cofactor sum formula
 * a.k.a. 'brute force' or the 'traditional method' */
double calculate_determinant(struct DoubleMatrix* matrix_ptr);

/* Calculate the cofactor of a matrix element at row n and column m of a matrix.
 * This is in determinant.h rather than matrix_core.h as this operation requires
 * taking the determinant of the minor matrix. */
double calculate_cofactor(struct DoubleMatrix* matrix_ptr, int n, int m);

/* Calculate the LU decomposition for a given matrix. Since it is a bit
 * complicated to return two valeus (i.e. two matrices) this function is called
 * by first creating two empty matrices whose values should be filled with the
 * contents of the L and U matrix, respectively. */
void calculate_lu_matrices(struct DoubleMatrix* matrix_ptr, 
                           struct DoubleMatrix* L_buffer, struct DoubleMatrix* U_buffer);

/* Given a matrix m, calculate its determinant using the LU decomposition
 * method. */
double calculate_determinant_lu(struct DoubleMatrix* matrix_ptr);