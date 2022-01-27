/* 
 * PHYM0004 Assignment 1 / Jay Malhotra / 10/10/2021
 * 
 * test.c defines the unit tests for the program. It is executed by using the
 * 'test' make target in the makefile, which sidesteps main.c and uses this
 * alternative main() definition. It leverages cmocka as a unit testing
 * framework, which is used under the Apache License 2.0.
 * 
 * Please note that cmocka is not installed by default on most systems (it was
 * not on mine, but was in the repositories). Their website is here:
 * https://cmocka.org/
 */

#include <stdio.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>
#include <stdbool.h>
#include "matrix_core.h"
#include "determinant.h"

#define RELTOL 0.01

/* A test case that creates a matrix and checks its values are correctly
 * populated. */
static void test_create(void **state) {
    double matrix_data[4] = {3, 2, 1, 0};
    struct DoubleMatrix* matrix = matrix_factory(2, 2, matrix_data);

    assert_true(*matrix_get_element(matrix, 0, 0) == 3);
    assert_true(*matrix_get_element(matrix, 0, 1) == 2);
    assert_true(*matrix_get_element(matrix, 1, 0) == 1);
    assert_true(*matrix_get_element(matrix, 1, 1) == 0);

    matrix_free(matrix);
}

/* A test case that creates two identical matrices to test the check_equal
 * function. */
static void test_check_equal(void **state) {
    double matrix_data_1[4] = {1, 2, 3, 4};
    struct DoubleMatrix* matrix_1 = matrix_factory(2, 2, matrix_data_1);

    double matrix_data_2[4] = {1, 2, 3, 4};
    struct DoubleMatrix* matrix_2 = matrix_factory(2, 2, matrix_data_2);

    assert_true(matrix_check_equal(matrix_1, matrix_2, RELTOL));

    matrix_free(matrix_1);
    matrix_free(matrix_2);
}

/* A test case that multiplies two 2x2 matrices to check the multiply
 * function. */
static void test_multiply_1(void **state) {
    double matrix_data_1[4] = {1, 2, 3, 4};
    struct DoubleMatrix* matrix_1 = matrix_factory(2, 2, matrix_data_1);

    double matrix_data_2[4] = {4, 3, 2, 1};
    struct DoubleMatrix* matrix_2 = matrix_factory(2, 2, matrix_data_2);

    struct DoubleMatrix* matrix_result = matrix_multiply(matrix_1, matrix_2);
    
    double matrix_data_exp[4] = {8, 5, 20, 13};
    struct DoubleMatrix* matrix_expected = matrix_factory(2, 2, matrix_data_exp);

    assert_true(matrix_check_equal(matrix_expected, matrix_result, RELTOL));

    matrix_free(matrix_1);
    matrix_free(matrix_2);
    matrix_free(matrix_result);
    matrix_free(matrix_expected);
}

static void test_multiply_2(void **state) {
    double matrix_data_1[3] = {7, 6, 2};
    struct DoubleMatrix* matrix_1 = matrix_factory(3, 1, matrix_data_1);

    double matrix_data_2[3] = {12, 40, 8};
    struct DoubleMatrix* matrix_2 = matrix_factory(1, 3, matrix_data_2);

    struct DoubleMatrix* matrix_result = matrix_multiply(matrix_1, matrix_2);
    
    double matrix_data_exp[9] = {84, 280, 56, 72, 240, 48, 24, 80, 16};
    struct DoubleMatrix* matrix_expected = matrix_factory(3, 3, matrix_data_exp);

    /* matrix_print(matrix_result);
    matrix_print(matrix_expected); */

    assert_true(matrix_check_equal(matrix_expected, matrix_result, RELTOL));

    matrix_free(matrix_1);
    matrix_free(matrix_2);
    matrix_free(matrix_result);
    matrix_free(matrix_expected);
}

static void test_multiply_scalar(void **state) {
    double matrix_data[6] = {3, 1, 4, 1, 5, 9};
    struct DoubleMatrix* matrix_1 = matrix_factory(3, 2, matrix_data);
    
    struct DoubleMatrix* matrix_result = matrix_multiply_scalar(2.0, matrix_1);

    double matrix_data_exp[6] = {6, 2, 8, 2, 10, 18};
    struct DoubleMatrix* matrix_expected = matrix_factory(3, 2, matrix_data_exp);


    assert_true(matrix_check_equal(matrix_expected, matrix_result, RELTOL));

    matrix_free(matrix_1);
    matrix_free(matrix_result);
    matrix_free(matrix_expected);
}

static void test_add(void **state) {
    double matrix_data_1[9] = {33, 25, 69, 7, 59, 80, 37, 36, 59};
    struct DoubleMatrix* matrix_1 = matrix_factory(3, 3, matrix_data_1);

    double matrix_data_2[9] = {5.4, 5.2, 3.5, 0.4, 9.6, 2.6, 4.5, 4.6, 8.7};
    struct DoubleMatrix* matrix_2 = matrix_factory(3, 3, matrix_data_2);

    struct DoubleMatrix* matrix_result = matrix_add(matrix_1, matrix_2);

    double matrix_data_exp[9] = {38.4, 30.2, 72.5, 7.4, 68.6, 82.6, 41.5, 40.6, 67.7};
    struct DoubleMatrix* matrix_expected = matrix_factory(3, 3, matrix_data_exp);

    assert_true(matrix_check_equal(matrix_expected, matrix_result, RELTOL));

    matrix_free(matrix_1);
    matrix_free(matrix_2);
    matrix_free(matrix_result);
    matrix_free(matrix_expected);
}

static void test_subtract(void **state) {
    double matrix_data_1[9] = {3, 72, 3, 12, 83, 3, 17, 43, 14};
    struct DoubleMatrix* matrix_1 = matrix_factory(3, 3, matrix_data_1);

    double matrix_data_2[9] = {0.26, 0.70, 0.94, 0.91, 0.46, 0.97, 0.15, 0.60, .91};
    struct DoubleMatrix* matrix_2 = matrix_factory(3, 3, matrix_data_2);

    struct DoubleMatrix* matrix_result = matrix_subtract(matrix_1, matrix_2);

    double matrix_data_exp[9] = {2.74, 71.3, 2.06, 11.09, 82.54, 2.03, 16.85, 42.4, 13.09};
    struct DoubleMatrix* matrix_expected = matrix_factory(3, 3, matrix_data_exp);

    assert_true(matrix_check_equal(matrix_expected, matrix_result, RELTOL));

    matrix_free(matrix_1);
    matrix_free(matrix_2);
    matrix_free(matrix_result);
    matrix_free(matrix_expected);
}

static void test_minor(void **state) {
    double matrix_data[9] = {6, 4, 3, 9, 2, 5, 1, 7, 8};
    struct DoubleMatrix* matrix = matrix_factory(3, 3, matrix_data);

    // Minor of central 2
    struct DoubleMatrix* matrix_result = matrix_minor(matrix, 1, 1);

    double matrix_data_exp[4] = {6, 3, 1, 8};
    struct DoubleMatrix* matrix_expected = matrix_factory(2, 2, matrix_data_exp);

    assert_true(matrix_check_equal(matrix_expected, matrix_result, RELTOL));

    /* matrix_print(matrix);
     * matrix_print(matrix_result); */

    matrix_free(matrix);
    matrix_free(matrix_result);
    matrix_free(matrix_expected);
}

// Test 'regular' determinant
static void test_determinant_1(void **state) {
    double matrix_data[4] = {1, 2, 3, 4};
    struct DoubleMatrix* matrix = matrix_factory(2, 2, matrix_data);

    double result = calculate_determinant(matrix);
    double expected_result = -2;

    assert_true(abs(1 - result / expected_result) < RELTOL);
    matrix_free(matrix);
}


static void test_determinant_2(void **state) {
    double matrix_data[9] = {2, 3, 5, 20, 9, 5, 3, 9, 15};
    struct DoubleMatrix* matrix = matrix_factory(3, 3, matrix_data);

    double result = calculate_determinant(matrix);
    double expected_result = 90;

    assert_true(abs(1 - result / expected_result) < RELTOL);
    matrix_free(matrix);
}

static void test_lu_decomp_1(void **state) {
    double matrix_data[9] = {1, 4, 3, 2, 4, 5, 1, 3, 4};
    struct DoubleMatrix* matrix = matrix_factory(3, 3, matrix_data);
    
    struct DoubleMatrix* L_matrix = matrix_factory(3, 3, NULL);
    struct DoubleMatrix* U_matrix = matrix_factory(3, 3, NULL);

    calculate_lu_matrices(matrix, L_matrix, U_matrix);

    // Easiest way to check it's correct is to multiply them, which should give
    // back the original matrix
    struct DoubleMatrix* result = matrix_multiply(L_matrix, U_matrix);

    assert_true(matrix_check_equal(matrix, result, RELTOL));

    matrix_free(matrix);
    matrix_free(L_matrix);
    matrix_free(U_matrix);
    matrix_free(result);
}

static void test_lu_decomp_2(void **state) {
    double matrix_data[16] = {10, 5, 4, 2, 8, 5, 9, 2, 14, 7, 6, 13, 4, 7, 18, 9};
    struct DoubleMatrix* matrix = matrix_factory(4, 4, matrix_data);
    
    struct DoubleMatrix* L_matrix = matrix_factory(4, 4, NULL);
    struct DoubleMatrix* U_matrix = matrix_factory(4, 4, NULL);

    calculate_lu_matrices(matrix, L_matrix, U_matrix);

    struct DoubleMatrix* result = matrix_multiply(L_matrix, U_matrix);

    assert_true(matrix_check_equal(matrix, result, RELTOL));

    matrix_free(matrix);
    matrix_free(L_matrix);
    matrix_free(U_matrix);
    matrix_free(result);
}

static void test_determinant_lu_1(void **state) {
    double matrix_data[4] = {1, 2, 3, 4};
    struct DoubleMatrix* matrix = matrix_factory(2, 2, matrix_data);

    double result = calculate_determinant_lu(matrix);
    double expected_result = -2;

    assert_true(abs(1 - result / expected_result) < RELTOL);
    matrix_free(matrix);
}

static void test_determinant_lu_2(void **state) {
    double matrix_data[9] = {2, 3, 5, 20, 9, 5, 3, 9, 15};
    struct DoubleMatrix* matrix = matrix_factory(3, 3, matrix_data);

    double result = calculate_determinant_lu(matrix);
    double expected_result = 90;

    assert_true(abs(1 - result / expected_result) < RELTOL);
    matrix_free(matrix);
}

static void test_leading_matrix(void **state) {
    double matrix_data[9] = {13, 4, 5, 10, 48, 2, 9, 5, 17};
    struct DoubleMatrix* matrix = matrix_factory(3, 3, matrix_data);

    struct DoubleMatrix* leading_matrix = matrix_leading(matrix, 2);

    double matrix_data_exp[4] = {13, 4, 10, 48};
    struct DoubleMatrix* matrix_expected = matrix_factory(2, 2, matrix_data_exp);

    assert_true(matrix_check_equal(leading_matrix, matrix_expected, RELTOL));

    matrix_free(matrix);
    matrix_free(leading_matrix);
    matrix_free(matrix_expected);
}

int main(void) {
    const struct CMUnitTest tests[] = {
        cmocka_unit_test(test_create),
        cmocka_unit_test(test_check_equal),
        cmocka_unit_test(test_multiply_1),
        cmocka_unit_test(test_multiply_2),
        cmocka_unit_test(test_multiply_scalar),
        cmocka_unit_test(test_add),
        cmocka_unit_test(test_subtract),
        cmocka_unit_test(test_minor),
        cmocka_unit_test(test_determinant_1),
        cmocka_unit_test(test_determinant_2),
        cmocka_unit_test(test_lu_decomp_1),
        cmocka_unit_test(test_lu_decomp_2),
        cmocka_unit_test(test_determinant_lu_1),
        cmocka_unit_test(test_determinant_lu_2),
        cmocka_unit_test(test_leading_matrix)
    };
    return cmocka_run_group_tests(tests, NULL, NULL);
}