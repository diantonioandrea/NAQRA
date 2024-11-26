/**
 * @file Interface.h
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief NAQRA interface.
 * @date 2024-11-25
 * 
 * @copyright Copyright (c) 2024
 * 
 * Different naming convention than standard NAQRA methods and variables.
 * No differences in implementation.
 */

#ifndef NAQRA_INTERFACE_H
#define NAQRA_INTERFACE_H

// Assertions.
#include <assert.h>

// Matrices.
#include "./Matrix.h"


// Complex.

/**
 * @brief Complex constructor.
 * 
 * @param real Real part.
 * @param imaginary Imaginary part.
 * @return Complex Complex number.
 */
static inline Complex NewComplex(const Real real, const Real imaginary) {
    const register Real R0 = real;
    const register Real R1 = imaginary;

    return C_RR_C(R0, R1);
}


// Vectors.

// Complex vector.
typedef struct {

    // Entries.
    Complex* Cv0;

    // Size.
    Natural N0;

} Vector;

/**
 * @brief Vector constructor.
 * 
 * @param size Vector's size. 
 * @return Vector* Vector.
 */
[[nodiscard]] static inline Vector* NewVector(const Natural size) {
    #ifndef NDEBUG // Integrity.
    assert(size > 0);
    #endif

    const register Natural N0 = size;

    Vector* V0 = (Vector*) malloc(sizeof(Vector));

    V0->Cv0 = (Complex*) calloc(N0, sizeof(Complex));
    V0->N0 = N0;

    return V0;
}

/**
 * @brief Vector getter.
 * 
 * @param vector Vector.
 * @param index Index.
 * @return Complex Complex number.
 */
static inline Complex GetVectorAt(const Vector* vector, const Natural index) {
    const register Natural N0 = index;

    #ifndef NDEBUG // Integrity.
    assert(N0 < vector->N0);
    #endif

    return vector->Cv0[N0];
}

/**
 * @brief Vector setter.
 * 
 * @param vector Vector.
 * @param index Index.
 * @param complex Complex number.
 */
static inline void SetVectorAt(Vector* vector, const Natural index, const Complex complex) {
    const register Natural N0 = index;
    const register Complex C0 = complex;

    #ifndef NDEBUG // Integrity.
    assert(N0 < vector->N0);
    #endif

    vector->Cv0[N0] = C0;
}


// Matrices.

// Complex matrix.
typedef struct {

    // Entries.
    Complex* Cm0;

    // Rows.
    Natural N0;

    // Columns.
    Natural N1;

} Matrix;

/**
 * @brief Matrix constructor.
 * 
 * @param rows Matrix's rows.
 * @param columns Matrix's columns.
 * @return Matrix* Matrix.
 */
[[nodiscard]] static inline Matrix* NewMatrix(const Natural rows, const Natural columns) {
    const register Natural N0 = rows;
    const register Natural N1 = columns;

    #ifndef NDEBUG // Integrity.
    assert(N0 > 0);
    assert(N1 > 0);
    #endif

    Matrix* M0 = (Matrix*) malloc(sizeof(Matrix));

    M0->Cm0 = (Complex*) calloc(N0 * N1, sizeof(Complex));
    M0->N0 = N0;
    M0->N1 = N1;

    return M0;
}

/**
 * @brief Matrix getter.
 * 
 * @param matrix Matrix.
 * @param row Row index.
 * @param column Column index.
 * @return Complex Complex number.
 */
static inline Complex GetMatrixAt(const Matrix* matrix, const Natural row, const Natural column) {
    const register Natural N0 = row;
    const register Natural N1 = column;

    #ifndef NDEBUG // Integrity.
    assert(N0 < matrix->N0);
    assert(N1 < matrix->N1);
    #endif

    return matrix->Cm0[N1 * matrix->N0 + N0];
}

/**
 * @brief Matrix setter.
 * 
 * @param matrix Matrix.
 * @param row Row index.
 * @param column Column index.
 * @param complex Complex number.
 */
static inline void SetMatrixAt(Matrix* matrix, const Natural row, const Natural column, const Complex complex) {
    const register Natural N0 = row;
    const register Natural N1 = column;
    const register Complex C0 = complex;

    #ifndef NDEBUG // Integrity.
    assert(N0 < matrix->N0);
    assert(N1 < matrix->N1);
    #endif

    matrix->Cm0[N1 * matrix->N0 + N0] = C0;
}


// Eigenvalues.

[[nodiscard]] Vector* Eigenvalues(const Matrix*);


// Output.

inline void PrintRowVector(const Vector* vector) { Pn_CrvN_0(vector->Cv0, vector->N0); }
inline void PrintColumnVector(const Vector* vector) { Pn_CcvN_0(vector->Cv0, vector->N0); }
inline void PrintMatrix(const Matrix* matrix) { Pn_CmNN_0(matrix->Cm0, matrix->N0, matrix->N1); }

#endif