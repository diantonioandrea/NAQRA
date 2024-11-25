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


// Vector structure.

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


// Matrix structure.

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
    #ifndef NDEBUG // Integrity.
    assert(rows > 0);
    assert(columns > 0);
    #endif

    const register Natural N0 = rows;
    const register Natural N1 = columns;

    Matrix* M0 = (Matrix*) malloc(sizeof(Matrix));

    M0->Cm0 = (Complex*) calloc(N0 * N1, sizeof(Complex));
    M0->N0 = N0;
    M0->N1 = N1;

    return M0;
}


// Eigenvalues.

[[nodiscard]] Vector* Eigenvalues(const Matrix*);


// Output.

inline void PrintRowVector(const Vector* vector) { Pn_CrvN_0(vector->Cv0, vector->N0); }
inline void PrintColumnVector(const Vector* vector) { Pn_CcvN_0(vector->Cv0, vector->N0); }
inline void PrintMatrix(const Matrix* matrix) { Pn_CmNN_0(matrix->Cm0, matrix->N0, matrix->N1); }

#endif