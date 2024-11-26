/**
 * @file NAQRA_Interface.c
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief include/Interface.h implementation.
 * @date 2024-11-25
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "../include/Interface.h"

// Eigenvalues.

/**
 * @brief Matrix eigenvalues.
 * 
 * @param matrix Square matrix.
 * @return Vector* Vector.
 */
[[nodiscard]] Vector* Eigenvalues(const Matrix* matrix) {
    #ifndef NDEBUG // Integrity check.
    assert(matrix->N0 == matrix->N1);
    #endif

    const register Natural N0 = matrix->N0;
    register Natural N1 = 0;

    Complex* Cm0 = (Complex*) malloc(N0 * N0 * sizeof(Complex)); // Matrix copy.
    Cp_CvtCvN_0(Cm0, matrix->Cm0, N0 * N0);

    Hsn_CqtN_0(Cm0, N0); // Hessenberg.
    Eig_ChsnqtN_0(Cm0, N0); // Eigenvalues.

    Vector* V0 = NewVector(N0);

    for(; N1 < N0; ++N1) // Eigenvalues copy.
        V0->Cv0[N1] = Cm0[N1 * (N0 + 1)];

    free(Cm0);

    return V0;
}