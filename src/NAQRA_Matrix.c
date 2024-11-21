/**
 * @file NAQRA_Matrix.c
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief include/Matrix.h implementation.
 * @date 2024-11-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <stdlib.h>
#include "../include/Matrix.h"

// Householder products.

/**
 * @brief Householder Left [Hsl] 
 * 
 * @param Cqt0 Complex Sqaure Matrix [Cq], Target [t].
 * @param Cv0 Complex Vector [Cv].
 * @param N0 Rows and Columns [N].
 * @param N1 Entries [N].
 */
void Hsl_CqtCvNN_0(Complex* Cqt0, const Complex* Cv0, const Natural N0, const Natural N1) {
    register Natural N2 = 0, N3, N4;
    const register Natural N5 = N0 - N1;
    register Complex C0 = {0.0, 0.0};

    if(N0 < 5)
        for(; N2 < N0; ++N2) {
            N3 = N5;
            N4 = N2 * N0 + N5;

            for(C0 = C_R_C(0.0); N3 < N0; ++N3, ++N4)
                C0 = A_CC_C(C0, M_CcjC_C(Cv0[N3], Cqt0[N4]));

            C0 = M_CR_C(C0, 2.0);

            for(N3 = N5, N4 = N2 * N0 + N5; N3 < N0; ++N3, ++N4)
                Cqt0[N4] = S_CC_C(Cqt0[N4], M_CC_C(Cv0[N3], C0));
        }
    else {
        register Complex C1 = {0.0, 0.0};
        register Complex C2 = {0.0, 0.0};
        register Complex C3 = {0.0, 0.0};
        register Complex C4 = {0.0, 0.0};

        for(; N2 < N0; ++N2) {
            N3 = N5;
            N4 = N2 * N0 + N5;

            C1 = M_CcjC_C(Cv0[N3], Cqt0[N4]);
            C2 = M_CcjC_C(Cv0[N3 + 1], Cqt0[N4 + 1]);
            C3 = M_CcjC_C(Cv0[N3 + 2], Cqt0[N4 + 2]);
            C4 = M_CcjC_C(Cv0[N3 + 3], Cqt0[N4 + 3]);

            N3 += 4;
            N4 += 4;

            for(; N3 + 3 < N0; N3 += 4, N4 += 4) {
                C1 = A_CC_C(C1, M_CcjC_C(Cv0[N3], Cqt0[N4]));
                C2 = A_CC_C(C2, M_CcjC_C(Cv0[N3 + 1], Cqt0[N4 + 1]));
                C3 = A_CC_C(C3, M_CcjC_C(Cv0[N3 + 2], Cqt0[N4 + 2]));
                C4 = A_CC_C(C4, M_CcjC_C(Cv0[N3 + 3], Cqt0[N4 + 3]));
            }

            for(; N3 < N0; ++N3, ++N4)
                C1 = A_CC_C(C1, M_CcjC_C(Cv0[N3], Cqt0[N4]));

            C0 = M_CR_C(A_CC_C(A_CC_C(C1, C2), A_CC_C(C3, C4)), 2.0);

            for(N3 = N5, N4 = N2 * N0 + N5; N3 + 3 < N0; N3 += 4, N4 += 4) {
                Cqt0[N4] = S_CC_C(Cqt0[N4], M_CC_C(Cv0[N3], C0));
                Cqt0[N4 + 1] = S_CC_C(Cqt0[N4 + 1], M_CC_C(Cv0[N3 + 1], C0));
                Cqt0[N4 + 2] = S_CC_C(Cqt0[N4 + 2], M_CC_C(Cv0[N3 + 2], C0));
                Cqt0[N4 + 3] = S_CC_C(Cqt0[N4 + 3], M_CC_C(Cv0[N3 + 3], C0));
            }

            for(N3 = N5, N4 = N2 * N0 + N5; N3 < N0; ++N3, ++N4)
                Cqt0[N4] = S_CC_C(Cqt0[N4], M_CC_C(Cv0[N3], C0));
        }
    }
}

/**
 * @brief Householder Right [Hsr] 
 * 
 * @param Cqt0 Complex Sqaure Matrix [Cq], Target [t].
 * @param Cv0 Complex Vector [Cv].
 * @param N0 Rows and Columns [N].
 * @param N1 Entries [N].
 */
void Hsr_CqtCvNN_0(Complex* Cqt0, const Complex* Cv0, const Natural N0, const Natural N1) {
    register Natural N2 = 0, N3;
    const register Natural N4 = N0 - N1;

    for(; N2 < N0; ++N2) {
        register Complex C0 = {0.0, 0.0};

        for(N3 = N4; N3 < N0; ++N3)
            C0 = A_CC_C(C0, M_CC_C(Cqt0[N3 * N0 + N2], Cv0[N3]));

        C0 = M_CR_C(C0, 2.0);

        for(N3 = N4; N3 < N0; ++N3)
            Cqt0[N3 * N0 + N2] = S_CC_C(Cqt0[N3 * N0 + N2], M_CCcj_C(C0, Cv0[N3]));
    }
}

// Givens products.

/**
 * @brief Givens Left [Gvl].
 * 
 * @param Chsnqt0 Chsnqt0 Complex Hessenberg Square Matrix [Chsnq], Target [t].
 * @param C0 Complex Number [C].
 * @param C1 Complex Number [C].
 * @param N0 Rows and Columns [N].
 * @param N1 Index [N].
 */
void Gvl_ChsnqtCCNN_0(Complex* Chsnqt0, const Complex C0, const Complex C1, const Natural N0, const Natural N1) {
    register Natural N3 = (N1 + 1) * (N0 + 1) - 1;
    register Natural N4 = N3 + 1;

    for(; N3 < N0 * N0; N3 += N0, N4 += N0) {
        const register Complex C2 = Chsnqt0[N3];
        const register Complex C3 = Chsnqt0[N4];

        Chsnqt0[N3] = A_CC_C(M_CcjC_C(C0, C2), M_CcjC_C(C1, C3));
        Chsnqt0[N4] = S_CC_C(M_CC_C(C0, C3), M_CC_C(C1, C2));
    }
}

/**
 * @brief Hermitian Givens Right [Gvrtr].
 * 
 * @param Chsnqt0 Chsnqt0 Complex Hessenberg Square Matrix [Chsnq], Target [t].
 * @param C0 Complex Number [C].
 * @param C1 Complex Number [C].
 * @param N0 Rows and Columns [N].
 * @param N1 Index [N].
 */
void Gvrhr_ChsnqtCCNN_0(Complex* Chsnqt0, const Complex C0, const Complex C1, const Natural N0, const Natural N1) {
    register Natural N3 = N1 * N0;
    register Natural N4 = (N1 + 1) * N0;

    for(; N3 < N1 * (N0 + 1) + 2; ++N3, ++N4) {
        const register Complex C2 = Chsnqt0[N3];
        const register Complex C3 = Chsnqt0[N4];

        Chsnqt0[N3] = A_CC_C(M_CC_C(C2, C0), M_CC_C(C3, C1));
        Chsnqt0[N4] = S_CC_C(M_CCcj_C(C3, C0), M_CCcj_C(C2, C1));
    }
}

// Hessenberg form.

/**
 * @brief Hessenberg form [Hsn].
 * 
 * @param Cqt0 Complex Sqaure Matrix [Cq], Target [t].
 * @param N0 Rows and Columns [N].
 */
void Hsn_CqtN_0(Complex* Cqt0, const Natural N0) {
    register Natural N1 = 0;
    register Complex* Cv0 = (Complex*) calloc(N0, sizeof(Complex));
    register Complex* Cv1 = (Complex*) calloc(N0 - 2, sizeof(Complex)); // Zeros.

    for(; N1 < N0 - 2; ++N1) {

        const register Natural N3 = N0 - N1 - 1; // Entries.
        
        // Householder vector.

        Cp_CvtCvN_0(Cv0 + N1 + 1, Cqt0 + N1 * (N0 + 1) + 1, N3); // Copy.
        Cv0[N1 + 1] = S_CC_C(Cv0[N1 + 1], M_CR_C(Nzd2_C_C(Cv0[N1 + 1]), N2_CvN_R(Cv0 + N1 + 1, N3))); // Direction.
        Nz2_CvN_0(Cv0 + N1 + 1, N3); // Normalization.

        // Householder products.

        Hsl_CqtCvNN_0(Cqt0, Cv0, N0, N3);
        Hsr_CqtCvNN_0(Cqt0, Cv0, N0, N3);

        // Zeroing.

        Cp_CvtCvN_0(Cqt0 + N1 * (N0 + 1) + 2, Cv1, N0 - N1 - 2);
    }

    free(Cv0);
    free(Cv1);
}

// QR algorithm.

/**
 * @brief Eigenvalue [Eig].
 * 
 * @param C0 Complex Number [C], Matrix element.
 * @param C1 Complex Number [C], Matrix element.
 * @param C2 Complex Number [C], Matrix element.
 * @param C3 Complex Number [C], Matrix element.
 * @return Complex* Complex Vector [Cv].
 */
[[nodiscard]] Complex* Eig_CCCC_C(const Complex C0, const Complex C1, const Complex C2, const Complex C3) {
    const register Complex C4 = A_CC_C(C0, C3);
    const register Complex C5 = S_CC_C(M_CC_C(C0, C3), M_CC_C(C1, C2));
    const register Complex C6 = S_CC_C(Sq_C_C(C4), M_CR_C(C5, 4.0));
    const register Complex C7 = D_CR_C(C4, 2.0);

    register Complex* Cv0 = Nrt_C_Cv(C6, 2);

    Cv0[0] = A_CC_C(C7, D_CR_C(Cv0[0], 2.0));
    Cv0[1] = A_CC_C(C7, D_CR_C(Cv0[1], 2.0));

    return Cv0;
}

/**
 * @brief Eigenvalues [Eig].
 * 
 * @param Chsnqt0 Complex Hessenberg Square Matrix [Chsnq], Target [t].
 * @param N0 Rows and Columns [N].
 */
void Eig_ChsnqtN_0(Complex *Chsnqt0, const Natural N0) {
    register Complex* Cv0 = (Complex*) calloc(2 * (N0 - 1), sizeof(Complex));
    register Complex* Cv1 = (Complex*) calloc(2, sizeof(Complex));

    register Natural N1 = 0, N2, N3 = N0 - 1, N4;
    register Real R0;

    #ifdef VERBOSE
    printf("--- QR Algorithm\n");
    #endif

    for(; N1 < ITM0; ++N1) {

        if(N3 == 0) break; // Stop.
        if(N2_C_R(Chsnqt0[(N3 - 1) * (N0 + 1) + 1]) <= TOL0) { --N3; continue; } // Deflation.

        *Cv1 = *Eig_CCCC_C(Chsnqt0[(N3 - 1) * (N0 + 1)], Chsnqt0[N3 * N0 + N3 - 1], Chsnqt0[(N3 - 1) * (N0 + 1) + 1], Chsnqt0[N3 * (N0 + 1)]); // Double Wilkinson's shift.

        for(N4 = 0; N4 < 2; ++N4) { // Double shift.

            for(N2 = 0; N2 < N3 + 1; ++N2) // Shift (-).
                Chsnqt0[N2 * (N0 + 1)] = S_CC_C(Chsnqt0[N2 * (N0 + 1)], Cv1[N4]);

            for(N2 = 0; N2 < N3; ++N2) { // QR.
                Cp_CvtCvN_0(Cv0 + 2 * N2, Chsnqt0 + N2 * (N0 + 1), 2); // Copy.

                R0 = N2_CvN_R(Cv0 + 2 * N2, 2); // Norm.
                D_CvR_0(Cv0 + 2 * N2, R0, 2); // Normalization.

                // First product.
                Chsnqt0[N2 * (N0 + 1)] = C_R_C(R0);
                Chsnqt0[N2 * (N0 + 1) + 1] = C_R_C(0.0);

                // Other products.
                Gvl_ChsnqtCCNN_0(Chsnqt0, Cv0[2 * N2], Cv0[2 * N2 + 1], N0, N2);
            }

            for(N2 = 0; N2 < N3; ++N2) // RQ.
                Gvrhr_ChsnqtCCNN_0(Chsnqt0, Cv0[2 * N2], Cv0[2 * N2 + 1], N0, N2);

            for(N2 = 0; N2 < N3 + 1; ++N2) // Shift (+).
                Chsnqt0[N2 * (N0 + 1)] = A_CC_C(Chsnqt0[N2 * (N0 + 1)], Cv1[N4]);
        }
    }

    #ifdef VERBOSE
    printf("Exited after %zu iterations.\n", N1 + 1);
    printf("---\n");
    #endif

    free(Cv0);
    free(Cv1);
}

// Output.

/**
 * @brief Print with new line [Pn].
 * 
 * @param Cm0 Complex Matrix [Cm].
 * @param N0 Rows [N].
 * @param N1 Columns [N].
 */
void Pn_CmNN_0(const Complex* Cm0, const Natural N0, const Natural N1) {
    register Natural N2 = 0, N3;

    printf("--- Matrix\n");

    for(; N2 < N0; ++N2) {
        for(N3 = 0; N3 < N1 - 1; ++N3)
            P_C_0(Cm0[N2 + N3 * N0]);

        Pn_C_0(Cm0[N2 + (N1 - 1) * N0]);
    }

    printf("---\n");
}