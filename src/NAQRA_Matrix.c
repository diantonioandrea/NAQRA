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
    register Natural N2 = 0, N3;
    const register Natural N4 = N0 - N1;

    for(; N2 < N0; ++N2) {
        register Complex C0 = {0.0, 0.0};

        for(N3 = N4; N3 < N0; ++N3)
            C0 = A_CC_C(C0, M_CcjC_C(Cv0[N3], Cqt0[N2 * N0 + N3]));

        C0 = M_CR_C(C0, 2.0);

        for(N3 = N4; N3 < N0; ++N3)
            Cqt0[N2 * N0 + N3] = S_CC_C(Cqt0[N2 * N0 + N3], M_CC_C(Cv0[N3], C0));
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
 * @param Cqt0 Cqt0 Complex Sqaure Matrix [Cq], Target [t].
 * @param C0 Complex Number [C].
 * @param C1 Complex Number [C].
 * @param N0 Rows and Columns [N].
 * @param N1 Index [N].
 * @param N2 Offset [N].
 */
void Gvl_CqtCCNNN_0(Complex* Cqt0, const Complex C0, const Complex C1, const Natural N0, const Natural N1, const Natural N2) {
    register Natural N3 = N1;
    register Natural N4 = N1 + N2;

    for(; N3 < N0 * N0; N3 += N0, N4 += N0) {
        const register Complex C2 = Cqt0[N3];
        const register Complex C3 = Cqt0[N4];

        Cqt0[N3] = A_CC_C(M_CC_C(C0, C2), M_CC_C(C1, C3));
        Cqt0[N4] = S_CC_C(M_CC_C(C0, C3), M_CC_C(C1, C2));
    }
}

/**
 * @brief Transposed Givens Right [Gvrtr].
 * 
 * @param Cqt0 Cqt0 Complex Sqaure Matrix [Cq], Target [t].
 * @param C0 Complex Number [C].
 * @param C1 Complex Number [C].
 * @param N0 Rows and Columns [N].
 * @param N1 Index [N].
 * @param N2 Offset [N].
 */
void Gvrtr_CqtCCNNN_0(Complex* Cqt0, const Complex C0, const Complex C1, const Natural N0, const Natural N1, const Natural N2) {
    register Natural N3 = N1 * N0;
    register Natural N4 = (N1 + N2) * N0;

    for(; N3 < (N1 + 1) * N0; ++N3, ++N4) {
        const register Complex C2 = Cqt0[N3];
        const register Complex C3 = Cqt0[N4];

        Cqt0[N3] = A_CC_C(M_CC_C(C2, C0), M_CC_C(C3, C1));
        Cqt0[N4] = S_CC_C(M_CC_C(C3, C0), M_CC_C(C2, C1));
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

    for(; N1 < N0 - 2; ++N1) {
        
        // Householder vector.

        const register Natural N2 = N0 - N1 - 1; // Entries.

        Cp_CvtCvN_0(Cv0 + N1 + 1, Cqt0 + N1 * (N0 + 1) + 1, N2); // Copy.
        Cv0[N1 + 1] = S_CC_C(Cv0[N1 + 1], M_CR_C(Nzd2_C_C(Cv0[N1 + 1]), N2_CvN_R(Cv0 + N1 + 1, N2))); // Direction.
        Nz2_CvN_0(Cv0 + N1 + 1, N2); // Normalization.

        // Householder products.

        Hsl_CqtCvNN_0(Cqt0, Cv0, N0, N2);
        Hsr_CqtCvNN_0(Cqt0, Cv0, N0, N2);
    }

    free(Cv0);
}

// QR algorithm.

/**
 * @brief Eigenvalues [Eig].
 * 
 * @param Chsnqt0 Complex Hessenberg Square Matrix [Chsnq], Target [t].
 * @param N0 Rows and Columns [N].
 */
void Eig_ChsnqtN_0(Complex *Chsnqt0, const Natural N0) {
    register Natural N1 = 0, N2;
    register Complex* Cv0 = (Complex*) calloc(2 * (N0 - 1), sizeof(Complex));

    for(; N1 < IT0; ++N1) {
        for(N2 = 0; N2 < N0 - 1; ++N2) {
            Cv0[2 * N2] = Chsnqt0[N2 * (N0 + 1)];
            Cv0[2 * N2 + 1] = Chsnqt0[N2 * (N0 + 1) + 1];

            Nz2_CvN_0(Cv0 + 2 * N2, 2);

            Gvl_CqtCCNNN_0(Chsnqt0, Cv0[2 * N2], Cv0[2 * N2 + 1], N0, N2, 1);
        }

        for(N2 = 0; N2 < N0 - 1; ++N2)
            Gvrtr_CqtCCNNN_0(Chsnqt0, Cv0[2 * N2], Cv0[2 * N2 + 1], N0, N2, 1);
    }

    free(Cv0);
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