/**
 * @file Test_Symmetric.c
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief QR Algorithm test on a symmetric matrix.
 * @date 2024-11-26
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "./Test.h"

int main(int argc, char **argv) {
    if(argc < 4) {
        printf("Usage: %s N (Rows and Columns) R R (Range) * (Output, optional)\n", argv[0]);
        return -1;
    }

    // Size.
    const register Natural N0 = (Natural) atoi(argv[1]);

    // Indices.
    register Natural N1 = 0, N2;

    // Interval.
    const register Real R0 = (Real) atof(argv[2]), R1 = (Real) atof(argv[3]);
    const register Real R2 = R1 - R0;

    // Random value.
    register Real R3;

    // Complex Matrix.
    Complex* Cm0 = (Complex*) malloc((N0 * N0) * sizeof(Complex));

    srand(time(NULL));
    for(; N1 < N0; ++N1) {
        R3 = R0 + R2 * (Real) rand() / RAND_MAX;

        Cm0[N1 * (N0 + 1)] = C_R_C(R3);

        for(N2 = 0; N2 < N1; ++N2) {
            R3 = R0 + R2 * (Real) rand() / RAND_MAX;
            const register Complex C0 = C_R_C(R3);

            Cm0[N2 * N0 + N1] = C0;
            Cm0[N1 * N0 + N2] = C0;
        }
    }

    #ifndef NVERBOSE
    printf("Testing on a symmetric %zu x %zu matrix.\n", N0, N0);
    printf("Coefficients generated in [%.1f, %.1f].\n\n", R0, R1);
    #endif

    #ifndef NVERBOSE
    if(argc > 4)
        Pn_CmNN_0(Cm0, N0, N0);
    #endif

    Hsn_CqtN_0(Cm0, N0); // Hessenberg.

    #ifndef NVERBOSE
    if(argc > 4)
        Pn_CmNN_0(Cm0, N0, N0);
    #endif

    Eig_ChsnqtN_0(Cm0, N0); // Eigenvalues.

    #ifndef NVERBOSE
    if(argc > 4)
        Pn_CmNN_0(Cm0, N0, N0);
    #endif

    free(Cm0);
    return 0;
}