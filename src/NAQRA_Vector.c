/**
 * @file NAQRA_Vector.c
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief include/Vector.h implementation.
 * @date 2024-11-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <math.h>
#include "../include/Vector.h"

// Dot product.

/**
 * @brief Dot [D].
 * 
 * @param Crv0 Complex Row Vector [Crv].
 * @param Ccv0 Complex Column Vector [Ccv].
 * @param N0 Entries [N].
 * @return Complex Complex [C].
 */
Complex Dot_CrvCcvN_C(const Complex* Crv0, const Complex* Ccv0, const Natural N0) {
    register Natural N1 = 0;

    register Complex C0 = {0.0, 0.0};
    register Complex C1 = {0.0, 0.0};
    register Complex C2 = {0.0, 0.0};
    register Complex C3 = {0.0, 0.0};

    for(; N1 + 3 < N0; N1 += 4) {
        C0 = vaddq_f64(C0, M_CCcj_C(Crv0[N1], Ccv0[N1]));
        C1 = vaddq_f64(C1, M_CCcj_C(Crv0[N1 + 1], Ccv0[N1 + 1]));
        C2 = vaddq_f64(C2, M_CCcj_C(Crv0[N1 + 2], Ccv0[N1 + 2]));
        C3 = vaddq_f64(C3, M_CCcj_C(Crv0[N1 + 3], Ccv0[N1 + 3]));
    }

    for(; N1 < N0; ++N1)
        C0 = vaddq_f64(C0, M_CCcj_C(Crv0[N1], Ccv0[N1]));

    return vaddq_f64(vaddq_f64(C0, C1), vaddq_f64(C2, C3));
}

// Norms.

/**
 * @brief Norm 2 [N2].
 * 
 * @param Cv0 Complex Vector [Cv].
 * @param N0 Entries [N].
 * @return Real Real [R].
 */
Real N2_CvN_R(const Complex* Cv0, const Natural N0) {
    register Natural N1 = 0;

    register Real2 Rt0 = {0.0, 0.0};
    register Real2 Rt1 = {0.0, 0.0};
    register Real2 Rt2 = {0.0, 0.0};
    register Real2 Rt3 = {0.0, 0.0};

    for(; N1 + 3 < N0; N1 += 4) {
        Rt0 = vaddq_f64(Rt0, vmulq_f64(Cv0[N1], Cv0[N1]));
        Rt1 = vaddq_f64(Rt1, vmulq_f64(Cv0[N1 + 1], Cv0[N1 + 1]));
        Rt2 = vaddq_f64(Rt2, vmulq_f64(Cv0[N1 + 2], Cv0[N1 + 2]));
        Rt3 = vaddq_f64(Rt3, vmulq_f64(Cv0[N1 + 3], Cv0[N1 + 3]));
    }

    for(; N1 < N0; ++N1)
        Rt0 = vaddq_f64(Rt0, vmulq_f64(Cv0[N1], Cv0[N1]));

    return sqrt(vaddvq_f64(Rt0) + vaddvq_f64(Rt1) + vaddvq_f64(Rt2) + vaddvq_f64(Rt3));
}

/**
 * @brief Normalize 2 [Nz2].
 * 
 * @param Cvt0 Cv0 Complex Vector [Cv], Target [t].
 * @param N0 Entries [N].
 */
void Nz2_CvN_0(Complex* Cvt0, const Natural N0) {
    register Natural N1 = 0;

    register Real2 Rt1 = {0.0, 0.0};
    register Real2 Rt2 = {0.0, 0.0};
    register Real2 Rt3 = {0.0, 0.0};
    register Real2 Rt0 = {0.0, 0.0};

    for(; N1 + 3 < N0; N1 += 4) {
        Rt0 = vaddq_f64(Rt0, vmulq_f64(Cvt0[N1], Cvt0[N1]));
        Rt1 = vaddq_f64(Rt1, vmulq_f64(Cvt0[N1 + 1], Cvt0[N1 + 1]));
        Rt2 = vaddq_f64(Rt2, vmulq_f64(Cvt0[N1 + 2], Cvt0[N1 + 2]));
        Rt3 = vaddq_f64(Rt3, vmulq_f64(Cvt0[N1 + 3], Cvt0[N1 + 3]));
    }

    for(; N1 < N0; ++N1)
        Rt0 = vaddq_f64(Rt0, vmulq_f64(Cvt0[N1], Cvt0[N1]));

    const register Real2 Rt4 = vdupq_n_f64(vaddvq_f64(Rt0) + vaddvq_f64(Rt1) + vaddvq_f64(Rt2) + vaddvq_f64(Rt3));

    for(N1 = 0; N1 + 3 < N0; N1 += 4) {
        Cvt0[N1] = vdivq_f64(Cvt0[N1], Rt4);
        Cvt0[N1 + 1] = vdivq_f64(Cvt0[N1 + 1], Rt4);
        Cvt0[N1 + 2] = vdivq_f64(Cvt0[N1 + 2], Rt4);
        Cvt0[N1 + 3] = vdivq_f64(Cvt0[N1 + 3], Rt4);
    }

    for(; N1 < N0; ++N1)
        Cvt0[N1] = vdivq_f64(Cvt0[N1], Rt4);
}

// Output.

/**
 * @brief Print with new line [Pn].
 * 
 * @param Crv0 Complex Row Vector [Crv].
 * @param N0 Entries [N].
 */
void Pn_CrvN_0(const Complex* Crv0, const Natural N0) {
    register Natural N1 = 0;

    printf("--- Row Vector\n");

    for(; N1 < N0; ++N1)
        P_C_0(Crv0[N1]);

    printf("\n---\n");
}

/**
 * @brief Print with new line [Pn].
 * 
 * @param Ccv0 Complex Column Vector [Ccv].
 * @param N0 Entries [N].
 */
void Pn_CcvN_0(const Complex* Crc0, const Natural N0) {
    register Natural N1 = 0;

    printf("--- Column Vector\n");

    for(; N1 < N0; ++N1)
        Pn_C_0(Crc0[N1]);

    printf("---\n");
}