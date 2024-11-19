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

    register Real2 Rp0 = {0.0, 0.0};
    register Real2 Rp1 = {0.0, 0.0};
    register Real2 Rp2 = {0.0, 0.0};
    register Real2 Rp3 = {0.0, 0.0};

    for(; N1 + 3 < N0; N1 += 4) {
        Rp0 = vaddq_f64(Rp0, vmulq_f64(Cv0[N1], Cv0[N1]));
        Rp1 = vaddq_f64(Rp1, vmulq_f64(Cv0[N1 + 1], Cv0[N1 + 1]));
        Rp2 = vaddq_f64(Rp2, vmulq_f64(Cv0[N1 + 2], Cv0[N1 + 2]));
        Rp3 = vaddq_f64(Rp3, vmulq_f64(Cv0[N1 + 3], Cv0[N1 + 3]));
    }

    for(; N1 < N0; ++N1)
        Rp0 = vaddq_f64(Rp0, vmulq_f64(Cv0[N1], Cv0[N1]));

    return sqrt(vaddvq_f64(Rp0) + vaddvq_f64(Rp1) + vaddvq_f64(Rp2) + vaddvq_f64(Rp3));
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