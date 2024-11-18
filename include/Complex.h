/**
 * @file Complex.h
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Neon complex numbers.
 * @date 2024-11-18
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef NAQRA_COMPLEX_H
#define NAQRA_COMPLEX_H

// Neon.
#include <arm_neon.h>

// Complex number.
typedef float64_t Real;
typedef float64x2_t Complex;

// Complex "constructors".
inline Complex C_R_C(Real R0) {
    Complex C0;

    C0 = vsetq_lane_f64(R0, C0, 0);
    C0 = vsetq_lane_f64(0.0, C0, 1);

    return C0;
}

inline Complex C_RR_C(Real R0, Real R1) {
    Complex C0;

    C0 = vsetq_lane_f64(R0, C0, 0);
    C0 = vsetq_lane_f64(R1, C0, 1);

    return C0;
}

// Complex-Complex arithmetic.
inline Complex A_CC_C(Complex C0, Complex C1) { return vaddq_f64(C0, C1); }
inline Complex S_CC_C(Complex C0, Complex C1) { return vsubq_f64(C0, C1); }
inline Complex M_CC_C(Complex C0, Complex C1) {}
inline Complex D_CC_C(Complex C0, Complex C1) {}

// Complex-Real arithmetic.
inline Complex M_CR_C(Complex C0, Real C1) {}
inline Complex D_CR_C(Complex C0, Real C1) {}

// Output.
inline void P_C_0(Complex C0) {}

#endif