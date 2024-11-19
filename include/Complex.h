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

// Types.
#include <stddef.h>

// Output.
#include <stdio.h>

// Neon.
#include <arm_neon.h>


typedef size_t Natural; // Natural numbers.
typedef ptrdiff_t Integer; // Integer numbers.
typedef float64_t Real; // Real numbers.
typedef float64x2_t Complex; // Complex numbers.

typedef float64x1_t Real1; // Real numbers (singlet).
typedef float64x2_t Real2; // Real numbers (double).


// Complex "constructors".

static inline Complex C_R_C(const Real R0) { Complex C0 = {R0, 0.0}; return C0; }
static inline Complex C_RR_C(const Real R0, const Real R1) { Complex C0 = {R0, R1}; return C0; }

// Complex-Complex arithmetic.

static inline Complex A_CC_C(const Complex C0, const Complex C1) { return vaddq_f64(C0, C1); }
static inline Complex S_CC_C(const Complex C0, const Complex C1) { return vsubq_f64(C0, C1); }

static inline Complex M_CC_C(const Complex C0, const Complex C1) {
    const register Real R0 = vgetq_lane_f64(C0, 0), R1 = vgetq_lane_f64(C0, 1); // C0.
    const register Real R2 = vgetq_lane_f64(C1, 0), R3 = vgetq_lane_f64(C1, 1); // C1.

    Complex C2 = {R0 * R2 - R1 * R3, R0 * R3 + R1 * R2}; return C2;
}

static inline Complex M_CCcj_C(const Complex C0, const Complex C1) {
    const register Real R0 = vgetq_lane_f64(C0, 0), R1 = vgetq_lane_f64(C0, 1); // C0.
    const register Real R2 = vgetq_lane_f64(C1, 0), R3 = vgetq_lane_f64(C1, 1); // C1.

    Complex C2 = {R0 * R2 + R1 * R3, R1 * R2 - R0 * R3}; return C2;
}

static inline Complex D_CC_C(const Complex C0, const Complex C1) {
    const register Real R0 = vgetq_lane_f64(C0, 0), R1 = vgetq_lane_f64(C0, 1); // C0.
    const register Real R2 = vgetq_lane_f64(C1, 0), R3 = vgetq_lane_f64(C1, 1); // C1.
    const register Real R4 = R2 * R2 + R3 * R3; // C1.

    Complex C2 = {(R0 * R2 + R1 * R3) / R4, (R1 * R2 - R0 * R3) / R4}; return C2;
}

// Complex-Real arithmetic.

static inline Complex A_CR_C(const Complex C0, const Real R1) { const register Complex C1 = {R1, 0.0}; return vaddq_f64(C0, C1); }
static inline Complex S_CR_C(const Complex C0, const Real R1) { const register Complex C1 = {R1, 0.0}; return vsubq_f64(C0, C1); }
static inline Complex M_CR_C(const Complex C0, const Real R0) { return vmulq_f64(C0, vdupq_n_f64(R0)); }
static inline Complex D_CR_C(const Complex C0, const Real R0) { return vdivq_f64(C0, vdupq_n_f64(R0)); }

// Complex methods.

static inline Complex Cj_C_C(const Complex C0) { return vsetq_lane_f64(-vgetq_lane_f64(C0, 1), C0, 1); }

// Output.

static inline void P_C_0(const Complex C0) { 
    const register Real R0 = vgetq_lane_f64(C0, 0), R1 = vgetq_lane_f64(C0, 1); // C0.

    if(R0 >= 0.0) printf(" "); printf("%.2e", R0);
    if(R1 >= 0.0) printf("+"); printf("%.2ei ", R1);
}

static inline void Pn_C_0(const Complex C0) { P_C_0(C0); printf("\n"); }

#endif