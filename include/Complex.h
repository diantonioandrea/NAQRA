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

// Base.
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

// Math.
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <math.h>

// Neon.
#include <arm_neon.h>


// Types.

typedef size_t Natural; // Natural numbers.
typedef ptrdiff_t Integer; // Integer numbers.
typedef float64_t Real; // Real numbers.
typedef float64x2_t Complex; // Complex numbers.

typedef float64x1_t Real1; // Real numbers (unit).
typedef float64x2_t Real2; // Real numbers (pair).


// Constants.

#ifndef TOL0

// Almost zero numbers.
#define TOL0 1.0E-12
#endif

#ifndef TOL1

// Definitely zero numbers.
#define TOL1 1.0E-24
#endif

#ifndef ITM0

// Maximum number of iterations.
#define ITM0 1E6
#endif


// Complex "constructors".

/**
 * @brief Construct [C].
 * 
 * @param R0 Real Number [N].
 * @return Complex Complex Number [C].
 */
static inline Complex C_R_C(const Real R0) { Complex C0 = {R0, 0.0}; return C0; }

/**
 * @brief Construct [C].
 * 
 * @param R0 Real Number [N].
 * @param R1 Real Number [N].
 * @return Complex Complex Number [C].
 */
static inline Complex C_RR_C(const Real R0, const Real R1) { Complex C0 = {R0, R1}; return C0; }

// Complex-Complex arithmetic.

/**
 * @brief Add [A].
 * 
 * @param C0 Complex Number [C].
 * @param C1 Complex Number [C].
 * @return Complex Complex Number [C].
 */
static inline Complex A_CC_C(const Complex C0, const Complex C1) { return vaddq_f64(C0, C1); }

/**
 * @brief Subtract [S].
 * 
 * @param C0 Complex Number [C].
 * @param C1 Complex Number [C].
 * @return Complex Complex Number [C].
 */
static inline Complex S_CC_C(const Complex C0, const Complex C1) { return vsubq_f64(C0, C1); }

/**
 * @brief Multiply [M].
 * 
 * @param C0 Complex Number [C].
 * @param C1 Complex Number [C].
 * @return Complex Complex Number [C].
 */
static inline Complex M_CC_C(const Complex C0, const Complex C1) {
    const register Real R0 = vgetq_lane_f64(C0, 0), R1 = vgetq_lane_f64(C0, 1); // C0.
    const register Real R2 = vgetq_lane_f64(C1, 0), R3 = vgetq_lane_f64(C1, 1); // C1.

    Complex C2 = {R0 * R2 - R1 * R3, R0 * R3 + R1 * R2}; return C2;
}

/**
 * @brief Multiply [M].
 * 
 * @param C0 Complex Number [C].
 * @param C1 Conjugate Complex Number [Ccj].
 * @return Complex Complex Number [C].
 */
static inline Complex M_CCcj_C(const Complex C0, const Complex C1) {
    const register Real R0 = vgetq_lane_f64(C0, 0), R1 = vgetq_lane_f64(C0, 1); // C0.
    const register Real R2 = vgetq_lane_f64(C1, 0), R3 = vgetq_lane_f64(C1, 1); // C1.

    Complex C2 = {R0 * R2 + R1 * R3, R1 * R2 - R0 * R3}; return C2;
}

/**
 * @brief Multiply [M].
 * 
 * @param C0 Conjugate Complex Number [Ccj].
 * @param C1 Complex Number [C].
 * @return Complex Complex Number [C].
 */
static inline Complex M_CcjC_C(const Complex C0, const Complex C1) {
    const register Real R0 = vgetq_lane_f64(C0, 0), R1 = vgetq_lane_f64(C0, 1); // C0.
    const register Real R2 = vgetq_lane_f64(C1, 0), R3 = vgetq_lane_f64(C1, 1); // C1.

    Complex C2 = {R0 * R2 + R1 * R3, R0 * R3 - R1 * R2}; return C2;
}

/**
 * @brief Divide [D].
 * 
 * @param C0 Complex Number [C].
 * @param C1 Complex Number [C].
 * @return Complex Complex Number [C].
 */
static inline Complex D_CC_C(const Complex C0, const Complex C1) {
    const register Real R0 = vgetq_lane_f64(C0, 0), R1 = vgetq_lane_f64(C0, 1); // C0.
    const register Real R2 = vgetq_lane_f64(C1, 0), R3 = vgetq_lane_f64(C1, 1); // C1.
    const register Real R4 = R2 * R2 + R3 * R3; // C1.

    Complex C2 = {(R0 * R2 + R1 * R3) / R4, (R1 * R2 - R0 * R3) / R4}; return C2;
}

// Complex-Real arithmetic.

/**
 * @brief Add [A].
 * 
 * @param C0 Complex Number [C].
 * @param R0 Real Number [R].
 * @return Complex Complex Number [C].
 */
static inline Complex A_CR_C(const Complex C0, const Real R1) { const register Complex C1 = {R1, 0.0}; return vaddq_f64(C0, C1); }

/**
 * @brief Subtract [S].
 * 
 * @param C0 Complex Number [C].
 * @param R0 Real Number [R].
 * @return Complex Complex Number [C].
 */
static inline Complex S_CR_C(const Complex C0, const Real R1) { const register Complex C1 = {R1, 0.0}; return vsubq_f64(C0, C1); }

/**
 * @brief Multiply [M].
 * 
 * @param C0 Complex Number [C].
 * @param R0 Real Number [R].
 * @return Complex Complex Number [C].
 */
static inline Complex M_CR_C(const Complex C0, const Real R0) { return vmulq_f64(C0, vdupq_n_f64(R0)); }

/**
 * @brief Divide [D].
 * 
 * @param C0 Complex Number [C].
 * @param R0 Real Number [R].
 * @return Complex Complex Number [C].
 */
static inline Complex D_CR_C(const Complex C0, const Real R0) { return vdivq_f64(C0, vdupq_n_f64(R0)); }

// Complex methods.

/**
 * @brief Conjugate [A].
 * 
 * @param C0 Complex Number [C].
 * @return Complex Complex Number [C].
 */
static inline Complex Cj_C_C(const Complex C0) { return vsetq_lane_f64(-vgetq_lane_f64(C0, 1), C0, 1); }

/**
 * @brief Square [Sq].
 * 
 * @param C0 Complex Number [C].
 * @return Complex Complex Number [C].
 */
static inline Complex Sq_C_C(const Complex C0) {
    const register Real R0 = vgetq_lane_f64(C0, 0), R1 = vgetq_lane_f64(C0, 1); // C0.

    Complex C2 = {R0 * R0 - R1 * R1, 2.0 * R0 * R1}; return C2;
}

/**
 * @brief N-th Roots [Nrt].
 * 
 * @param C0 Complex Number [C].
 * @param N0 Natural Number [N].
 * @return Complex* Complex Vector [Cv].
 */
[[nodiscard]] static inline Complex* Nrt_C_Cv(const Complex C0, const Natural N0) {
    register Natural N1 = 0;
    register Complex* Cv0 = (Complex*) calloc(N0, sizeof(Complex));

    const register Real R0 = vgetq_lane_f64(C0, 0), R1 = vgetq_lane_f64(C0, 1); // C0.
    const register Real R2 = 2.0 * M_PI / (Real) N0;

    register Real R3 = R0 * R0 + R1 * R1; R3 = pow(R3, 0.5 / (Real) N0);
    register Real R4 = atan2(R1, R0); if(R4 < 0.0) R4 += M_PI; R4 /= (Real) N0;

    for(; N1 < N0; ++N1, R4 += R2)
        Cv0[N1] = C_RR_C(R3 * cos(R4), R3 * sin(R4));

    return Cv0;
}

// Norms.

/**
 * @brief Norm2 [N2].
 * 
 * @param C0 Complex Number [N].
 * @return Real Real Number [R].
 */
static inline Real N2_C_R(const Complex C0) { const register Real R0 = vgetq_lane_f64(C0, 0), R1 = vgetq_lane_f64(C0, 1); return sqrt(R0 * R0 + R1 * R1); }

/**
 * @brief Normalized 2 [Nzd2].
 * 
 * @param C0 Complex Number [N].
 * @return Complex Complex Number [C].
 */
static inline Complex Nzd2_C_C(const Complex C0) {
    const register Real R0 = vgetq_lane_f64(C0, 0), R1 = vgetq_lane_f64(C0, 1); 
    return vdivq_f64(C0, vdupq_n_f64(sqrt(R0 * R0 + R1 * R1)));
}

// Output.

/**
 * @brief Print [P].
 * 
 * @param C0 Complex Number [C].
 */
static inline void P_C_0(const Complex C0) { 
    const register Real R0 = vgetq_lane_f64(C0, 0), R1 = vgetq_lane_f64(C0, 1); // C0.

    if(fabs(R0) <= TOL1) { printf("\x1b[2m"); printf(" %.3e", 0.0); printf("\033[0m"); }
    else { if(fabs(R0) <= TOL0) printf("\x1b[2m"); if(R0 >= 0.0) printf(" "); printf("%.3e", R0); printf("\033[0m"); }

    if(fabs(R1) <= TOL1) { printf("\x1b[2m"); printf("+%.3ei ", 0.0); printf("\033[0m"); }
    else { if(fabs(R1) <= TOL0) printf("\x1b[2m"); if(R1 >= 0.0) printf("+"); printf("%.3ei ", R1); printf("\033[0m"); }
}

/**
 * @brief Print with new line [Pn].
 * 
 * @param C0 Complex Number [C].
 */
static inline void Pn_C_0(const Complex C0) { P_C_0(C0); printf("\n"); }

#endif