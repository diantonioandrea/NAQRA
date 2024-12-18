/**
 * @file Vector.h
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Complex vectors.
 * @date 2024-11-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef NAQRA_VECTOR_H
#define NAQRA_VECTOR_H

// Complex numbers.
#include "./Complex.h"

// Copy.

void Cp_CvtCvN_0(Complex*, const Complex*, Natural);

// Complex-Real arithmetic.

void M_CvR_0(Complex* , const Real, const Natural);
void D_CvR_0(Complex* , const Real, const Natural);

// Dot product.

Complex Dot_CrvCcvN_C(const Complex*, const Complex*, const Natural);

// Norms.

Real N2_CvN_R(const Complex*, const Natural);
void Nz2_CvN_0(Complex*, const Natural);

// Output.

void Pn_CrvN_0(const Complex*, const Natural);
void Pn_CcvN_0(const Complex*, const Natural);

#endif