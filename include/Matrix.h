/**
 * @file Matrix.h
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Complex matrices, column-major storage.
 * @date 2024-11-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef NAQRA_MATRIX_H
#define NAQRA_MATRIX_H

// Vectors.
#include "./Vector.h"

// Householder products.

void Hsl_CqtCvNN_0(Complex*, const Complex*, const Natural, const Natural);
void Hsr_CqtCvNN_0(Complex*, const Complex*, const Natural, const Natural);

// Givens products.

void Gvl_ChsnqtCCNN_0(Complex*, const Complex, const Complex, const Natural, const Natural);
void Gvrhr_ChsnqtCCNN_0(Complex*, const Complex, const Complex, const Natural, const Natural);

// Hessenberg form.

void Hsn_CqtN_0(Complex*, const Natural);

// QR algorithm.

void Eig_ChsnqtN_0(Complex *, const Natural);

// Output.

void Pn_CmNN_0(const Complex*, const Natural, const Natural);

#endif