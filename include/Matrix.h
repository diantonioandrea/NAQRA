/**
 * @file Matrix.h
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Complex matrices.
 * @date 2024-11-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef NAQRA_MATRIX_H
#define NAQRA_MATRIX_H

// Complex numbers.
#include "./Complex.h"

// Hessenberg form.

void Hsn_CqtN_0(Complex*, const Natural);

// Output.

void Pn_CmNN_0(const Complex*, const Natural, const Natural);

#endif