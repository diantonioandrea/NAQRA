/**
 * @file NAQRA_Matrix.c
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief include/Matrix.h implementation.
 * @date 2024-11-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "../include/Matrix.h"

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

    for(; N2 < N0; ++N2) {
        for(N3 = 0; N3 < N1 - 1; ++N3)
            P_C_0(Cm0[N2 + N3 * N0]);

        Pn_C_0(Cm0[N2 + (N1 - 1) * N0]);
    }
}