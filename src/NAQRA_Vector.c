/**
 * @file NAQRA_Vector.c
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief include/Vector.h implementation.
 * @date 2024-11-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "../include/Vector.h"

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