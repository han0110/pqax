#ifndef POLY_H
#define POLY_H

#include <stdint.h>

#include "params.h"
#include "reduce.h"
#include "type.h"

void mul_small_schoolbook(Fq *h, const Fq *f, const small *g) {
    Fq f_reduced[NTRU_LPRIME_P];
    Fq fg[2 * NTRU_LPRIME_P - 1];
    int32_t result;
    int i, j;

    for (i = 0; i < NTRU_LPRIME_P; ++i) {
        f_reduced[i] = barrett_q(f[i]);
    }
    for (i = 0; i < NTRU_LPRIME_P; ++i) {
        result = 0;
        for (j = 0; j <= i; ++j)
            result = result + f_reduced[j] * (int32_t)g[i - j];
        fg[i] = barrett_q(result);
    }
    for (i = NTRU_LPRIME_P; i < 2 * NTRU_LPRIME_P - 1; ++i) {
        result = 0;
        for (j = i - NTRU_LPRIME_P + 1; j < NTRU_LPRIME_P; ++j)
            result = result + f_reduced[j] * (int32_t)g[i - j];
        fg[i] = barrett_q(result);
    }

    for (i = 2 * NTRU_LPRIME_P - 2; i >= NTRU_LPRIME_P; --i) {
        fg[i - NTRU_LPRIME_P] = barrett_q(fg[i - NTRU_LPRIME_P] + fg[i]);
        fg[i - NTRU_LPRIME_P + 1] =
            barrett_q(fg[i - NTRU_LPRIME_P + 1] + fg[i]);
    }

    for (i = 0; i < NTRU_LPRIME_P; ++i) h[i] = fg[i];
}

// TODO: use good's trick
void mul_small(Fq *h, const Fq *f, const small *g) {
    mul_small_schoolbook(h, f, g);
}

#endif
