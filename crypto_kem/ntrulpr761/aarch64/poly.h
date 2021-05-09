#ifndef POLY_H
#define POLY_H

#include <stdint.h>

#include "params.h"
#include "reduce.h"
#include "type.h"

void mul_small_schoolbook(Fq *h, const Fq *f, const small *g) {
    Fq fg[NTRU_LPRIME_P + NTRU_LPRIME_P - 1];
    Fq result;
    int i, j;

    for (i = 0; i < NTRU_LPRIME_P; ++i) {
        result = 0;
        for (j = 0; j <= i; ++j)
            result = barrett_reduce_q(result + f[j] * (int32_t)g[i - j]);
        fg[i] = result;
    }
    for (i = NTRU_LPRIME_P; i < NTRU_LPRIME_P + NTRU_LPRIME_P - 1; ++i) {
        result = 0;
        for (j = i - NTRU_LPRIME_P + 1; j < NTRU_LPRIME_P; ++j)
            result = barrett_reduce_q(result + f[j] * (int32_t)g[i - j]);
        fg[i] = result;
    }

    for (i = NTRU_LPRIME_P + NTRU_LPRIME_P - 2; i >= NTRU_LPRIME_P; --i) {
        fg[i - NTRU_LPRIME_P] = barrett_reduce_q(fg[i - NTRU_LPRIME_P] + fg[i]);
        fg[i - NTRU_LPRIME_P + 1] =
            barrett_reduce_q(fg[i - NTRU_LPRIME_P + 1] + fg[i]);
    }

    for (i = 0; i < NTRU_LPRIME_P; ++i) h[i] = fg[i];
}

void mul_small(Fq *h, const Fq *f, const small *g) {
    mul_small_schoolbook(h, f, g);
    // TODO:
}

#endif
