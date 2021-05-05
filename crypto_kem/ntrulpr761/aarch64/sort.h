#ifndef SORT_H
#define SORT_H

#include <stdint.h>

#include "params.h"

static inline void minmax(uint32_t *x, uint32_t *y) {
    uint32_t xv = *x, xy = *y;
    uint32_t mask = -((xy - xv) >> 31);  // 00...00 or 11...11 in bits
    mask &= xv ^ xy;                     // 0 or x ^ y
    *x = xv ^ mask;                      // x or y (x ^ x ^ y)
    *y = xy ^ mask;                      // y or x (y ^ y ^ x)
}

// sort algorithm from ntrulpr761/subroutines/crypto_sort_uint32.c of NTRU Prime
// Round 3 submission
void sort_xp(uint32_t *src) {
    const int n = NTRU_LPRIME_P, top = 512;
    int p, q, i;

    for (p = top; p > 0; p >>= 1) {
        for (i = 0; i < n - p; ++i)
            if (!(i & p)) minmax(src + i, src + i + p);
        for (q = top; q > p; q >>= 1)
            for (i = 0; i < n - q; ++i)
                if (!(i & p)) minmax(src + i + p, src + i + q);
    }
}

#endif
