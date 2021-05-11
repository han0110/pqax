#ifndef SORT_H
#define SORT_H

#include <stdint.h>

#include "params.h"

#define minmax(x, y)                               \
    do {                                           \
        mask = y - x;                              \
        mask ^= (x ^ y) & (mask ^ y ^ 0x80000000); \
        mask = -(mask >> 31);                      \
        mask &= x ^ y;                             \
        x = x ^ mask;                              \
        y = y ^ mask;                              \
    } while (0)

// sort algorithm from ntrulpr761/subroutines/crypto_sort_uint32.c of NTRU
// Prime Round 3 submission
static void sort_xp(uint32_t *src) {
    const int n = NTRU_LPRIME_P, top = 512;
    int p, q, i;
    uint32_t mask;

    for (p = top; p > 0; p >>= 1) {
        for (i = 0; i < n - p; ++i)
            if (!(i & p)) minmax(src[i], src[i + p]);
        for (q = top; q > p; q >>= 1)
            for (i = 0; i < n - q; ++i)
                if (!(i & p)) minmax(src[i + p], src[i + q]);
    }
}

void sort_to_short(small *dst, uint32_t *L) {
    int i;
    for (i = 0; i < NTRU_LPRIME_W; ++i)
        L[i] = L[i] & (uint32_t)-2;  // non-zero part (1 or -1)
    for (i = NTRU_LPRIME_W; i < NTRU_LPRIME_P; ++i)
        L[i] = (L[i] & (uint32_t)-3) | 1;  // zero part

    sort_xp(L);

    for (i = 0; i < NTRU_LPRIME_P; ++i) dst[i] = (L[i] & 3) - 1;
}

#endif
