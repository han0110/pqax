#ifndef RAND_H
#define RAND_H

#include "params.h"
#include "randombytes.h"
#include "sort.h"
#include "type.h"

// TODO: optmize with neon
static uint32_t random_u32() {
    uint8_t c[4];
    randombytes(c, 4);
    return (uint32_t)c[0] + ((uint32_t)c[1] << 8) + ((uint32_t)c[2] << 16) +
           ((uint32_t)c[3] << 24);
}

void random_u8_xn(uint8_t *dst, size_t n) { randombytes(dst, n); }

// TODO: optmize with neon
void random_small_xp(small *dst) {
    uint32_t L[NTRU_LPRIME_P];
    int i;

    for (i = 0; i < NTRU_LPRIME_P; i++) L[i] = random_u32();
    for (i = 0; i < NTRU_LPRIME_W; i++)
        L[i] = L[i] & (uint32_t)-2;  // non-zero part (1 or -1)
    for (; i < NTRU_LPRIME_P; i++)
        L[i] = (L[i] & (uint32_t)-3) | 1;  // zero part

    sort_xp(L);

    for (i = 0; i < NTRU_LPRIME_P; i++) dst[i] = (small)(L[i] & 3) - 1;
}

#endif
