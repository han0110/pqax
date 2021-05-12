#ifndef RAND_H
#define RAND_H

#include "kem.h"
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

void random_input(int8_t *dst) {
    uint8_t s[CRYPTO_BYTES];

    randombytes(s, CRYPTO_BYTES);
    for (int i = 0; i < 8 * NTRU_LPRIME_INPUT_BYTES; ++i)
        dst[i] = 1 & (s[i >> 3] >> (i & 7));
}

// TODO: optmize with neon
void random_short(small *dst) {
    uint32_t L[NTRU_LPRIME_P];

    for (int i = 0; i < NTRU_LPRIME_P; ++i) L[i] = random_u32();
    sort_to_short(dst, L);
}

#endif
