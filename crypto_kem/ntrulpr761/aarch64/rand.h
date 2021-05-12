#ifndef RAND_H
#define RAND_H

#include "kem.h"
#include "params.h"
#include "randombytes.h"
#include "sort.h"
#include "type.h"

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

    randombytes((unsigned char *)L, sizeof L);
    sort_to_short(dst, L);
}

#endif
