#ifndef DECODE_H
#define DECODE_H

#include <stdint.h>

#include "kem.h"
#include "params.h"
#include "type.h"

#define udivmod_x(N)                                                  \
    static void udivmod_##N(uint32_t *q, uint16_t *r, uint32_t src) { \
        const uint64_t M = ((uint32_t)1 << 31) / N;                   \
        uint32_t q2;                                                  \
        uint32_t mask;                                                \
                                                                      \
        *q = ((uint64_t)src * M) >> 31;                               \
        src -= *q * N;                                                \
        q2 = ((uint64_t)src * M) >> 31;                               \
        src -= q2 * N;                                                \
        *q += q2;                                                     \
        src -= N;                                                     \
        *q += 1;                                                      \
        mask = -(src >> 31);                                          \
        src += mask & N;                                              \
        *q += mask;                                                   \
        *r = src;                                                     \
    }

udivmod_x(150);
udivmod_x(304);
udivmod_x(367);
udivmod_x(593);
udivmod_x(625);
udivmod_x(1263);
udivmod_x(1280);
udivmod_x(1500);
udivmod_x(1526);
udivmod_x(1531);
udivmod_x(2188);
udivmod_x(3475);
udivmod_x(6232);
udivmod_x(6400);
udivmod_x(9097);
udivmod_x(9157);

#define decode_round_medium(n, divmod_m, divmod_m_last, dst, dst_last) \
    do {                                                               \
        src -= n;                                                      \
        for (i = 0; i < n - 1; ++i) {                                  \
            divmod_m(&r, dst++, *(src + i) + 256 * *(dst_last + i));   \
            divmod_m(&r, dst++, r);                                    \
        }                                                              \
        divmod_m(&r, dst++, *(src + i) + 256 * *(dst_last + i));       \
        divmod_m_last(&r, dst++, r);                                   \
        dst -= 2 * n;                                                  \
    } while (0)

#define decode_round_large(n, divmod_m, divmod_m_last, dst, dst_last)          \
    do {                                                                       \
        src -= 2 * n;                                                          \
        for (i = 0; i < n - 1; ++i) {                                          \
            divmod_m(&r, dst++,                                                \
                     *src++ + (*src++ << 8) + 65536 * *(dst_last + i));        \
            divmod_m(&r, dst++, r);                                            \
        }                                                                      \
        divmod_m(&r, dst++, *src++ + (*src++ << 8) + 65536 * *(dst_last + i)); \
        divmod_m_last(&r, dst++, r);                                           \
        src -= 2 * n;                                                          \
        dst -= 2 * n;                                                          \
    } while (0)

void decode_round(int16_t *dst, const unsigned char *src) {
    int i;
    int16_t *dst_hi = &dst[380];
    uint32_t r;

    // len = 1
    src += NTRU_LPRIME_ROUND_ENCODED_BYTES - 2;
    udivmod_3475(&r, dst, *src + (*(src + 1) << 8));
    // len = 2
    decode_round_medium(1, udivmod_593, udivmod_1500, dst_hi, dst);
    // len = 3
    decode_round_large(1, udivmod_6232, udivmod_6232, dst, dst_hi);
    *(dst + 2) = *(dst_hi + 1);
    // len = 6
    decode_round_medium(3, udivmod_1263, udivmod_304, dst_hi, dst);
    // len = 12
    decode_round_large(6, udivmod_9097, udivmod_2188, dst, dst_hi);
    // len = 24
    decode_round_medium(12, udivmod_1526, udivmod_367, dst_hi, dst);
    // len = 48
    decode_round_medium(24, udivmod_625, udivmod_150, dst, dst_hi);
    // len = 96
    decode_round_large(48, udivmod_6400, udivmod_1531, dst_hi, dst);
    // len = 191
    decode_round_medium(95, udivmod_1280, udivmod_1280, dst, dst_hi);
    *(dst + 190) = *(dst_hi + 95);
    // len = 381
    decode_round_large(190, udivmod_9157, udivmod_9157, dst_hi, dst);
    *(dst_hi + 380) = *(dst + 190);
    // len = 761
    decode_round_medium(380, udivmod_1531, udivmod_1531, dst, dst_hi);
    *(dst + 760) = *(dst_hi + 380);

    for (i = 0; i < NTRU_LPRIME_P; ++i)
        dst[i] = dst[i] * 3 - NTRU_LPRIME_Q_HALF_FLOOR;
}

void decode_small(small *dst, const unsigned char *src) {
    unsigned char x;

    for (int i = 0; i < NTRU_LPRIME_P / 4; ++i) {
        x = *src++;
        *dst++ = (small)(x & 3) - 1;
        x >>= 2;
        *dst++ = (small)(x & 3) - 1;
        x >>= 2;
        *dst++ = (small)(x & 3) - 1;
        x >>= 2;
        *dst++ = (small)(x & 3) - 1;
    }
    x = *src++;
    *dst++ = (small)(x & 3) - 1;
}

void decode_input(uint8_t *dst, const unsigned char *src) {
    for (int i = 0; i < 8 * NTRU_LPRIME_INPUT_BYTES; ++i)
        dst[i] = 1 & (src[i >> 3] >> (i & 7));
}

void decode_top(int8_t *T, const unsigned char *src) {
    for (int i = 0; i < 4 * NTRU_LPRIME_INPUT_BYTES; ++i) {
        T[2 * i] = src[i] & 15;
        T[2 * i + 1] = src[i] >> 4;
    }
}

#endif
