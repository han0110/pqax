#ifndef ENCODE_H
#define ENCODE_H

#include "kem.h"
#include "params.h"
#include "type.h"

#define NTRU_LPRIME_Q_HALF_FLOOR (NTRU_LPRIME_Q - 1) / 2
#define NTRU_LPRIME_Q_THIRD_CEIL (NTRU_LPRIME_Q + 2) / 3
#define UBAR_3_K 15
#define UBAR_3_M ((1 << UBAR_3_K) + 1) / 3

#define for_range(n, task)        \
    do {                          \
        for (i = 0; i < n; ++i) { \
            task                  \
        }                         \
    } while (0)

#define encode_round_first()                                                 \
    r0 = ((R[2 * i] + NTRU_LPRIME_Q_HALF_FLOOR) * UBAR_3_M) >> UBAR_3_K;     \
    r1 = ((R[2 * i + 1] + NTRU_LPRIME_Q_HALF_FLOOR) * UBAR_3_M) >> UBAR_3_K; \
    r2 = r0 + r1 * (uint32_t)NTRU_LPRIME_Q_THIRD_CEIL;                       \
    *dst++ = r2;                                                             \
    r2 >>= 8;                                                                \
    R[i] = r2;

#define encode_round_small_r(r) \
    r0 = R[2 * i];              \
    r1 = R[2 * i + 1];          \
    r2 = r0 + r1 * (uint32_t)r; \
    *dst++ = r2;                \
    r2 >>= 8;                   \
    R[i] = r2;

#define encode_round_large_r(r) \
    r0 = R[2 * i];              \
    r1 = R[2 * i + 1];          \
    r2 = r0 + r1 * (uint32_t)r; \
    *dst++ = r2;                \
    r2 >>= 8;                   \
    *dst++ = r2;                \
    r2 >>= 8;                   \
    R[i] = r2;

#define encode_round_last() \
    r0 = R[0];              \
    *dst++ = r0;            \
    r0 >>= 8;               \
    *dst++ = r0;            \
    r0 >>= 8;

// TODO: optmize with neon
void encode_round(uint8_t* dst, Fq* R) {
    int i;
    uint16_t r0, r1;
    uint32_t r2;

    // len = 761
    for_range(380, encode_round_first());
    R[380] = ((R[760] + NTRU_LPRIME_Q_HALF_FLOOR) * UBAR_3_M) >> UBAR_3_K;
    // len = 381; 9157 = (1531 * 1531 + 255) >> 8
    for_range(190, encode_round_large_r(9157));
    R[190] = R[380];
    // len = 191; 1280 = (((9157 * 9157 + 255) >> 8) + 255) >> 8
    for_range(95, encode_round_small_r(1280));
    R[95] = R[190];
    // len = 96; 6400 = (1280 * 1280 + 255) >> 8
    for_range(48, encode_round_large_r(6400));
    // len = 48; 625 = (((6400 * 6400 + 255) >> 8) + 255) >> 8
    for_range(24, encode_round_small_r(625));
    // len = 24; 1526 = (625 * 625 + 255) >> 8
    for_range(12, encode_round_small_r(1526));
    // len = 12; 9097 = (1526 * 1526 + 255) >> 8
    for_range(6, encode_round_large_r(9097));
    // len = 6; 1263 = (((9097 * 9097 + 255) >> 8) + 255) >> 8
    for_range(3, encode_round_small_r(1263));
    // len = 3; 6232 = (1263 * 1263 + 255) >> 8
    i = 0;
    encode_round_large_r(6232);
    R[1] = R[2];
    // len = 2; 593 = (((6232 * 6232 + 255) >> 8) + 255) >> 8
    i = 0;
    encode_round_small_r(593);
    // len = 1
    encode_round_last();
}

void encode_small(uint8_t* dst, small* src) {
    small x;

    for (int i = 0; i < NTRU_LPRIME_P / 4; ++i) {
        x = *src++ + 1;
        x += (*src++ + 1) << 2;
        x += (*src++ + 1) << 4;
        x += (*src++ + 1) << 6;
        *dst++ = x;
    }
    x = *src++ + 1;
    *dst++ = x;
}

void encode_input(unsigned char* dst, const uint8_t* r) {
    int i;

    for (i = 0; i < CRYPTO_BYTES; ++i) dst[i] = 0;
    for (i = 0; i < 8 * CRYPTO_BYTES; ++i) dst[i >> 3] |= r[i] << (i & 7);
}

void encode_top(unsigned char* s, const int8_t* T) {
    for (int i = 0; i < 128; ++i) s[i] = T[2 * i] + (T[2 * i + 1] << 4);
}

#endif
