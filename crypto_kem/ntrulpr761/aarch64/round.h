#ifndef ROUND_H
#define ROUND_H

#include <arm_neon.h>

#include "params.h"

static int16_t round3(int16_t src) {
    const int k = 15;
    const int16_t r = (int16_t)1 << (k - 1);
    const int16_t m = (r * 2 + 1) / 3;

    return 3 * ((src * m + r) >> k);
}

static void round3_x8(int16_t* dst, const int16_t* src) {
    const int k = 15;
    const int16_t r = (int16_t)1 << (k - 1);
    const int16_t m = (r * 2 + 1) / 3;

    int16x8_t mx8 = vdupq_n_s16(m);
    int16x8_t x8 = vld1q_s16(src);

    x8 = vqrdmulhq_s16(x8, mx8);
    x8 = vmulq_n_s16(x8, 3);

    vst1q_s16(dst, x8);
}

void round3_xp(int16_t* dst, const int16_t* src) {
    for (int i = 0; i < NTRU_LPRIME_P / 8 * 8; i += 8) {
        round3_x8(dst + i, src + i);
    }
    dst[NTRU_LPRIME_P - 1] = round3(src[NTRU_LPRIME_P - 1]);
}

#endif
