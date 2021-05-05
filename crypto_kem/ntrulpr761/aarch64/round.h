#ifndef ROUND_H
#define ROUND_H

#include <arm_neon.h>

#include "params.h"

#define THREE 3
#define BAR_3_K 17

static int16_t round3(int16_t src) {
    const int32_t r = (int32_t)1 << (BAR_3_K - 1);
    const int32_t m = (r * 2 + THREE / 2) / THREE;

    // 3 * (⌈in * m⌋ >> 17)
    int32_t quotient = (src * m + r) >> BAR_3_K;
    return quotient * THREE;
}

static void round3_x4(int16_t* dst, const int16_t* src) {
    const int32_t r = (int32_t)1 << (BAR_3_K - 1);
    const int32_t m = (r * 2 + THREE / 2) / THREE;
    const int32x4_t mx4 = {m, m, m, m};

    int16x4_t in16x4 = vld1_s16(src);
    int32x4_t in32x4 = vmovl_s16(in16x4);

    // 3 * (⌈in * m⌋ >> 17)
    int32x4_t tmp = vmulq_s32(in32x4, mx4);
    tmp = vrshrq_n_s32(tmp, BAR_3_K);
    tmp = vmulq_n_s32(tmp, THREE);

    int16x4_t out16x4 = vmovn_s32(tmp);
    vst1_s16(dst, out16x4);
}

void round3_xp(int16_t* dst, const int16_t* src) {
    for (int i = 0; i < NTRU_LPRIME_P / 4 * 4; i += 4) {
        round3_x4(dst + i, src + i);
    }
    dst[NTRU_LPRIME_P - 1] = round3(src[NTRU_LPRIME_P - 1]);
}

#endif
