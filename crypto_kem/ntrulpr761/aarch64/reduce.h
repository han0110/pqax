#ifndef REDUCE_H
#define REDUCE_H

#include <stdint.h>

#include "params.h"

// TODO: optimize
// unsigned_barrett_q does barrett reduction twice to ensure -q/2 < output < q/2
int16_t unsigned_barrett_q(uint32_t src) {
    const int k = 31;
    const uint64_t m = ((uint64_t)1 << k) / NTRU_LPRIME_Q;

    src -= (uint32_t)(NTRU_LPRIME_Q * ((src * m) >> k));
    src -= (uint32_t)(NTRU_LPRIME_Q * ((src * m) >> k));
    src -= NTRU_LPRIME_Q;
    return src + (~(src >> k) + 1) & (uint32_t)NTRU_LPRIME_Q;
}

// barrett_q reduces to -q/2 to q/2 for -14000000 < src < 14000000
// if q in {4591, 4621, 5167}
static int16_t barrett_q(int32_t src) {
    const int k1 = 18, k2 = 27;
    const int32_t m1 = ((1 << k1) + NTRU_LPRIME_Q / 2) / NTRU_LPRIME_Q,
                  m2 = ((1 << k2) + NTRU_LPRIME_Q / 2) / NTRU_LPRIME_Q;

    src -= NTRU_LPRIME_Q * ((m1 * src) >> k1);
    return src - NTRU_LPRIME_Q * ((m2 * src + (1 << (k2 - 1))) >> k2);
}

#endif
