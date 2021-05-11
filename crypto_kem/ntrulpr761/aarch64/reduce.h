#ifndef REDUCE_H
#define REDUCE_H

#include <stdint.h>

#include "params.h"

#define BAR_Q_K 31

// TODO: optimize
// barrett_q does barrett reduction twice to ensure -q/2 < output < q/2
int32_t barrett_q(uint32_t src) {
    // r = 2^31
    const int64_t r = (int64_t)1 << (BAR_Q_K - 1);
    // m = (2^32 + q/2)/q
    const int64_t m = (2 * r + NTRU_LPRIME_Q / 2) / NTRU_LPRIME_Q;

    int32_t quotient = src * m >> BAR_Q_K;
    int32_t signed_src = src - quotient * NTRU_LPRIME_Q;

    quotient = (signed_src * m + r) >> BAR_Q_K;
    return signed_src - quotient * NTRU_LPRIME_Q;
}

// TODO: optimize
// unsigned_barrett_q does barrett reduction twice to ensure -q/2 < output < q/2
int32_t unsigned_barrett_q(int32_t src) {
    // r = 2^31
    const int64_t r = (int64_t)1 << (BAR_Q_K - 1);
    // m = (2^32 + q/2)/q
    const int64_t m = (2 * r + NTRU_LPRIME_Q / 2) / NTRU_LPRIME_Q;

    int32_t quotient = src * m >> BAR_Q_K;
    int32_t signed_src = src - quotient * NTRU_LPRIME_Q;

    quotient = (signed_src * m + r) >> BAR_Q_K;
    return signed_src - quotient * NTRU_LPRIME_Q;
}

#endif
