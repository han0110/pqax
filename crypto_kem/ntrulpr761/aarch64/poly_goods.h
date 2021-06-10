#ifndef POLY_GOODS_H
#define POLY_GOODS_H

#include <arm_neon.h>
#include <stdint.h>

#include "params.h"
#include "poly_goods_macro.h"
#include "reduce.h"
#include "type.h"

#define GOODS_N 512
#define GOODS_N_HALF 256
#define GOODS_P 1536
// 4547 * 1536 + 1
#define GOODS_Q 6984193
#define GOODS_Q_HALF 3492096
#define GOODS_MONT_R_POW 32
// pow(-6984193, -1, 1 << 32)
#define GOODS_MONT_MM_PRIME 2368115199
// pow(6984193, -1, 1 << 32)
#define GOODS_MONT_PM_PRIME 1926852097
// ((1 << 32) % 6984193) - 6984193
#define GOODS_MONT_R_MOD_Q -311399
// (pow(512, -1, 6984193) * 1 << 64) % 6984193
#define GOODS_MONT_R_SQUARE_INV_N_MOD_Q 2770689

// NOTE: python3 script to generate goods_omegas
// q, r, g  = 6984193, -311399, 3991943
// bit_rev = lambda i: int('{:08b}'.format(i)[::-1], 2)
// mod_q = lambda i: i % q if i % q < q / 2 else i % q - q
// goods_omegas = [mod_q(pow(g, bit_rev(i)) * r) for i in range(256)]
int32_t goods_omegas[GOODS_N_HALF] = {
    -311399,  3471433,  2894487,  -2592208, 2462969,  -1536046, -3050359,
    3290424,  -3196524, -2794208, 16826,    1356310,  -3421612, -922778,
    -1259013, -2252520, 804399,   -2052193, 3095387,  2009488,  -1238959,
    -1340819, -2891570, 1431001,  -3367043, -1642889, -1905608, 2340431,
    2867644,  -2977751, -1280757, -3307920, 3454899,  -1292838, -259914,
    2786644,  -2832551, 3486211,  1210248,  -137539,  -751477,  -1407403,
    250674,   -968737,  2000205,  1318306,  2447325,  1605297,  -721754,
    -2223407, -414065,  -679168,  -1100903, -1454521, -2571496, -2397960,
    -924236,  1454474,  2541407,  3381211,  2826879,  -2497269, 204037,
    -94891,   3323907,  859674,   2022528,  -581505,  -105697,  -1792351,
    -2188804, -2324210, 1510310,  -2378504, 1639254,  -3366174, -354333,
    77023,    2500516,  3346795,  1373015,  3283943,  258867,   3244798,
    -395580,  -1855625, -1088322, 181403,   1805696,  -1223284, 2214207,
    -181570,  1228765,  3251180,  2465850,  -848883,  1181496,  -2226884,
    -1451212, -10442,   -1666824, -1209711, -984462,  -3426788, 154192,
    -2907394, -521560,  2929892,  -1025914, -1430178, 2462096,  -2110328,
    2016364,  59786,    -1253072, 1959632,  1911041,  -758518,  176145,
    1773588,  2722679,  1495278,  3253206,  -2040069, 2608536,  -444921,
    -117323,  -1632019, 1517041,  -702754,  -828380,  1389288,  -757948,
    -912063,  -2572353, -689654,  2499044,  2874489,  1448837,  -1823902,
    1415089,  2720529,  960857,   531157,   -474160,  -2586175, -393484,
    -3156896, 2565613,  2968093,  2676380,  -1792652, 772986,   -1379888,
    -263639,  424245,   1918870,  415491,   -297131,  -415074,  -2864466,
    -1107849, 2854985,  1719191,  -1606365, -3247371, 2350237,  516032,
    -2772451, 868382,   2099121,  -2206084, -2235707, -985328,  1939566,
    -1370570, 2278654,  996196,   -2857757, 921439,   -1437016, -200209,
    -3348223, 1321134,  300404,   1152099,  -3057683, -919476,  2059375,
    2216813,  502038,   2212528,  -1947819, 2565716,  -1410605, 1413195,
    -347562,  475050,   830671,   2429855,  1746045,  -2646211, -3164056,
    -1372468, -3451694, -686943,  -3040100, -1569261, -1490851, -3021365,
    1112301,  -3295918, 209553,   -2386487, -2133996, 1369017,  -3419919,
    -2097142, 2999191,  2459030,  460056,   1932437,  -2936036, 2972966,
    -2751330, -1398124, 2030215,  -2204982, 2208840,  -2824097, -1583075,
    -1638985, 1761773,  -788160,  -968429,  1783987,  880438,   -2393162,
    -611880,  -3427476, -2435318, -1538998, -1225208, 3080817,  2570693,
    1333711,  -767312,  -2815827, -818724,  3045125,  1913461,  2257460,
    2566940,  1960976,  -714395,  -1750587};

// NOTE: python3 script to generate goods_inv_omegas
// q, r, g_inv  = 6984193, -311399, -1539202
// bit_rev = lambda i: int('{:08b}'.format(i)[::-1], 2)
// mod_q = lambda i: i % q if i % q < q / 2 else i % q - q
// goods_inv_omegas = [mod_q(pow(g_inv, bit_rev(i)) * r) for i in range(256)]
int32_t goods_inv_omegas[GOODS_N_HALF] = {
    -311399,  -3471433, 2592208,  -2894487, -3290424, 3050359,  1536046,
    -2462969, 2252520,  1259013,  922778,   3421612,  -1356310, -16826,
    2794208,  3196524,  3307920,  1280757,  2977751,  -2867644, -2340431,
    1905608,  1642889,  3367043,  -1431001, 2891570,  1340819,  1238959,
    -2009488, -3095387, 2052193,  -804399,  94891,    -204037,  2497269,
    -2826879, -3381211, -2541407, -1454474, 924236,   2397960,  2571496,
    1454521,  1100903,  679168,   414065,   2223407,  721754,   -1605297,
    -2447325, -1318306, -2000205, 968737,   -250674,  1407403,  751477,
    137539,   -1210248, -3486211, 2832551,  -2786644, 259914,   1292838,
    -3454899, 1632019,  117323,   444921,   -2608536, 2040069,  -3253206,
    -1495278, -2722679, -1773588, -176145,  758518,   -1911041, -1959632,
    1253072,  -59786,   -2016364, 2110328,  -2462096, 1430178,  1025914,
    -2929892, 521560,   2907394,  -154192,  3426788,  984462,   1209711,
    1666824,  10442,    1451212,  2226884,  -1181496, 848883,   -2465850,
    -3251180, -1228765, 181570,   -2214207, 1223284,  -1805696, -181403,
    1088322,  1855625,  395580,   -3244798, -258867,  -3283943, -1373015,
    -3346795, -2500516, -77023,   354333,   3366174,  -1639254, 2378504,
    -1510310, 2324210,  2188804,  1792351,  105697,   581505,   -2022528,
    -859674,  -3323907, 1750587,  714395,   -1960976, -2566940, -2257460,
    -1913461, -3045125, 818724,   2815827,  767312,   -1333711, -2570693,
    -3080817, 1225208,  1538998,  2435318,  3427476,  611880,   2393162,
    -880438,  -1783987, 968429,   788160,   -1761773, 1638985,  1583075,
    2824097,  -2208840, 2204982,  -2030215, 1398124,  2751330,  -2972966,
    2936036,  -1932437, -460056,  -2459030, -2999191, 2097142,  3419919,
    -1369017, 2133996,  2386487,  -209553,  3295918,  -1112301, 3021365,
    1490851,  1569261,  3040100,  686943,   3451694,  1372468,  3164056,
    2646211,  -1746045, -2429855, -830671,  -475050,  347562,   -1413195,
    1410605,  -2565716, 1947819,  -2212528, -502038,  -2216813, -2059375,
    919476,   3057683,  -1152099, -300404,  -1321134, 3348223,  200209,
    1437016,  -921439,  2857757,  -996196,  -2278654, 1370570,  -1939566,
    985328,   2235707,  2206084,  -2099121, -868382,  2772451,  -516032,
    -2350237, 3247371,  1606365,  -1719191, -2854985, 1107849,  2864466,
    415074,   297131,   -415491,  -1918870, -424245,  263639,   1379888,
    -772986,  1792652,  -2676380, -2968093, -2565613, 3156896,  393484,
    2586175,  474160,   -531157,  -960857,  -2720529, -1415089, 1823902,
    -1448837, -2874489, -2499044, 689654,   2572353,  912063,   757948,
    -1389288, 828380,   702754,   -1517041};

// idea from `fqmul` in ~/crypto_kem/kyber768/aarch64/neon_ntt.c
#define mont_mul_neon_n(z, x, y, t)                                         \
    do {                                                                    \
        t.val[0] = (int32x4_t)vmull_n_s32(vget_low_s32(x), y);              \
        t.val[1] = (int32x4_t)vmull_high_n_s32(x, y);                       \
        t.val[2] = vuzp1q_s32(t.val[0], t.val[1]);                          \
        t.val[3] = vuzp2q_s32(t.val[0], t.val[1]);                          \
        t.val[0] = vmulq_n_s32(t.val[2], GOODS_MONT_PM_PRIME);              \
        t.val[1] = (int32x4_t)vmull_n_s32(vget_low_s32(t.val[0]), GOODS_Q); \
        t.val[2] = (int32x4_t)vmull_high_n_s32(t.val[0], GOODS_Q);          \
        t.val[0] = vuzp2q_s32(t.val[1], t.val[2]);                          \
        z = vsubq_s32(t.val[3], t.val[0]);                                  \
    } while (0)

#define mont_mul_neon(z, x, y, t)                                           \
    do {                                                                    \
        t.val[0] = (int32x4_t)vmull_s32(vget_low_s32(x), vget_low_s32(y));  \
        t.val[1] = (int32x4_t)vmull_high_s32(x, y);                         \
        t.val[2] = vuzp1q_s32(t.val[0], t.val[1]);                          \
        t.val[3] = vuzp2q_s32(t.val[0], t.val[1]);                          \
        t.val[0] = vmulq_n_s32(t.val[2], GOODS_MONT_PM_PRIME);              \
        t.val[1] = (int32x4_t)vmull_n_s32(vget_low_s32(t.val[0]), GOODS_Q); \
        t.val[2] = (int32x4_t)vmull_high_n_s32(t.val[0], GOODS_Q);          \
        t.val[0] = vuzp2q_s32(t.val[1], t.val[2]);                          \
        z = vsubq_s32(t.val[3], t.val[0]);                                  \
    } while (0)

#define mont_mul_neon_x3(z, x1, y1, x2, y2, x3, y3, t1, t2)                   \
    do {                                                                      \
        t2.val[0] = vmull_s32(vget_low_s32(x1), vget_low_s32(y1));            \
        t2.val[1] = vmull_high_s32(x1, y1);                                   \
        t2.val[0] = vaddq_s64(t2.val[0],                                      \
                              vmull_s32(vget_low_s32(x2), vget_low_s32(y2))); \
        t2.val[1] = vaddq_s64(t2.val[1], vmull_high_s32(x2, y2));             \
        t2.val[0] = vaddq_s64(t2.val[0],                                      \
                              vmull_s32(vget_low_s32(x3), vget_low_s32(y3))); \
        t2.val[1] = vaddq_s64(t2.val[1], vmull_high_s32(x3, y3));             \
        t1.val[2] = vuzp1q_s32((int32x4_t)t2.val[0], (int32x4_t)t2.val[1]);   \
        t1.val[3] = vuzp2q_s32((int32x4_t)t2.val[0], (int32x4_t)t2.val[1]);   \
        t1.val[0] = vmulq_n_s32(t1.val[2], GOODS_MONT_PM_PRIME);              \
        t1.val[1] = (int32x4_t)vmull_n_s32(vget_low_s32(t1.val[0]), GOODS_Q); \
        t1.val[2] = (int32x4_t)vmull_high_n_s32(t1.val[0], GOODS_Q);          \
        t1.val[0] = vuzp2q_s32(t1.val[1], t1.val[2]);                         \
        z = vsubq_s32(t1.val[3], t1.val[0]);                                  \
    } while (0)

int32_t mont_mul(int32_t x, int32_t y) {
    int64_t z = (int64_t)x * y;
    int64_t l = ((z & UINT32_MAX) * GOODS_MONT_MM_PRIME) & UINT32_MAX;
    return (z + l * GOODS_Q) >> GOODS_MONT_R_POW;
}

int32_t centered_goods_q(int32_t x) {
    int32_t y = GOODS_Q_HALF - x;
    uint32_t mask = -((uint32_t)y >> 31);
    return x - (mask & GOODS_Q);
}

#define layer_permute1(z, x2, x1, a, b, c, d)                    \
    do {                                                         \
        x2[i].val[0] = v##z##1q_s32(x1[i].val[a], x1[i].val[b]); \
        x2[i].val[1] = v##z##2q_s32(x1[i].val[a], x1[i].val[b]); \
        x2[i].val[2] = v##z##1q_s32(x1[i].val[c], x1[i].val[d]); \
        x2[i].val[3] = v##z##2q_s32(x1[i].val[c], x1[i].val[d]); \
    } while (0)

#define layer_permute2(z, x2, x1, a, b, c, d)                    \
    do {                                                         \
        x2[i].val[0] = v##z##1q_s32(x1[i].val[a], x1[i].val[b]); \
        x2[i].val[1] = v##z##1q_s32(x1[i].val[c], x1[i].val[d]); \
        x2[i].val[2] = v##z##2q_s32(x1[i].val[a], x1[i].val[b]); \
        x2[i].val[3] = v##z##2q_s32(x1[i].val[c], x1[i].val[d]); \
    } while (0)

#define cooley_tukey_x4(x, omega, t1, t2)                          \
    do {                                                           \
        mont_mul_neon_n(t1.val[0], x[i + size].val[0], omega, t2); \
        mont_mul_neon_n(t1.val[1], x[i + size].val[1], omega, t2); \
        mont_mul_neon_n(t1.val[2], x[i + size].val[2], omega, t2); \
        mont_mul_neon_n(t1.val[3], x[i + size].val[3], omega, t2); \
        x[i + size].val[0] = vsubq_s32(x[i].val[0], t1.val[0]);    \
        x[i + size].val[1] = vsubq_s32(x[i].val[1], t1.val[1]);    \
        x[i + size].val[2] = vsubq_s32(x[i].val[2], t1.val[2]);    \
        x[i + size].val[3] = vsubq_s32(x[i].val[3], t1.val[3]);    \
        x[i].val[0] = vaddq_s32(x[i].val[0], t1.val[0]);           \
        x[i].val[1] = vaddq_s32(x[i].val[1], t1.val[1]);           \
        x[i].val[2] = vaddq_s32(x[i].val[2], t1.val[2]);           \
        x[i].val[3] = vaddq_s32(x[i].val[3], t1.val[3]);           \
    } while (0)

#define cooley_tukey_x2(mul, x, omega1, omega2, t1, t2)  \
    do {                                                 \
        mul(t1.val[0], x[i].val[1], omega1, t2);         \
        mul(t1.val[1], x[i].val[3], omega2, t2);         \
        x[i].val[1] = vsubq_s32(x[i].val[0], t1.val[0]); \
        x[i].val[3] = vsubq_s32(x[i].val[2], t1.val[1]); \
        x[i].val[0] = vaddq_s32(x[i].val[0], t1.val[0]); \
        x[i].val[2] = vaddq_s32(x[i].val[2], t1.val[1]); \
    } while (0)

#define gentleman_sande_x4(x, omega, t1, t2)                                  \
    do {                                                                      \
        t1 = x[i];                                                            \
        x[i].val[0] = vaddq_s32(t1.val[0], x[i + size].val[0]);               \
        x[i].val[1] = vaddq_s32(t1.val[1], x[i + size].val[1]);               \
        x[i].val[2] = vaddq_s32(t1.val[2], x[i + size].val[2]);               \
        x[i].val[3] = vaddq_s32(t1.val[3], x[i + size].val[3]);               \
        mont_mul_neon_n(x[i + size].val[0],                                   \
                        vsubq_s32(t1.val[0], x[i + size].val[0]), omega, t2); \
        mont_mul_neon_n(x[i + size].val[1],                                   \
                        vsubq_s32(t1.val[1], x[i + size].val[1]), omega, t2); \
        mont_mul_neon_n(x[i + size].val[2],                                   \
                        vsubq_s32(t1.val[2], x[i + size].val[2]), omega, t2); \
        mont_mul_neon_n(x[i + size].val[3],                                   \
                        vsubq_s32(t1.val[3], x[i + size].val[3]), omega, t2); \
    } while (0)

#define gentleman_sande_no_mul_x4(x, t1)                               \
    do {                                                               \
        t1 = x[i];                                                     \
        x[i].val[0] = vaddq_s32(t1.val[0], x[i + size].val[0]);        \
        x[i].val[1] = vaddq_s32(t1.val[1], x[i + size].val[1]);        \
        x[i].val[2] = vaddq_s32(t1.val[2], x[i + size].val[2]);        \
        x[i].val[3] = vaddq_s32(t1.val[3], x[i + size].val[3]);        \
        x[i + size].val[0] = vsubq_s32(t1.val[0], x[i + size].val[0]); \
        x[i + size].val[1] = vsubq_s32(t1.val[1], x[i + size].val[1]); \
        x[i + size].val[2] = vsubq_s32(t1.val[2], x[i + size].val[2]); \
        x[i + size].val[3] = vsubq_s32(t1.val[3], x[i + size].val[3]); \
    } while (0)

#define gentleman_sande_x2(mul, x, omega1, omega2, t1, t2)               \
    do {                                                                 \
        t1.val[0] = x[i].val[0];                                         \
        t1.val[2] = x[i].val[2];                                         \
        x[i].val[0] = vaddq_s32(t1.val[0], x[i].val[1]);                 \
        x[i].val[2] = vaddq_s32(t1.val[2], x[i].val[3]);                 \
        mul(x[i].val[1], vsubq_s32(t1.val[0], x[i].val[1]), omega1, t2); \
        mul(x[i].val[3], vsubq_s32(t1.val[2], x[i].val[3]), omega2, t2); \
    } while (0)

void goods_ntt_512_x3(int32x4x4_t x[3][32]) {
    int32x4x4_t tx[3][32];
    int32x4x4_t t1, t2;

    // layer 2-5
    for (int size = 8; size > 0; size >>= 1) {
        int chunks = 16 / size;

        for (int offset = 0; offset < chunks; ++offset) {
            int32_t omega = goods_omegas[offset];

            int begin = 2 * size * offset;
            for (int i = begin; i < begin + size; ++i) {
                cooley_tukey_x4(x[0], omega, t1, t2);
                cooley_tukey_x4(x[1], omega, t1, t2);
                cooley_tukey_x4(x[2], omega, t1, t2);
            }
        }
    }

    // layer 6
    for (int i = 0; i < 32; ++i) {
        layer_permute1(zip, tx[0], x[0], 0, 1, 2, 3);
        layer_permute1(zip, tx[1], x[1], 0, 1, 2, 3);
        layer_permute1(zip, tx[2], x[2], 0, 1, 2, 3);

        int32_t omega = goods_omegas[i];
        cooley_tukey_x2(mont_mul_neon_n, tx[0], omega, omega, t1, t2);
        cooley_tukey_x2(mont_mul_neon_n, tx[1], omega, omega, t1, t2);
        cooley_tukey_x2(mont_mul_neon_n, tx[2], omega, omega, t1, t2);
    }

    // layer 7
    for (int i = 0; i < 32; ++i) {
        layer_permute1(zip, x[0], tx[0], 0, 2, 1, 3);
        layer_permute1(zip, x[1], tx[1], 0, 2, 1, 3);
        layer_permute1(zip, x[2], tx[2], 0, 2, 1, 3);

        int32_t omega1 = goods_omegas[2 * i];
        int32_t omega2 = goods_omegas[2 * i + 1];
        cooley_tukey_x2(mont_mul_neon_n, x[0], omega1, omega2, t1, t2);
        cooley_tukey_x2(mont_mul_neon_n, x[1], omega1, omega2, t1, t2);
        cooley_tukey_x2(mont_mul_neon_n, x[2], omega1, omega2, t1, t2);
    }

    // layer 8
    for (int i = 0; i < 32; ++i) {
        layer_permute1(uzp, tx[0], x[0], 0, 1, 2, 3);
        layer_permute1(uzp, tx[1], x[1], 0, 1, 2, 3);
        layer_permute1(uzp, tx[2], x[2], 0, 1, 2, 3);

        int32x2x2_t omega12 = vld2_dup_s32(&goods_omegas[4 * i]);
        int32x2x2_t omega34 = vld2_dup_s32(&goods_omegas[4 * i + 2]);
        int32x4_t omega1122 = vcombine_s32(omega12.val[0], omega12.val[1]);
        int32x4_t omega3344 = vcombine_s32(omega34.val[0], omega34.val[1]);
        cooley_tukey_x2(mont_mul_neon, tx[0], omega1122, omega3344, t1, t2);
        cooley_tukey_x2(mont_mul_neon, tx[1], omega1122, omega3344, t1, t2);
        cooley_tukey_x2(mont_mul_neon, tx[2], omega1122, omega3344, t1, t2);
    }

    // layer 9
    for (int i = 0; i < 32; ++i) {
        layer_permute1(uzp, x[0], tx[0], 0, 2, 1, 3);
        layer_permute1(uzp, x[1], tx[1], 0, 2, 1, 3);
        layer_permute1(uzp, x[2], tx[2], 0, 2, 1, 3);

        int32x4_t omega1234 = vld1q_s32(&goods_omegas[8 * i]);
        int32x4_t omega5678 = vld1q_s32(&goods_omegas[8 * i + 4]);
        int32x4_t omega1357 = vuzp1q_s32(omega1234, omega5678);
        int32x4_t omega2368 = vuzp2q_s32(omega1234, omega5678);
        cooley_tukey_x2(mont_mul_neon, x[0], omega1357, omega2368, t1, t2);
        cooley_tukey_x2(mont_mul_neon, x[1], omega1357, omega2368, t1, t2);
        cooley_tukey_x2(mont_mul_neon, x[2], omega1357, omega2368, t1, t2);
    }
}

void goods_schoolbook(int32x4x4_t x[3][32], int32x4x4_t y[3][32]) {
    int32x4x4_t t1, t2;
    int64x2x2_t t3;
    for (int i = 0; i < 32; ++i) {
        for (int j = 0; j < 4; ++j) {
            mont_mul_neon_x3(t1.val[0],                       //
                             x[0][i].val[j], y[0][i].val[j],  //
                             x[1][i].val[j], y[2][i].val[j],  //
                             x[2][i].val[j], y[1][i].val[j],  //
                             t2, t3);
            mont_mul_neon_x3(t1.val[1],                       //
                             x[0][i].val[j], y[1][i].val[j],  //
                             x[1][i].val[j], y[0][i].val[j],  //
                             x[2][i].val[j], y[2][i].val[j],  //
                             t2, t3);
            mont_mul_neon_x3(t1.val[2],                       //
                             x[0][i].val[j], y[2][i].val[j],  //
                             x[2][i].val[j], y[0][i].val[j],  //
                             x[1][i].val[j], y[1][i].val[j],  //
                             t2, t3);
            x[0][i].val[j] = t1.val[0];
            x[1][i].val[j] = t1.val[1];
            x[2][i].val[j] = t1.val[2];
        }
    }
}

void goods_intt_512_x3(int32x4x4_t x[3][32]) {
    int32x4x4_t tx[3][32];
    int32x4x4_t t1, t2;

    // layer 1
    for (int i = 0; i < 32; ++i) {
        int32x4_t omega1234 = vld1q_s32(&goods_inv_omegas[8 * i]);
        int32x4_t omega5678 = vld1q_s32(&goods_inv_omegas[8 * i + 4]);
        int32x4_t omega1357 = vuzp1q_s32(omega1234, omega5678);
        int32x4_t omega2368 = vuzp2q_s32(omega1234, omega5678);
        gentleman_sande_x2(mont_mul_neon, x[0], omega1357, omega2368, t1, t2);
        gentleman_sande_x2(mont_mul_neon, x[1], omega1357, omega2368, t1, t2);
        gentleman_sande_x2(mont_mul_neon, x[2], omega1357, omega2368, t1, t2);
    }

    // layer 2
    for (int i = 0; i < 32; ++i) {
        layer_permute2(zip, tx[0], x[0], 0, 1, 2, 3);
        layer_permute2(zip, tx[1], x[1], 0, 1, 2, 3);
        layer_permute2(zip, tx[2], x[2], 0, 1, 2, 3);

        int32x2x2_t omega12 = vld2_dup_s32(&goods_inv_omegas[4 * i]);
        int32x2x2_t omega34 = vld2_dup_s32(&goods_inv_omegas[4 * i + 2]);
        int32x4_t omega1122 = vcombine_s32(omega12.val[0], omega12.val[1]);
        int32x4_t omega3344 = vcombine_s32(omega34.val[0], omega34.val[1]);
        gentleman_sande_x2(mont_mul_neon, tx[0], omega1122, omega3344, t1, t2);
        gentleman_sande_x2(mont_mul_neon, tx[1], omega1122, omega3344, t1, t2);
        gentleman_sande_x2(mont_mul_neon, tx[2], omega1122, omega3344, t1, t2);
    }

    // layer 3
    for (int i = 0; i < 32; ++i) {
        layer_permute1(zip, x[0], tx[0], 0, 1, 2, 3);
        layer_permute1(zip, x[1], tx[1], 0, 1, 2, 3);
        layer_permute1(zip, x[2], tx[2], 0, 1, 2, 3);

        int32_t omega1 = goods_inv_omegas[2 * i];
        int32_t omega2 = goods_inv_omegas[2 * i + 1];
        gentleman_sande_x2(mont_mul_neon_n, x[0], omega1, omega2, t1, t2);
        gentleman_sande_x2(mont_mul_neon_n, x[1], omega1, omega2, t1, t2);
        gentleman_sande_x2(mont_mul_neon_n, x[2], omega1, omega2, t1, t2);
    }

    // layer 4
    for (int i = 0; i < 32; ++i) {
        layer_permute2(uzp, tx[0], x[0], 0, 1, 2, 3);
        layer_permute2(uzp, tx[1], x[1], 0, 1, 2, 3);
        layer_permute2(uzp, tx[2], x[2], 0, 1, 2, 3);

        int32_t omega = goods_inv_omegas[i];
        gentleman_sande_x2(mont_mul_neon_n, tx[0], omega, omega, t1, t2);
        gentleman_sande_x2(mont_mul_neon_n, tx[1], omega, omega, t1, t2);
        gentleman_sande_x2(mont_mul_neon_n, tx[2], omega, omega, t1, t2);

        layer_permute1(uzp, x[0], tx[0], 0, 1, 2, 3);
        layer_permute1(uzp, x[1], tx[1], 0, 1, 2, 3);
        layer_permute1(uzp, x[2], tx[2], 0, 1, 2, 3);
    }

    // layer 5-8
    for (int size = 1; size < 16; size <<= 1) {
        int chunks = 16 / size;

        for (int offset = 0; offset < chunks; ++offset) {
            int32_t omega = goods_inv_omegas[offset];

            int begin = 2 * size * offset;
            for (int i = begin; i < begin + size; ++i) {
                gentleman_sande_x4(x[0], omega, t1, t2);
                gentleman_sande_x4(x[1], omega, t1, t2);
                gentleman_sande_x4(x[2], omega, t1, t2);
            }
        }
    }

    // layer 9
    for (int i = 0; i < 16; ++i) {
        int size = 16;
        gentleman_sande_no_mul_x4(x[0], t1);
        gentleman_sande_no_mul_x4(x[1], t1);
        gentleman_sande_no_mul_x4(x[2], t1);
    }

    // divide 512
    for (int i = 0; i < 32; ++i) {
        for (int j = 0; j < 4; ++j) {
            mont_mul_neon_n(x[0][i].val[j], x[0][i].val[j],
                            GOODS_MONT_R_SQUARE_INV_N_MOD_Q, t1);
            mont_mul_neon_n(x[1][i].val[j], x[1][i].val[j],
                            GOODS_MONT_R_SQUARE_INV_N_MOD_Q, t1);
            mont_mul_neon_n(x[2][i].val[j], x[2][i].val[j],
                            GOODS_MONT_R_SQUARE_INV_N_MOD_Q, t1);
        }
    }
}

// TODO: optimize barrett_q with neon
void mul_small_goods(Fq *h, const Fq *f, const small *g) {
    int32_t pf[3][GOODS_N] = {0}, pg[3][GOODS_N] = {0};

    // layer_permute1 to layer 1
    goods_permute_to_layer1(pf, f);
    goods_permute_to_layer1(pg, g);

    // load
    int32x4x4_t vpf[3][32], vpg[3][32];
    for (int i = 0; i < 32; ++i) {
        vpf[0][i] = vld4q_s32(&pf[0][i * 16]);
        vpf[1][i] = vld4q_s32(&pf[1][i * 16]);
        vpf[2][i] = vld4q_s32(&pf[2][i * 16]);
        vpg[0][i] = vld4q_s32(&pg[0][i * 16]);
        vpg[1][i] = vld4q_s32(&pg[1][i * 16]);
        vpg[2][i] = vld4q_s32(&pg[2][i * 16]);
    }

    // ntt
    goods_ntt_512_x3(vpf);
    goods_ntt_512_x3(vpg);

    // goods_schoolbook
    goods_schoolbook(vpf, vpg);

    // intt
    goods_intt_512_x3(vpf);

    // store
    for (int i = 0; i < 32; ++i) {
        vst4q_s32(&pf[0][i * 16], vpf[0][i]);
        vst4q_s32(&pf[1][i * 16], vpf[1][i]);
        vst4q_s32(&pf[2][i * 16], vpf[2][i]);
    }

    // unpermutation
    int32_t hh[GOODS_P];
    goods_unpermute(hh, pf);

    // mod poly
    for (int i = GOODS_P - 1; i >= NTRU_LPRIME_P; --i) {
        hh[i - NTRU_LPRIME_P] =
            mont_mul(hh[i - NTRU_LPRIME_P] + hh[i], GOODS_MONT_R_MOD_Q);
        hh[i - NTRU_LPRIME_P + 1] =
            mont_mul(hh[i - NTRU_LPRIME_P + 1] + hh[i], GOODS_MONT_R_MOD_Q);
    }

    // mod q
    for (int i = 0; i < NTRU_LPRIME_P; ++i) {
        h[i] = barrett_q(centered_goods_q(hh[i]));
    }
}

#endif
