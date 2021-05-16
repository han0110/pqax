#ifndef POLY_H
#define POLY_H

#include <stdint.h>

#include "params.h"
#include "reduce.h"
#include "type.h"

#define GOODS_N 512
#define GOODS_N_HALF 256
#define GOODS_P 1536
// 4547 * 1536 + 1
#define GOODS_Q 6984193
#define GOODS_MONT_R_POW 32
// pow(-6984193, -1, 1 << 32)
#define GOODS_MONT_M_PRIME 2368115199
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

void mul_small_schoolbook(Fq *h, const Fq *f, const small *g) {
    Fq f_reduced[NTRU_LPRIME_P];
    Fq fg[2 * NTRU_LPRIME_P - 1];
    int32_t result;
    int i, j;

    for (i = 0; i < NTRU_LPRIME_P; ++i) {
        f_reduced[i] = barrett_q(f[i]);
    }
    for (i = 0; i < NTRU_LPRIME_P; ++i) {
        result = 0;
        for (j = 0; j <= i; ++j)
            result = result + f_reduced[j] * (int32_t)g[i - j];
        fg[i] = barrett_q(result);
    }
    for (i = NTRU_LPRIME_P; i < 2 * NTRU_LPRIME_P - 1; ++i) {
        result = 0;
        for (j = i - NTRU_LPRIME_P + 1; j < NTRU_LPRIME_P; ++j)
            result = result + f_reduced[j] * (int32_t)g[i - j];
        fg[i] = barrett_q(result);
    }

    for (i = 2 * NTRU_LPRIME_P - 2; i >= NTRU_LPRIME_P; --i) {
        fg[i - NTRU_LPRIME_P] = barrett_q(fg[i - NTRU_LPRIME_P] + fg[i]);
        fg[i - NTRU_LPRIME_P + 1] =
            barrett_q(fg[i - NTRU_LPRIME_P + 1] + fg[i]);
    }

    for (i = 0; i < NTRU_LPRIME_P; ++i) h[i] = fg[i];
}

int32_t mont_mul(int32_t x, int32_t y) {
    int64_t z = (int64_t)x * y;
    int64_t l = ((int64_t)(z & UINT32_MAX) * GOODS_MONT_M_PRIME) & UINT32_MAX;
    return (z + l * GOODS_Q) >> GOODS_MONT_R_POW;
}

int32_t centered_goods_q(int32_t x) {
    int32_t y = GOODS_Q / 2 - x;
    uint32_t mask = -((uint32_t)y >> 31);
    return x - (mask & GOODS_Q);
}

int32_t mont_mul_sum_x3(int32_t x1, int32_t y1, int32_t x2, int32_t y2,
                        int32_t x3, int32_t y3) {
    int64_t z = (int64_t)x1 * y1;
    z += (int64_t)x2 * y2;
    z += (int64_t)x3 * y3;
    int64_t l = ((int64_t)(z & UINT32_MAX) * GOODS_MONT_M_PRIME) & UINT32_MAX;
    return (z + l * GOODS_Q) >> GOODS_MONT_R_POW;
}

// TODO: optimize bufferfly and mont_mul with neon
void ntt_512_x3(int32_t x[3][GOODS_N]) {
    for (int size = GOODS_N_HALF; size > 0; size >>= 1) {
        int chunks = GOODS_N_HALF / size;

        for (int offset = 0; offset < chunks; ++offset) {
            int begin = 2 * size * offset;

            for (int i = begin; i < begin + size; ++i) {
                int32_t omega = goods_omegas[offset], ox;

                ox = mont_mul(x[0][i + size], omega);
                x[0][i] += ox;
                x[0][i + size] = x[0][i] - (ox << 1);

                ox = mont_mul(x[1][i + size], omega);
                x[1][i] += ox;
                x[1][i + size] = x[1][i] - (ox << 1);

                ox = mont_mul(x[2][i + size], omega);
                x[2][i] += ox;
                x[2][i + size] = x[2][i] - (ox << 1);
            }
        }
    }
}

// TODO: optimize bufferfly and mont_mul with neon
void intt_512_x3(int32_t x[3][GOODS_N]) {
    for (int size = 1; size < GOODS_N; size <<= 1) {
        int chunks = GOODS_N_HALF / size;

        for (int offset = 0; offset < chunks; ++offset) {
            int begin = 2 * size * offset;

            for (int i = begin; i < begin + size; ++i) {
                int32_t inv_omega = goods_inv_omegas[offset];

                x[0][i] += x[0][i + size];
                x[0][i + size] =
                    mont_mul((x[0][i] - (x[0][i + size] << 1)), inv_omega);

                x[1][i] += x[1][i + size];
                x[1][i + size] =
                    mont_mul((x[1][i] - (x[1][i + size] << 1)), inv_omega);

                x[2][i] += x[2][i + size];
                x[2][i + size] =
                    mont_mul((x[2][i] - (x[2][i + size] << 1)), inv_omega);
            }
        }
    }
}

// TODO: unroll permutation with first round
// TODO: unroll unpermutation
// TODO: optimize mont_mul with neon
void mul_small_goods(Fq *h, const Fq *f, const small *g) {
    // permutation
    int32_t pf[3][GOODS_N], pg[3][GOODS_N];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < GOODS_N; ++j) {
            int k = (1024 * i + 513 * j) % GOODS_P;
            pf[i][j] = k < NTRU_LPRIME_P ? f[k] : 0;
            pg[i][j] = k < NTRU_LPRIME_P ? g[k] : 0;
        }
    }

    // ntt
    ntt_512_x3(pf);
    ntt_512_x3(pg);

    // schoolbook
    for (int i = 0; i < GOODS_N; ++i) {
        int32_t s0 = mont_mul_sum_x3(pf[0][i], pg[0][i], pf[1][i], pg[2][i],
                                     pf[2][i], pg[1][i]);
        int32_t s1 = mont_mul_sum_x3(pf[0][i], pg[1][i], pf[1][i], pg[0][i],
                                     pf[2][i], pg[2][i]);
        int32_t s2 = mont_mul_sum_x3(pf[0][i], pg[2][i], pf[2][i], pg[0][i],
                                     pf[1][i], pg[1][i]);
        pf[0][i] = s0;
        pf[1][i] = s1;
        pf[2][i] = s2;
    }

    // intt
    intt_512_x3(pf);

    // unpermutation
    int32_t hh[GOODS_P] = {0};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < GOODS_N; ++j) {
            hh[(1024 * i + 513 * j) % GOODS_P] = pf[i][j];
        }
    }

    // mod poly
    for (int i = GOODS_P - 1; i >= NTRU_LPRIME_P; --i) {
        hh[i - NTRU_LPRIME_P] =
            mont_mul(hh[i - NTRU_LPRIME_P] + hh[i], GOODS_MONT_R_MOD_Q);
        hh[i - NTRU_LPRIME_P + 1] =
            mont_mul(hh[i - NTRU_LPRIME_P + 1] + hh[i], GOODS_MONT_R_MOD_Q);
    }

    // mod q
    for (int i = 0; i < NTRU_LPRIME_P; ++i) {
        h[i] = barrett_q(
            centered_goods_q(mont_mul(hh[i], GOODS_MONT_R_SQUARE_INV_N_MOD_Q)));
    }
}

void mul_small(Fq *h, const Fq *f, const small *g) { mul_small_goods(h, f, g); }

#endif
