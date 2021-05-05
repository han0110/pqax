#include "kem.h"

#include <stdint.h>

#include "encode.h"
#include "params.h"
#include "poly.h"
#include "rand.h"
#include "reduce.h"
#include "round.h"
#include "type.h"

#define SEED_BYTES 32
#define SECRET_KEY_BYTES (NTRU_LPRIME_P + 3) / 4
#define NTRU_LPRIME_Q_HALF (NTRU_LPRIME_Q - 1) / 2

static void expand(uint32_t *L, const uint8_t *k) {
    // TODO:
}

static void generator(Fq *G, const uint8_t *S) {
    uint32_t L[NTRU_LPRIME_P];
    expand(L, S);
    for (int i = 0; i < NTRU_LPRIME_P; i++)
        G[i] = barrett_reduce_q(L[i] - NTRU_LPRIME_Q_HALF);
}

static void key_gen_a(Fq *A, small *a, const Fq *G) {
    Fq aG[NTRU_LPRIME_P];

    random_small_xp(a);
    mul_small(aG, G, a);
    round3_xp(A, aG);
}

static void key_gen_s(uint8_t *S, Fq *A, small *a) {
    Fq G[NTRU_LPRIME_P];

    random_u8_xn(S, SEED_BYTES);
    generator(G, S);
    key_gen_a(A, a, G);
}

int crypto_kem_keypair(uint8_t *pk, uint8_t *sk) {
    Fq A[NTRU_LPRIME_P];
    small a[NTRU_LPRIME_P];

    key_gen_s(pk += SEED_BYTES, A, a);
    encode_round(pk, A);
    encode_small(sk += SECRET_KEY_BYTES, a);

    for (int i = 0; i < CRYPTO_PUBLICKEYBYTES; i++) *sk++ = pk[i];
    random_u8_xn(sk += CRYPTO_BYTES, CRYPTO_BYTES);
    // TODO: hash_prefix
}

int crypto_kem_enc(uint8_t *ct, uint8_t *ss, const uint8_t *pk) {
    // TODO:
}

int crypto_kem_dec(uint8_t *ss, const uint8_t *ct, const uint8_t *sk) {
    // TODO:
}
