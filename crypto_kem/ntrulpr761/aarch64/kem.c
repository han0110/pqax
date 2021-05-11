#include "kem.h"

#include <stdint.h>

#include "../../../mupq/common/aes.h"
#include "decode.h"
#include "encode.h"
#include "params.h"
#include "poly.h"
#include "rand.h"
#include "reduce.h"
#include "round.h"
#include "sha2.h"
#include "sort.h"
#include "type.h"

#define AES256_CTR_NONCE_BYTES 12
#define HASH_BYTES 32

static const unsigned char iv[AES256_CTR_NONCE_BYTES] = {0};

static void expand(uint32_t *L, const uint8_t *k) {
    aes256ctx ctx;
    aes256_ctr_keyexp(&ctx, k);
    aes256_ctr(L, 4 * NTRU_LPRIME_P, iv, &ctx);
    aes256_ctx_release(&ctx);

    for (int i = 0; i < NTRU_LPRIME_P; ++i) {
        uint32_t L1 = ((unsigned char *)L)[4 * i + 1];
        uint32_t L0 = ((unsigned char *)L)[4 * i];
        uint32_t L2 = ((unsigned char *)L)[4 * i + 2];
        uint32_t L3 = ((unsigned char *)L)[4 * i + 3];
        L[i] = L0 + (L1 << 8) + (L2 << 16) + (L3 << 24);
    }
}

static void generator(Fq *G, const uint8_t *S) {
    uint32_t L[NTRU_LPRIME_P];

    expand(L, S);
    for (int i = 0; i < NTRU_LPRIME_P; ++i)
        G[i] = barrett_q(L[i] - NTRU_LPRIME_Q_HALF_FLOOR);
}

static int8_t top(Fq C) {
    return (NTRU_LPRIME_TAU_1 * (int32_t)(C + NTRU_LPRIME_TAU_0) + 16384) >> 15;
}

static Fq right(int8_t T) {
    return barrett_q(NTRU_LPRIME_TAU_3 * (int32_t)T - NTRU_LPRIME_TAU_2);
}

static void hash_prefix(unsigned char *dst, int prefix,
                        const unsigned char *src, int srclen) {
    unsigned char x[srclen + 1], h[64];

    x[0] = prefix;
    for (int i = 0; i < srclen; ++i) x[i + 1] = src[i];
    sha512(h, x, srclen + 1);
    for (int i = 0; i < 32; ++i) dst[i] = h[i];
}

static void hash_short(small *dst, const int8_t *r) {
    unsigned char s[NTRU_LPRIME_INPUT_BYTES];
    unsigned char h[HASH_BYTES];
    uint32_t L[NTRU_LPRIME_P];

    encode_input(s, r);
    hash_prefix(h, 5, s, sizeof s);
    expand(L, h);
    sort_to_short(dst, L);
}

static void hash_confirm(unsigned char *h, const unsigned char *r,
                         const unsigned char *pk,
                         const unsigned char *pk_hash) {
    int i;
    unsigned char x[NTRU_LPRIME_INPUT_BYTES + HASH_BYTES];

    for (i = 0; i < NTRU_LPRIME_INPUT_BYTES; ++i) x[i] = r[i];
    for (i = 0; i < HASH_BYTES; ++i)
        x[NTRU_LPRIME_INPUT_BYTES + i] = pk_hash[i];

    hash_prefix(h, 2, x, sizeof x);
}

static void key_gen_a(Fq *A, small *a, const Fq *G) {
    Fq aG[NTRU_LPRIME_P];

    random_short(a);
    mul_small(aG, G, a);
    round3_xp(A, aG);
}

static void key_gen_s(uint8_t *S, Fq *A, small *a) {
    Fq G[NTRU_LPRIME_P];

    random_u8_xn(S, NTRU_LPRIME_SEED_BYTES);
    generator(G, S);
    key_gen_a(A, a, G);
}

static void encrypt(Fq *B, int8_t *T, const int8_t *r, const unsigned char *S,
                    const Fq *A) {
    small b[NTRU_LPRIME_P];
    Fq G[NTRU_LPRIME_P], bG[NTRU_LPRIME_P], bA[NTRU_LPRIME_P];

    generator(G, S);
    hash_short(b, r);
    mul_small(bG, G, b);
    round3_xp(B, bG);
    mul_small(bA, A, b);
    for (int i = 0; i < 8 * NTRU_LPRIME_INPUT_BYTES; ++i)
        T[i] = top(barrett_q(bA[i] + r[i] * NTRU_LPRIME_Q_HALF_FLOOR));
}

#define neg_mask(x) -(int)(x >> 15)

static void decrypt(int8_t *r, const unsigned char *ct, unsigned char *sk) {
    small a[NTRU_LPRIME_P];
    Fq B[NTRU_LPRIME_P], aB[NTRU_LPRIME_P];
    int8_t T[8 * NTRU_LPRIME_INPUT_BYTES];

    decode_small(a, sk);
    decode_round(B, ct);
    ct += NTRU_LPRIME_ROUND_ENCODED_BYTES;
    decode_top(T, ct);

    mul_small(aB, B, a);
    for (int i = 0; i < 8 * NTRU_LPRIME_INPUT_BYTES; ++i)
        r[i] = neg_mask(barrett_q(right(T[i]) - aB[i] + 4 * NTRU_LPRIME_W + 1));
}

static void hide(uint8_t *ct, unsigned char *r_encoded, const int8_t *r,
                 const uint8_t *pk, unsigned char *pk_hash) {
    Fq A[NTRU_LPRIME_P];
    Fq B[NTRU_LPRIME_P];
    int8_t T[8 * CRYPTO_CIPHERTEXTBYTES];

    encode_input(r_encoded, r);

    decode_round(A, pk);
    pk += NTRU_LPRIME_SEED_BYTES;
    encrypt(B, T, r, pk, A);
    encode_round(ct, B);
    ct += NTRU_LPRIME_ROUND_ENCODED_BYTES;
    encode_top(ct, T);
    ct += 4 * NTRU_LPRIME_INPUT_BYTES;

    hash_confirm(ct, r_encoded, pk, pk_hash);
}

static void hash_session(unsigned char *ss, int b,
                         const unsigned char *r_encoded,
                         const unsigned char *ct) {
    int i;
    unsigned char x[CRYPTO_CIPHERTEXTBYTES + CRYPTO_BYTES];

    for (i = 0; i < CRYPTO_BYTES; ++i) x[i] = r_encoded[i];
    for (i = CRYPTO_BYTES; i < CRYPTO_CIPHERTEXTBYTES + CRYPTO_BYTES; ++i)
        x[i] = ct[i];
    hash_prefix(ss, b, x, CRYPTO_CIPHERTEXTBYTES + CRYPTO_BYTES);
}

static int ciphertext_diff_mask(const unsigned char *ct,
                                const unsigned char *ct_ref) {
    int len = CRYPTO_CIPHERTEXTBYTES;
    uint16_t mask = 0;

    while (len-- > 0) mask |= (*ct++) ^ (*ct_ref++);
    return (1 & ((mask - 1) >> 8)) - 1;
}

int crypto_kem_keypair(uint8_t *pk, uint8_t *sk) {
    Fq A[NTRU_LPRIME_P];
    small a[NTRU_LPRIME_P];

    key_gen_s(pk, A, a);
    encode_round(pk + NTRU_LPRIME_SEED_BYTES, A);
    encode_small(sk, a);
    sk += NTRU_LPRIME_SECRET_KEY_BYTES;

    for (int i = 0; i < CRYPTO_PUBLICKEYBYTES; ++i) *sk++ = pk[i];
    random_u8_xn(sk, CRYPTO_BYTES);
    sk += CRYPTO_BYTES;
    hash_prefix(sk, 4, pk, CRYPTO_PUBLICKEYBYTES);

    return 0;
}

int crypto_kem_enc(uint8_t *ct, uint8_t *ss, const uint8_t *pk) {
    int8_t r[8 * NTRU_LPRIME_INPUT_BYTES];
    unsigned char r_encoded[NTRU_LPRIME_INPUT_BYTES];
    unsigned char pk_hash[HASH_BYTES];

    hash_prefix(pk_hash, 4, pk, CRYPTO_PUBLICKEYBYTES);
    random_input(r);
    hide(ct, r_encoded, r, pk, pk_hash);
    hash_session(ss, 1, r_encoded, ct);

    return 0;
}

int crypto_kem_dec(uint8_t *ss, const uint8_t *ct, const uint8_t *sk) {
    const unsigned char *pk = sk + NTRU_LPRIME_SECRET_KEY_BYTES;
    const unsigned char *rho = pk + CRYPTO_PUBLICKEYBYTES;
    const unsigned char *pk_hash = rho + NTRU_LPRIME_INPUT_BYTES;

    int8_t r[8 * NTRU_LPRIME_INPUT_BYTES];
    unsigned char r_encoded[NTRU_LPRIME_INPUT_BYTES];
    unsigned char ct_ref[CRYPTO_CIPHERTEXTBYTES];

    decrypt(r, ct, sk);
    hide(ct_ref, r_encoded, r, pk, pk_hash);
    int mask = ciphertext_diff_mask(ct, ct_ref);
    for (int i = 0; i < NTRU_LPRIME_INPUT_BYTES; ++i)
        r_encoded[i] ^= mask & (r_encoded[i] & rho[i]);
    hash_session(ss, 1 + mask, r_encoded, ct);

    return 0;
}
