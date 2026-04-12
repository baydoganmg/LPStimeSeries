/*
 * Per-thread random number generator for OpenMP parallelization.
 * Uses xoshiro256** (public domain by David Blackman and Sebastiano Vigna).
 * See: https://prng.di.unimi.it/
 */
#ifndef RNG_H
#define RNG_H

#include <stdint.h>

typedef struct {
    uint64_t s[4];
} rng_state_t;

static inline uint64_t rng_rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

/* xoshiro256** generator */
static inline uint64_t rng_next(rng_state_t *state) {
    const uint64_t result = rng_rotl(state->s[1] * 5, 7) * 9;
    const uint64_t t = state->s[1] << 17;

    state->s[2] ^= state->s[0];
    state->s[3] ^= state->s[1];
    state->s[1] ^= state->s[2];
    state->s[0] ^= state->s[3];

    state->s[2] ^= t;
    state->s[3] = rng_rotl(state->s[3], 45);

    return result;
}

/* Returns a uniform double in [0, 1) */
static inline double rng_uniform(rng_state_t *state) {
    return (rng_next(state) >> 11) * 0x1.0p-53;
}

/* Seed the RNG state from a single 64-bit seed using splitmix64 */
static inline void rng_seed(rng_state_t *state, uint64_t seed) {
    uint64_t z;
    int i;
    for (i = 0; i < 4; i++) {
        z = (seed += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        state->s[i] = z ^ (z >> 31);
    }
}

#endif /* RNG_H */
