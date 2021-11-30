#ifndef PORTABLE_RAND_INCLUDED__
#define PORTABLE_RAND_INCLUDED__

#include <stdint.h>

static unsigned long int portable_next = 1;

static inline int portable_rand(void) {
    portable_next = portable_next * 1103515245 + 12345;
    return (uint64_t)(portable_next/65536) % 32768;
}

static inline void portable_srand(uint64_t seed) {
    portable_next = seed;
}

static inline int portable_int2randint(uint64_t randseed, int n_iters) {
    portable_srand(randseed);
    for (int i = 0; i < n_iters; i++) {
        randseed = (uint64_t)portable_rand();
    }
    return randseed;
}

#endif

