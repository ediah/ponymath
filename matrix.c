#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef AVX
#include <xmmintrin.h>
#include <emmintrin.h>
#include <immintrin.h>
#endif

void debug_int(int x) {
    #ifdef DEBUG
    printf("debug\t%d\n", x);
    #endif
}

void debug_point(double * p) {
    #ifdef DEBUG
    printf("debug\t%p\n", p);
    #endif
}

double * realign(double * src, size_t n) {
    double * dst = _mm_malloc(n * sizeof(double), 32);
    memcpy(dst, src, n * sizeof(double));
    return dst;
}

extern double * mtxmul (double* a, double* b, int m, int n, int k)
{
    #ifdef DEBUG
    printf("multiplying\n");
    #endif
    double * c = calloc(m * k, sizeof(double));

    #ifdef AVX
    double * prod = _mm_malloc(4 * sizeof(double), 32);

    const double *p1, *p2;
    __m256d v1, v2;
    int M = n - (n % 4);
    #endif

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            c[i * k + j] = 0;
            #ifdef AVX
            p1 = &(a[i * n]);
            p2 = &(b[j * n]);
            v1 = _mm256_setzero_pd();

            for (int run = 0; run < M; run += 4) {
                v2 = _mm256_mul_pd(
                    _mm256_loadu_pd(p1),
                    _mm256_loadu_pd(p2)
                );
                v1 = _mm256_add_pd(v1, v2);
                p1 += 4;
                p2 += 4;
            }
            _mm256_store_pd(prod, v1);

            for (int run = 0; run < n % 4; run++) {
                prod[run] += a[i * n + M + run] * b[j * n + M + run];
            }

            for (int run = 0; run < 4; run++) {
                c[i * k + j] += prod[run];
                prod[run] = 0;
            }
            #else
            for (int run = 0; run < n; run++)
                c[i * k + j] += a[i * n + run] * b[j * n + run];
            #endif
        }
    }

    #ifdef DEBUG
    printf("done\n");
    #endif

    #ifdef AVX
    free(prod);
    #endif

    return c;
}
