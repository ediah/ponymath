#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef AVX
#include <xmmintrin.h>
#include <emmintrin.h>
#include <immintrin.h>
#endif

void debug_int(int x) {
    printf("debug\t%d\n", x);
}

void debug_point(double * p) {
    printf("debug\t%p\n", p);
}

extern double * mtxmul (double* a, double* b, int m, int n, int k)
{
    printf("multiplying\n");
    double * c = calloc(m * k, sizeof(double));

    #ifdef AVX
    double * prod = _mm_malloc(4*sizeof(double), 256);

    __m256d *p1, *p2, v1, v2;
    int M = n - (n % 4);
    #endif

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            c[i * k + j] = 0;
            #ifdef AVX
            p1 = (__m256d*) &(a[i * n]);
            p2 = (__m256d*) &(b[j * n]);
            v1 = _mm256_setzero_pd();

            for (int run = 0; run < M; run += 4) {
                v2 = _mm256_mul_pd(*p1, *p2);
                v1 = _mm256_add_pd(v1, v2);
                p1++;
                p2++;
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
    printf("done\n");
    return c;
}

double * realign(double * src, size_t n) {
    double * dst = aligned_alloc(256, n * sizeof(double));
    memcpy(dst, src, n * sizeof(double));
    return dst;
}
