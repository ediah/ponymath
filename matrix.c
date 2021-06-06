#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef AVX
#include <xmmintrin.h>
#include <emmintrin.h>
#include <immintrin.h>
#endif

extern double * mtxmul (double* a, double* b, int m, int n, int k)
{
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

    return c;
}

double * aligned_copy(double * src, size_t n) {
    double * dst = aligned_alloc(256, n * sizeof(double));
    memcpy(dst, src, n * sizeof(double));
    return dst;
}

extern void mtxmulpar
(double* a, double* b, double* c, int m, int n, int k,
 int l, int threads)
{
    int block = m * k / threads;
    int shift = block * (l - 1);
    double* from = c + shift;
    if (l == threads) block += (m * k) % threads;
    a = aligned_copy(a, m * n);
    b = aligned_copy(b, n * k);
    #ifdef AVX
    double * prod = aligned_alloc(64, 4 * sizeof(double));

    __m256d *p1, *p2, v1, v2;
    int M = n - (n % 4);
    #endif

    for (int i = 0; i < block; i++) {
        from[i] = 0;
        int ax = ((i + shift) / k) * n;
        int bx = ((i + shift) % k) * n;

        #ifdef AVX
        p1 = (__m256d*) &(a[ax]);
        p2 = (__m256d*) &(b[bx]);
        v1 = _mm256_setzero_pd();

        for (int run = 0; run < M; run += 4) {
            v2 = _mm256_mul_pd(*p1, *p2);
            v1 = _mm256_add_pd(v1, v2);
            p1++;
            p2++;
        }
        _mm256_store_pd(prod, v1);

        for (int run = 0; run < n % 4; run++) {
            prod[run] += a[ax + M + run] * b[bx + M + run];
        }

        for (int run = 0; run < 4; run++) {
            from[i] += prod[run];
            prod[run] = 0;
        }
        #else
        for (int run = 0; run < n; run++)
            from[i] += a[ax + run] * b[bx + run];
        #endif
    }

    free(a);
    free(b);
}
//clang -fprofile-instr-generate a.c && LLVM_PROFILE_FILE=a.profraw ./a.out
//llvm-profdata merge -output=a.profdata a.profraw
//clang -fprofile-instr-use=a.profdata a.c
