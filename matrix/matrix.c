#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if defined(AVX) || defined(AVX2)
#include <xmmintrin.h>
#include <emmintrin.h>
#include <immintrin.h>
#endif

#include "matrix.h"

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

double * mtxmul (double* a, double* b, int m, int n, int k)
{
    double * c = calloc(m * k, sizeof(double));

    int block_size = LL_CACHE_SIZE / 2 / sizeof(double) / n;
    int row = 0, col = 0;

    for (row = 0; row < m; row += block_size) {
        for (col = 0; col < k; col += block_size) {
            mtxmul_l2(&a[row * n], &b[col * n], &c[row * k + col], m, n, k,
                      MIN(m - row, block_size), MIN(k - col, block_size));
        }
    }
    
    return c;
}

void mtxmul_l2
(double* a, double* b, double* c, int m, int n, int k, int rblock, int cblock)
{
    int block_size = L2_CACHE_SIZE / 2 / sizeof(double) / n;
    int row = 0, col = 0;

    for (row = 0; row < rblock; row += block_size) {
        for (col = 0; col < cblock; col += block_size) {
            mtxmul_l1(&a[row * n], &b[col * n], &c[row * k + col], m, n, k,
                         MIN(rblock - row, block_size),
                         MIN(cblock - col, block_size));
        }
    }
}

void mtxmul_l1
(double* a, double* b, double* c, int m, int n, int k, int rblock, int cblock)
{
    int block_size = L1_CACHE_SIZE / 2 / sizeof(double) / n;
    int row = 0, col = 0;

    for (row = 0; row < rblock; row += block_size) {
        for (col = 0; col < cblock; col += block_size) {
            mtxmul_micro(&a[row * n], &b[col * n], &c[row * k + col], m, n, k,
                         MIN(rblock - row, block_size),
                         MIN(cblock - col, block_size));
        }
    }
}

inline void mtxmul_micro (
    double* a,  // Left matrix A
    double* b,  // Right matrix B (transposed)
    double* c,  // Result C
    int m,      // Rows in A
    int n,      // Columns in A = Columns in B
    int k,      // Rows in B = Columns in C
    int rblock, // Size of micro core, Rows
    int cblock  // Size of micro core, Columns
) {
    #if defined(AVX) || defined(AVX2)
    double prod[4];
    const double *p1, *p2;
    __m256d v1, v2;
    int M = n - (n % 4);
    #endif

    for (int i = 0; i < MIN(m, rblock); i++) {
        for (int j = 0; j < MIN(k, cblock); j++) {
            c[i * k + j] = 0;
            #if defined(AVX) || defined(AVX2)
            p1 = &(a[i * n]);
            p2 = &(b[j * n]);
            v1 = _mm256_setzero_pd();

            for (int run = 0; run < M; run += 4) {
                #if !defined(AVX2)
                v2 = _mm256_mul_pd(
                    _mm256_loadu_pd(p1),
                    _mm256_loadu_pd(p2)
                );
                v1 = _mm256_add_pd(v1, v2);
                #else
                v1 = _mm256_fmadd_pd(
                    _mm256_loadu_pd(p1),
                    _mm256_loadu_pd(p2),
                    v1
                );
                #endif
                p1 += 4;
                p2 += 4;
            }
            _mm256_storeu_pd(prod, v1);

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
}
