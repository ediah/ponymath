#ifndef PONYMATH_MATRIX_H
#define PONYMATH_MATRIX_H

#define L1_CACHE_SIZE 128 * 1024
#define LL_CACHE_SIZE 4 * 1024 * 1024

#define MIN(x, y) ((x) < (y) ? (x) : (y))

double * realign
(double * src, size_t n);

double * mtxmul
(double* a, double* b, int m, int n, int k);

void mtxmul_l1
(double* a, double* b, double* c, int m, int n, int k, int rblock, int cblock);

void mtxmul_micro
(double* a, double* b, double* c, int m, int n, int k, int rblock, int cblock);

#endif
