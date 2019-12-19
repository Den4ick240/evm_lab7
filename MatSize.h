#pragma once

//#define TESTING
#define BLAS_INSTALLED

#ifdef TESTING

#define WRITE_MAT
#define READ_MAT

const unsigned int N = 4;
const unsigned int  M = 10000000;

#else

const unsigned int N = 2048;
const unsigned int  M = 4;

#endif

#define N2 (N * N)

