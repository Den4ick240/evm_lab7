#pragma once


class IMatOperations {
public:
    virtual void getI(float *I) = 0;

    virtual void trans(float *mat, float *res) = 0;

    virtual void div(float *mat, float v) = 0;

    virtual float maxI(float *A) = 0;

    virtual float maxJ(float *A) = 0;

    virtual void mul(float *A, float *B, float *C) = 0;

    virtual void add(float *a, float *b) = 0;

    virtual void sub(float *A, float *B, float *R) = 0;

    virtual void inverseMat(float *A, float *BB, float *R, float *R1, float *R2, float *I) = 0;
};