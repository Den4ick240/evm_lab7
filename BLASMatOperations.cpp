#include "BLASMatOperations.h"
#include "MatGenerator.h"

#ifdef BLAS_INSTALLED

namespace {
    float II[N2];
}

void BLASMatOperations::inverseMat(float *A, float *BB, float *R, float *R1, float *R2, float *I) {
    float *K[2] = {R1, R2};

    this->getI(I);


//    this->trans(A, BB);
//    this->div(BB, this->maxI(A) * this->maxJ(A));
    cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                N, N, N,
                1 / (this->maxI(A) * this->maxJ(A)), A, N,
                II, N,
                0.0, BB, N);

//    this->mul(BB, A, R);

    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                N, N, N,
                1.0, BB, N,
                A, N,
                0.0, R, N);

//    this->sub(I, R, R);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                N, N, N,
                1.0, I, N,
                I, N,
                -1.0, R, N);


//    this->add(I, R);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                N, N, N,
                1.0, R, N,
                II, N,
                1.0, I, N);

    bool w = false;

//    this->mul(R, R, K[w]);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                N, N, N,
                1.0, R, N,
                R, N,
                0, K[w], N);

//    this->add(I, K[w]);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                N, N, N,
                1.0, K[w], N,
                II, N,
                1.0, I, N);

    for (int i = 2; i < M; i++) {
//        this->mul(K[w], R, K[!w]);
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    N, N, N,
                    1.0, K[w], N,
                    R, N,
                    0, K[!w], N);

        w = !w;
    //    this->add(I, K[w]);

        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    N, N, N,
                    1.0, K[w], N,
                    II, N,
                    1, I, N);

    }



//    this->mul(I, BB, R);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                N, N, N,
                1.0, I, N,
                BB, N,
                0, R, N);
}

BLASMatOperations::BLASMatOperations() {
    this->getI(II);
}

#endif