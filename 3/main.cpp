#include <iostream>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

using namespace std;

void powerMethodMin(gsl_matrix *M) {
    int size = M->size1;

    gsl_vector *y, *z, *temp;
    gsl_permutation *p;
    int s;

    y = gsl_vector_alloc(size);
    temp = gsl_vector_alloc(size);

    gsl_vector_set_zero(y);
    // aby || y || = 1
    gsl_vector_set(y, 0, 1);


    z = gsl_vector_alloc(size);
    p = gsl_permutation_alloc(size);

    // LU
    gsl_linalg_LU_decomp(M, p, &s);

    while (true) {
        // poprzednia norma aby porownac i zakonczyc w odpowiednim momencie
        double z_prev_norm = gsl_blas_dnrm2(z);

        // rozwiazywanie Az = y
        gsl_linalg_LU_solve(M, p, y, z);

        // wyliczanie kolejnego y
        double norm = gsl_blas_dnrm2(z);
        gsl_vector_memcpy(temp, z);
        gsl_vector_scale(temp, 1 / norm);
        gsl_vector_memcpy(y, temp);


        if (abs(norm - z_prev_norm) < 10e-16)
            break;
    }
    cout << 1 / gsl_blas_dnrm2(z) << " dla:" << endl;
    gsl_vector_fprintf(stdout, y, "%f");
    cout << endl;


    // wyznaczanie wektora prostopadlego do wektora e1
    gsl_vector *e1;
    e1 = gsl_vector_alloc(size);
    gsl_vector_memcpy(e1, y);
    gsl_vector_set_all(y, 1);
    double sum = gsl_blas_dasum(e1) - gsl_vector_get(e1, size - 1);
    double buf = -(sum / gsl_vector_get(e1, size - 1));
    gsl_vector_set(y, size - 1, buf);
    double normY = gsl_blas_dnrm2(y);
    gsl_blas_dscal(1 / normY, y);


    while (true) {
        // poprzednia norma aby porownac i zakonczyc w odpowiednim momencie
        double z_prev_norm = gsl_blas_dnrm2(z);

        // rozwiazywanie Az = y
        gsl_linalg_LU_solve(M, p, y, z);


        // reortogonalizacja
        double data = 0;
        double *e1z = &data;
        gsl_blas_ddot(e1, z, e1z);
        gsl_vector *e1_e1z;
        e1_e1z = gsl_vector_alloc(size);
        gsl_vector_memcpy(e1_e1z, e1);
        gsl_blas_dscal(*e1z, e1_e1z);
        gsl_vector_sub(z, e1_e1z);

        // wyliczanie kolejnego y
        double norm = gsl_blas_dnrm2(z);
        gsl_vector_memcpy(temp, z);
        gsl_vector_scale(temp, 1 / norm);
        gsl_vector_memcpy(y, temp);


        if (abs(norm - z_prev_norm) < 10e-16)
            break;
    }
    cout << 1 / gsl_blas_dnrm2(z) << " dla:" << endl;
    gsl_vector_fprintf(stdout, y, "%f");
}


void powerMethodMax(gsl_matrix *M) {
    int size = M->size1;

    gsl_vector *y, *z, *temp;

    y = gsl_vector_alloc(size);
    temp = gsl_vector_alloc(size);
    gsl_vector_set_zero(y);
    gsl_vector_set(y, 0, 1);


    z = gsl_vector_alloc(size);


    while (true) {
        // poprzednia norma aby porownac i zakonczyc w odpowiednim momencie
        double z_prev_norm = gsl_blas_dnrm2(z);

        // wyliczanie Ay = z
        gsl_blas_dgemv(CblasNoTrans, 1.0, M, y, 0.0, z);

        // wyliczanie kolejnego y
        double norm = gsl_blas_dnrm2(z);
        gsl_vector_memcpy(temp, z);
        gsl_vector_scale(temp, 1 / norm);
        gsl_vector_memcpy(y, temp);


        if (abs(norm - z_prev_norm) < 10e-16)
            break;
    }
    cout << gsl_blas_dnrm2(z) << " dla:" << endl;
    gsl_vector_fprintf(stdout, y, "%f");
    cout << endl;

    // wyznaczanie wektora prostopadlego do wektora e1
    gsl_vector *e1;
    e1 = gsl_vector_alloc(size);
    gsl_vector_memcpy(e1, y);
    gsl_vector_set_all(y, 1);
    double sum = gsl_blas_dasum(e1) - gsl_vector_get(e1, size - 1);
    double buf = -(sum / gsl_vector_get(e1, size - 1));
    gsl_vector_set(y, size - 1, buf);
    double normY = gsl_blas_dnrm2(y);
    gsl_blas_dscal(1 / normY, y);


    while (true) {
        // poprzednia norma aby porownac i zakonczyc w odpowiednim momencie
        double z_prev_norm = gsl_blas_dnrm2(z);

        // wyliczanie Ay = z
        gsl_blas_dgemv(CblasNoTrans, 1.0, M, y, 0.0, z);

        // reortogonalizacja
        double data = 0;
        double *e1z = &data;
        gsl_blas_ddot(e1, z, e1z);
        gsl_vector *e1_e1z;
        e1_e1z = gsl_vector_alloc(size);
        gsl_vector_memcpy(e1_e1z, e1);
        gsl_blas_dscal(*e1z, e1_e1z);
        gsl_vector_sub(z, e1_e1z);

        // wyliczanie kolejnego y
        double norm = gsl_blas_dnrm2(z);
        gsl_vector_memcpy(temp, z);
        gsl_vector_scale(temp, 1 / norm);
        gsl_vector_memcpy(y, temp);


        if (abs(norm - z_prev_norm) < 10e-16)
            break;
    }
    cout << gsl_blas_dnrm2(z) << " dla:" << endl;
    gsl_vector_fprintf(stdout, y, "%f");
}


int main() {
    vector<vector<double>> A_ = {{19.0 / 12.0,  13.0 / 12.0,  5.0 / 6.0,  5.0 / 6.0,  13.0 / 12.0,  -17.0 / 12.0},
                                 {13.0 / 12.0,  13.0 / 12.0,  5.0 / 6.0,  5.0 / 6.0,  -11.0 / 12.0, 13.0 / 12.0},
                                 {5.0 / 6.0,    5.0 / 6.0,    5.0 / 6.0,  -1.0 / 6.0, 5.0 / 6.0,    5.0 / 6.0},
                                 {5.0 / 6.0,    5.0 / 6.0,    -1.0 / 6.0, 5.0 / 6.0,  5.0 / 6.0,    5.0 / 6.0},
                                 {13.0 / 12.0,  -11.0 / 12.0, 5.0 / 6.0,  5.0 / 6.0,  13.0 / 12.0,  13.0 / 12.0},
                                 {-17.0 / 12.0, 13.0 / 12.0,  5.0 / 6.0,  5.0 / 6.0,  13.0 / 12.0,  19.0 / 12.0}};
    int A_size = A_.size();
    vector<vector<double>> B_ = {{1, 2, 1, 1, 3, 4},
                                 {2, 2, 3, 5, 5, 6},
                                 {1, 3, 4, 4, 5, 5},
                                 {1, 5, 4, 5, 4, 6},
                                 {3, 5, 5, 4, 6, 2},
                                 {4, 6, 5, 6, 2, 0}};

    int B_size = B_.size();

    gsl_matrix *A;
    A = gsl_matrix_alloc(A_size, A_size);

    for (int i = 0; i < A_size; ++i) {
        for (int j = 0; j < A_size; ++j) {
            gsl_matrix_set(A, i, j, A_[i][j]);
        }
    }

    gsl_matrix *B;
    B = gsl_matrix_alloc(B_size, B_size);

    for (int i = 0; i < A_size; ++i) {
        for (int j = 0; j < A_size; ++j) {
            gsl_matrix_set(B, i, j, B_[i][j]);
        }
    }
    cout << "Najwieksze wartosci wlasne A: " << endl;
    powerMethodMax(A);
    cout << "\nNajmniejsze wartosci wlasne A: " << endl;
    powerMethodMin(A);
    cout << "\n\nNajwieksze wartosci wlasne B: " << endl;
    powerMethodMax(B);
    cout << "\nNajmniejsze wartosci wlasne B: " << endl;
    powerMethodMin(B);


    return 0;
}