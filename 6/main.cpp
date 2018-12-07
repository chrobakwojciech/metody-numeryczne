#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

using namespace std;

int main() {
    const int size = 5;

    gsl_vector *diag, *e;
    diag = gsl_vector_alloc(size);
    // wektor na diagonali
    gsl_vector_set_all(diag, 2.0);

    // odejmujemy 1 od pierwszego i ostatniego (Sherman - Morrison)
    gsl_vector_set(diag, 0, 1.0);
    gsl_vector_set(diag, 4, 1.0);
    gsl_vector_set(diag, 2, 1.0);


    // odejmujemy wartosc wlasna na diagonali
    gsl_vector_add_constant(diag, -0.38197);

    // wektor nad i pod diagonala
    e = gsl_vector_alloc(size - 1);
    gsl_vector_set_all(e, 1.0);
    gsl_vector_set(e, 0, -1.0);
    gsl_vector_set(e, 3, -1.0);


    // wektory potrzebne do obliczenia Az = y ze wzoru Shermana - Morrisona
    gsl_vector *u, *v, *x, *q;
    u = gsl_vector_alloc(size);
    gsl_vector_set_all(u, 0);
    gsl_vector_set(u, 0, 1);
    gsl_vector_set(u, 4, 1);

    v = gsl_vector_alloc(size);
    gsl_vector_set_all(v, 0);
    gsl_vector_set(v, 0, 1);
    gsl_vector_set(v, 4, 1);

    x = gsl_vector_alloc(size);
    q = gsl_vector_alloc(size);


    // wektory potrzebne do wyliczenia wektora wlasnego
    gsl_vector *y, *temp;
    y = gsl_vector_alloc(size);
    temp = gsl_vector_alloc(size);

    gsl_vector_set_zero(y);
    // aby || y || = 1
    gsl_vector_set(y, 0, 1);


    // szukanie wektora wlasnego dla 0.38197
    while(true) {
        double x_prev_norm = gsl_blas_dnrm2(x);
        // Sherman - Morrison
        // Ax = y
        gsl_linalg_solve_symm_tridiag(diag, e, y, x);
        // Aq = u
        gsl_linalg_solve_symm_tridiag(diag, e, u, q);

        double data1 = 0;
        double data2 = 0;
        double *vx = &data1;
        double *vq = &data2;
        //vTx
        gsl_blas_ddot(v, x, vx);
        // vTq
        gsl_blas_ddot(v, q, vq);
        // 1 + vTq
        *vq += 1;

        gsl_vector_scale(q, *vx / *vq);
        gsl_vector_sub(x, q);
        // obliczony x

        // liczymy kolejny wektor y
        double norm = gsl_blas_dnrm2(x);
        gsl_vector_memcpy(temp, x);
        gsl_vector_scale(temp, 1 / norm);
        gsl_vector_memcpy(y, temp);


        if (abs(norm - x_prev_norm) < 10e-16)
            break;
    }

    cout << "Wektor wlasny dla wartosci 0.38197: " << endl;
    gsl_vector_fprintf(stdout, y, "%f");

    return 0;
}