#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

using namespace std;

const int size = 5;

int main() {
    double x[] = {-1.2300, -1.1900, -0.7400, 0.1100, 2.5600};
    double y[] = {1.5129, 1.4161, 0.5476, 0.0121, 6.5536};

    gsl_matrix *A;
    A = gsl_matrix_alloc(size, size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            double temp = pow(x[i],j);
            gsl_matrix_set(A,i,size - j -1,temp);
        }
    }

    gsl_vector *fx;
    fx = gsl_vector_alloc(size);
    for (int k = 0; k < size; ++k) {
        gsl_vector_set(fx, k, y[k]);
    }
    gsl_vector *a;
    a = gsl_vector_alloc(size);


    int s;

    gsl_permutation * p = gsl_permutation_alloc(size);

    gsl_linalg_LU_decomp (A, p, &s);

    gsl_linalg_LU_solve (A, p, fx, a);

    printf ("Wspolczynniki a = \n");
    gsl_vector_fprintf (stdout, a, "%g");


    return 0;
}