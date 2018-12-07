#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <iomanip>

using namespace std;

const int size = 8;

int main() {
    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(4);

    double x[] = {0.062500, 0.187500, 0.312500, 0.437500, 0.562500, 0.687500, 0.812500, 0.935700};
    double y[] = {0.687959, 0.073443, -0.517558, -1.077264, -1.600455, -2.080815, -2.507266, -2.860307};

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

    for (int l = 0; l < size; ++l) {
        cout << gsl_vector_get(a,l) << endl;
    }


    return 0;
}