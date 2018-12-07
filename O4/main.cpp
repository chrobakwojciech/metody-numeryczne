#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <fstream>

using namespace std;

int main() {
    ofstream file, file2;
    file2.open("input.dat"); // plik do wezlow i jego wartosci
    file.open("splines.dat"); // plik do wartosci splajnu
    int i;
    double xi, yi, x[65], y[65];

    // wypelnianie tablic wezlow i jego wartosci
    for (i = 0; i < 65; i++) {
        x[i] = -1 + i / 32.0;
        y[i] = 1 / (1 + 5 * x[i] * x[i]);
        file2 << x[i] << " " << y[i] << endl;
    }

    cout << endl << endl;
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, 65);

    gsl_spline_init(spline, x, y, 65);

    // obliczanie wartosci splajnu
    for (double j = -1.0; j <= 1.0; j = j + 0.01) {
        yi = gsl_spline_eval(spline, j, acc);
        file << j << " " << yi << endl;
    }

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    file.close();
    file2.close();
    return 0;
}