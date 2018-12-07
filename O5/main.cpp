#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>

using namespace std;

double f(double x) {
    return sin(M_PI*(1+sqrt(x))/(1+x*x))*exp(-x);
}

int findMax () {
    int B = 0;
    double res = exp(-B);
    double res_prev = res;
    B++;
    while(true) {
        res_prev = res;
        res = exp(-B);
        if(res_prev - res < 1e-7)
            break;
        B++;
    }
    return B;
}


int main() {
    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(10);
    int A = 0;
    int B = findMax();

    // METODA ROMBERGA
    vector<vector<double>> R;
    double h = B-A;

    R.push_back(vector<double>());

    // liczymy R[0][0]
    R[0].push_back(h/2.0 * (f(A)+f(B)));

    double sumF = 0.5 * f(A) + 0.5 * f(B);

    int row = 1;
    while(true) {
        R.push_back(vector<double>());
        h = h/2.0;

        // liczymy pierwsza kolumne
        for (int k = 1; k <= pow(2,row-1); ++k) {
            sumF += f(A+(2*k-1)*h);
        }
        R[row].push_back(sumF * h);

        // liczymy kolejne wartosci w wierszach
        for (int col = 1; col <= row; ++col) {
            R[row].push_back(((pow(4,col)*R[row][col-1] - R[row-1][col-1]))/(pow(4,col)-1));
        }

        // sprawdzamy czy koniec algorytmu
        if(abs(R[row-1][row-1] - R[row][row]) < 1e-7) {
            break;
        }

        row++;
    }

    for (int i = 0; i < R.size(); ++i) {
        cout << R[i][R[i].size()-1] << "\\\\" << endl;
    }

    // METODA TRAPEZOW
    vector<double> T;
    sumF = 0.5 * f(A) + 0.5 * f(B);
    h = B - A;
    T.push_back(h/2.0 * (f(A)+f(B)));
    int counter = 1;
    while(true) {
        h = h/2.0;
        for (int k = 1; k <= pow(2,counter-1); ++k) {
            sumF += f(A+(2*k-1)*h);
        }
        T.push_back(sumF * h);

        if(abs(T[counter-1] - T[counter]) < 1e-7) {
            break;
        }
        counter++;
    }

    for (int j = 0; j < T.size(); ++j) {
        cout << T[j] << "\\\\" << endl;
    }

    return 0;
}