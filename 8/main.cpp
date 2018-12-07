#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>

using namespace std;

const int SIZE = 65;

double lagrange(vector<double> x, vector<double> y, double val ){
    double t;
    double res = 0.0;

    for(int k = 0; k< x.size(); k++){
        t = 1.0;
        for(int j = 0; j < x.size() ; j++){
            if(j != k ){
                t=t*((val-x[j])/(x[k]-x[j]));
            }
        }
        res += t*y[k];
    }
    return res;
}

int main() {
    ofstream lgr;
    lgr.open("lagrange.dat");

    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(16);
    vector<double> x;
    vector<double> y;

    for (int i = 0; i < SIZE; ++i) {
        x.push_back(-1 + i / 32.0);
        y.push_back(1 / (1 + 5 * x[i] * x[i]));
    }

    // wyliczanie wartosci wielomianu Lagrange dla przedzialu -2 do 2
    for (double j = -2; j < 2; j = j + 0.01) {
        lgr << j << " " << lagrange(x,y,j) << endl;
    }

    lgr.close();
    return 0;
}