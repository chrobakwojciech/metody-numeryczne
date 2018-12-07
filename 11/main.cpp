#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>

using namespace std;

double f(double x) {
    return (x * x - 1) * pow(sinh(x), 3);
}

double lagrange(vector<double> x, vector<double> y, double val) {
    double t;
    double res = 0.0;

    for (int k = 0; k < x.size(); k++) {
        t = 1.0;
        for (int j = 0; j < x.size(); j++) {
            if (j != k) {
                t = t * ((val - x[j]) / (x[k] - x[j]));
            }
        }
        res += t * y[k];
    }
    return res;
}

vector<double> interpolacjaOdwrotna(double x1, double x2, double x3) {
    double x0 = x1;
    vector<double> x, y;
    x.resize(3);
    y.resize(3);

    vector<double> res;

    int counter = 1;
    while (true) {
        x[0] = f(x1);
        x[1] = f(x2);
        x[2] = f(x3);

        y[0] = x1;
        y[1] = x2;
        y[2] = x3;
        x0 = lagrange(x, y, 0.0);

        x1 = x2;
        x2 = x3;
        x3 = x0;

        if(isnan(x0)){
            cout << "[INTERPOLACJA ODWROTNA] Dzielenie przez 0" << endl;
            break;
        }
        if (abs(x0) < 1e-8)
            break;

        counter++;
    }


    res.push_back(x0);
    res.push_back((double) counter);

    return res;
}

vector<double> sieczne(double x1, double x2) {
    double x0;
    vector<double> res;

    if (abs(f(x1) - f(x2)) < 1e-8) {
        cout << "Wartosci w punktach startowych sa takie same! Blad.\n";
        return res;
    }

    double xm;
    int counter = 1;
    do {
        x0 = (x1 * f(x2) - x2 * f(x1)) / (f(x2) - f(x1));

        double c = f(x1) * f(x0);

        x1 = x2;
        x2 = x0;

        if (c == 0) {
            cout << "[SIECZNE]               Dzielenie przez 0" << endl;
            return res;
        }

        xm = (x1 * f(x2) - x2 * f(x1)) / (f(x2) - f(x1));

        counter++;

    } while (abs(xm - x0) >= 1e-8);


    res.push_back(x0);
    res.push_back((double) counter);

    return res;
}

int main() {
    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(16);
    srand(time(NULL));

    double x1, x2, x3;
    vector<double> temp;
    vector<vector<double>> punkty;
    for (int j = 0; j < 10; ++j) {
        temp.clear();

        x1 = (rand() % 100) / 100.0;
        do {
            x2 = (rand() % 100) / 100.0;
        } while(x2 == x1);

        do {
            x3 = (rand() % 100) / 100.0;
        } while(x3 == x2);


        temp.push_back(x1);
        temp.push_back(x2);
        temp.push_back(x3);

        punkty.push_back(temp);
    }

    for (int i = 0; i < 10; ++i) {
        double t = sieczne(punkty[i][0], punkty[i][1])[0];
        double k = sieczne(punkty[i][0], punkty[i][1])[1];
        cout << setprecision(2) << punkty[i][0] << " & " << punkty[i][1] << " & " << setprecision(16)<< t << " = "  << setprecision(2) << t << setprecision(0) << " & " << k << " \\\\" << endl;
    }

    cout << endl << endl;
    for (int i = 0; i < 10; ++i) {
        double t = interpolacjaOdwrotna(punkty[i][0], punkty[i][1], punkty[i][2])[0];
        double k = interpolacjaOdwrotna(punkty[i][0], punkty[i][1], punkty[i][2])[1];
        cout << setprecision(2) << punkty[i][0] << " & " << punkty[i][1] << " & " << punkty[i][2] << " & " << setprecision(16) << t << " = " << setprecision(2) <<  t << setprecision(0) << " & " << k << " \\\\" << endl;

    }

    return 0;
}