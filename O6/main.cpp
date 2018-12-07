#include <vector>
#include <iostream>
#include <iomanip>
#include <complex>

using namespace std;

typedef complex<double> zesp;

const double eps = 1e-8;

zesp horner(const vector<zesp> &P, zesp x) {
    int n = P.size() - 1;
    zesp res = P[n];
    for (int i = n - 1; i >= 0; i--) {
        res = res * x + P[i];
    }

    return res;
}

vector<zesp> pochodna(const vector<zesp> &P) {
    vector<zesp> result;
    result.resize(P.size() - 1);

    int n = P.size();
    for (double i = 1; i < n; ++i) {
        zesp temp = i;
        result[i - 1] = temp * P[i];
    }

    return result;
}

vector<zesp> deflacja(const vector<zesp> &P, zesp z0) {
    vector<zesp> result;
    result.resize(P.size() - 1);
    result[P.size() - 2] = P[P.size() - 1];
    for (int i = P.size() - 3; i >= 0; i--) {
        result[i] = P[i + 1] + result[i + 1] * z0;
    }

    return result;
}

zesp Laguerre(const vector<zesp> &P, zesp z0) {
    int n = P.size() - 1;
    zesp nz = zesp((double) n);
    vector<zesp> pochodna1 = pochodna(P);
    vector<zesp> pochodna2 = pochodna(pochodna1);

    for (int i = 0; i < 500000; ++i) {
        zesp wartosc = horner(P, z0);

        if (abs(wartosc) < eps)
            break;

        zesp wartosc_pochodna1 = horner(pochodna1, z0);
        zesp wartosc_pochodna2 = horner(pochodna2, z0);

        zesp G = wartosc_pochodna1 / wartosc;
        zesp H = G * G - wartosc_pochodna2 / wartosc;
        zesp R = sqrt((nz - 1.0) * (H * nz - G * G));

        zesp D1 = G + R;
        zesp D2 = G - R;

        zesp M;
        if (abs(D1) > abs(D2))
            M = D1;
        else
            M = D2;

        zesp a = nz / M;

        z0 -= a;


        if (abs(a) < eps)
            break;
    }

    return z0;
}

vector<zesp> szukaj(const vector<zesp> &P) {
    vector<zesp> res;
    vector<zesp> P_def = P;
    while (P_def.size() > 2) {
        zesp z = (rand() / double(RAND_MAX), rand() / double(RAND_MAX));
        z = Laguerre(P_def, z);
        z = Laguerre(P, z); // wygladzanie
        P_def = deflacja(P_def, z); // faktoryzacja
        res.push_back(z);
    }
    res.push_back(-P_def[0] / P_def[1]);
    return res;
}

void wypisz(vector<zesp> roots) {
    for (auto x : roots) {
        if (abs(x.imag()) < 1e-8)
            x.imag(0);
        if (abs(x.real()) < 1e-8)
            x.real(0);
        cout << x << endl;
    }
    cout << endl;
}

int main() {
    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(8);

    vector<zesp> P1{16.0, -72.0, -28.0, 558.0, -990.0, 783.0, -486.0, 243.0};
    vector<zesp> P2{-4.0, -4.0, -12.0, -8.0, -11.0, -3.0, -1.0, 2.0, 3.0, 1.0, 1.0};
    vector<zesp> P3{1.0, 0.0, -1.0, 0.0, 1.0};
    P3[1].imag(-1.0);
    P3[3].imag(1.0);

    vector<zesp> roots1 = szukaj(P1);
    wypisz(roots1);

    vector<zesp> roots2 = szukaj(P2);
    wypisz(roots2);

    vector<zesp> roots3 = szukaj(P3);
    wypisz(roots3);

    return 0;
}