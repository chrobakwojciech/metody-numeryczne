#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

void printG(vector<vector<double>> G) {
    cout << endl;
    for (int i = 0; i < G.size(); i++) {
        for (int j = 0; j < G.size(); j++) {
            cout << "[ " << fixed << setprecision(6) << showpos << G[i][j] << " ]  ";
        }
        cout << endl;
    }
    cout << endl;
}

void printA(vector<double> A) {
    for (int i = 0; i < A.size(); i++) {
        cout << "[ " << showpoint << A[i] << " ]  ";
    }
    cout << endl;
}

int main() {
    const int sizeMatrix = 7;

    vector<double> d_vector;
    d_vector.assign(sizeMatrix, 4);

    // wektor pod diagonala
    vector<double> e1_vector;
    e1_vector.assign(sizeMatrix - 1, 1);

    // wektor nad diagonala
    vector<double> e2_vector;
    e2_vector.assign(sizeMatrix - 1, 1);

    // wyrazy wolne
    vector<double> right = {1, 2, 3, 4, 5, 6, 7};

    // wynik
    vector<double> x;
    x.assign(sizeMatrix, 0);

    // wektor drugi nad diagonala
    vector<double> f_vector;
    f_vector.assign(sizeMatrix - 2, 0);


    cout << "Ogolna postac macierzy: " << endl;
    for (int i = 0; i < sizeMatrix; i++) {
        for (int j = 0; j < sizeMatrix; j++) {
            if (i == j)
                cout << " d  ";
            else if (j == i + 1)
                cout << "e2  ";
            else if (j == i - 1)
                cout << "e1  ";
            else
                cout << " 0  ";
        }
        cout << endl;
    }
    cout << endl;
    cout << "wektor D  = ";
    printA(d_vector);
    cout << "wektor E1 = ";
    printA(e1_vector);
    cout << "wektor E2 = ";
    printA(e2_vector);
    cout << endl;
    vector<vector<double>> G = {{1, 1},
                                {1, 1}};

    for (int i = 0; i < d_vector.size() - 1; i++) {
        // b element zerowany, a element nad nim
        double a = d_vector[i];
        double b = e1_vector[i];

        cout << "a = " << a << endl;
        cout << "b = " << b << endl;
        double cosX = a / sqrt(a * a + b * b);
        double sinX = -b / sqrt(a * a + b * b);
        cout << "cosX = " << cosX << endl;
        cout << "sinX = " << sinX << endl;
        G[0][0] = cosX;
        G[0][1] = -sinX;
        G[1][0] = sinX;
        G[1][1] = cosX;
        printG(G);

        vector<double> temp_a_diag = d_vector;
        vector<double> temp_a_pod = e1_vector;
        vector<double> temp_a_nad = e2_vector;
        vector<double> temp_a_nad2 = f_vector;

        // stosujemy obroty Givensa
        d_vector[i] = G[0][0] * temp_a_diag[i] + G[0][1] * temp_a_pod[i];
        d_vector[i + 1] = G[1][0] * temp_a_nad[i] + G[1][1] * temp_a_diag[i + 1];

        e2_vector[i] = G[0][0] * temp_a_nad[i] + G[0][1] * temp_a_diag[i+1];
        e2_vector[i + 1] = G[1][0] * temp_a_nad2[i] + G[1][1] * temp_a_nad[i + 1];

        f_vector[i] = G[0][0] * temp_a_nad2[i] + G[0][1] * temp_a_nad[i+1];

        e1_vector[i] = 0;

        cout << "wektor D  = ";
        printA(d_vector);
        cout << "wektor E1 = ";
        printA(e1_vector);
        cout << "wektor E2 = ";
        printA(e2_vector);
        cout << "wektor F  = ";
        printA(f_vector);

        // dzialamy macierza g na wektor wyrazow wolnych
        vector<double> temp_right = right;

        right[i] = G[0][0] * temp_right[i] + G[0][1] * temp_right[i + 1];
        right[i + 1] = G[1][0] * temp_right[i] + G[1][1] * temp_right[i + 1];


        cout << endl << endl << "---------------------------------------------------" << endl;
    }

    // rozwiazujemy Ax=b metoda backSubs
    cout << "wektor B  = ";
    printA(right);
    for (int r = d_vector.size() - 1; r >= 0; r--) {
        double val = 0;
        val = e2_vector[r] * x[r+1] + f_vector[r] * x[r+2];

        val = right[r] - val;
        x[r] = val / d_vector[r];
    }

    cout << "x  = ";
    printA(x);
    return 0;
}