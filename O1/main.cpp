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
        cout <<  A[i] << ", ";
    }
    cout << endl;
}

double multiVectorVector(vector<double> a, vector<double> b) {
    double sum = 0;
    for(int i = 0; i < a.size(); i++) {
        sum = sum + a[i] * b[i];
    }

    return sum;
}


int main() {
    const int sizeMatrix = 7;
    vector<double> u_vector = {1, 0, 0, 0, 0, 0, 1};
    vector<double> v_vector = {1, 0, 0, 0, 0, 0, 1};

    vector<double> d_vector;
    d_vector.assign(sizeMatrix, 4);

    // A1 = A + uvT
    d_vector[0] -= u_vector[0];
    d_vector[sizeMatrix-1] -= v_vector[sizeMatrix-1];


    // wektor pod diagonala
    vector<double> e1_vector;
    e1_vector.assign(sizeMatrix - 1, 1);

    // wektor nad diagonala
    vector<double> e2_vector;
    e2_vector.assign(sizeMatrix - 1, 1);

    // wyrazy wolne
    vector<double> right = {1, 2, 3, 4, 5, 6, 7};

    // z
    vector<double> z;
    z.assign(sizeMatrix, 0);
    // 1
    vector<double> q;
    q.assign(sizeMatrix, 0);
    // x
    vector<double> x;
    x.assign(sizeMatrix, 0);

    // wektor drugi nad diagonala
    vector<double> f_vector;
    f_vector.assign(sizeMatrix - 2, 0);

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


    for (int i = 0; i < d_vector.size()-1; i++) {
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

        e2_vector[i] = G[0][0] * temp_a_nad[i] + G[0][1] * temp_a_diag[i + 1];
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

        // dzialamy macierza g na wektor u
        vector<double> temp_u = u_vector;

        u_vector[i] = G[0][0] * temp_u[i] + G[0][1] * temp_u[i + 1];
        u_vector[i + 1] = G[1][0] * temp_u[i] + G[1][1] * temp_u[i + 1];


        cout << endl << endl << "---------------------------------------------------" << endl;
    }


    // rozwiazujemy Az=b metoda backSubs
    cout << "b = ";
    printA(right);
    for (int r = d_vector.size() - 1; r >= 0; r--) {
        double val = 0;
        val = e2_vector[r] * z[r+1] + f_vector[r] * z[r+2];

        val = right[r] - val;
        z[r] = val / d_vector[r];
    }

    cout << "z = ";
    printA(z);

//     rozwiazujemy Aq=u metoda backSubs
    cout << "u  = ";
    printA(u_vector);
    for (int r = d_vector.size() - 1; r >= 0; r--) {
        double val = 0;
        val = e2_vector[r] * q[r+1] + f_vector[r] * q[r+2];

        val = u_vector[r] - val;
        q[r] = val / d_vector[r];
    }

    cout << "q = ";
    printA(q);

    double l = multiVectorVector(v_vector,z);
    double m = 1 + multiVectorVector(v_vector,q);
    double temp = l/m;

    vector<double> temp2;
    temp2.assign(sizeMatrix, 0);

    for(int i = 0; i < sizeMatrix; i++) {
        temp2[i] = temp * q[i];
    }
    printA(temp2);
    for(int i = 0; i < sizeMatrix; i++) {
        x[i] = z[i] - temp2[i];
    }

    cout << "x = ";
    printA(x);

    return 0;
}