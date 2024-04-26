#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Function to calculate matrix inverse
vector<vector<double>> inverseMatrix(vector<vector<double>> A) {
    // Calculate determinant of A
    double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    // Check if inverse exists
    if (det == 0) {
        cout << "Matrix is singular, inverse does not exist." << endl;
        exit(1);
    }
    // Calculate inverse
    double invDet = 1.0 / det;
    vector<vector<double>> inverse(2, vector<double>(2));
    inverse[0][0] = A[1][1] * invDet;
    inverse[0][1] = -A[0][1] * invDet;
    inverse[1][0] = -A[1][0] * invDet;
    inverse[1][1] = A[0][0] * invDet;
    return inverse;
}

// Function to solve linear equations using LU decomposition
vector<double> solveLU(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<vector<double>> U(n, vector<double>(n, 0));

    // LU decomposition
    for (int i = 0; i < n; i++) {
        L[i][i] = 1;
        for (int j = i; j < n; j++) {
            double sum = 0;
            for (int k = 0; k < i; k++) {
                sum += L[i][k] * U[k][j];
            }
            U[i][j] = A[i][j] - sum;
        }
        for (int j = i + 1; j < n; j++) {
            double sum = 0;
            for (int k = 0; k < i; k++) {
                sum += L[j][k] * U[k][i];
            }
            L[j][i] = (A[j][i] - sum) / U[i][i];
        }
    }

    // Solve LY = b
    vector<double> Y(n);
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * Y[j];
        }
        Y[i] = (b[i] - sum) / L[i][i];
    }

    // Solve UX = Y
    vector<double> X(n);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += U[i][j] * X[j];
        }
        X[i] = (Y[i] - sum) / U[i][i];
    }
    return X;
}

// Function to solve linear equations using Gauss elimination
vector<double> solveGauss(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    // Forward Elimination
    for (int i = 0; i < n - 1; i++) {
        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }
    // Back Substitution
    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
    return x;
}

// Function to solve linear equations using Crout decomposition
vector<double> solveCrout(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<vector<double>> U(n, vector<double>(n, 0));

    // Crout decomposition
    for (int i = 0; i < n; i++) {
        U[i][i] = 1;
        for (int j = i; j < n; j++) {
            double sum = 0;
            for (int k = 0; k < i; k++) {
                sum += L[i][k] * U[k][j];
            }
            L[i][j] = A[i][j] - sum;
        }
        for (int j = i + 1; j < n; j++) {
            double sum = 0;
            for (int k = 0; k < i; k++) {
                sum += L[j][k] * U[k][i];
            }
            U[j][i] = (A[j][i] - sum) / L[i][i];
        }
    }

    // Solve LY = b
    vector<double> Y(n);
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * Y[j];
        }
        Y[i] = (b[i] - sum) / L[i][i];
    }

    // Solve UX = Y
    vector<double> X(n);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += U[i][j] * X[j];
        }
        X[i] = (Y[i] - sum) / U[i][i];
    }
    return X;
}
void test() {
    // Test matrices
    vector<vector<double>> A = {{2, 3}, {1, -2}};
    vector<double> b = {8, -4};
    // Test methods
    vector<double> x_inverse = solveGauss(A, b);
    vector<double> x_LU = solveLU(A, b);
    vector<double> x_Gauss = solveGauss(A, b);
    vector<double> x_Crout = solveCrout(A, b);
    // Output solutions
    cout << "Solution using matrix inverse method: ";
    for (double xi : x_inverse) cout << xi << " ";
    cout << endl;
    cout << "Solution using LU decomposition method: ";
    for (double xi : x_LU) cout << xi << " ";
    cout << endl;
    cout << "Solution using Gauss elimination method: ";
    for (double xi : x_Gauss) cout << xi << " ";
    cout << endl;
    cout << "Solution using Crout decomposition method: ";
    for (double xi : x_Crout) cout << xi << " ";
    cout << endl;
}

int main() {
    test(); // Run test
    return 0;
}
