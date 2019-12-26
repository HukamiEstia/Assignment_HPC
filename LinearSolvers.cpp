#include "LinearSolvers.h"
#include "Utils.h"

void LUFactorisation(vector<vector<double>> A, vector<vector<double>>& L, vector<vector<double>>& U);
vector<double> LUSolve(vector<vector<double>> L, vector<vector<double>> U, vector<double> b);

vector<double> Thomas(vector<vector<double>> A, vector<double> b) {
	double m, a, c, d;
	// Forward Elimination
	for (int k = 1; k < A.size(); k++) {
		m = A[k][k - 1] / A[k - 1][k - 1];
		A[k][k] = A[k][k] - m * A[k - 1][k];
		b[k] = b[k] - m * b[k - 1];
	}
	
	//Backwards Substitution
	vector<double> x = b;

	for (int n = 0; n < x.size(); n++) {
		x[n] = b[n] / A[n][n];
	}

	for (int k = (A.size() - 2); k >= 0; k--) {
		x[k] = (b[k] - A[k][k + 1] * x[k + 1]) / A[k][k];
	}
	//printVector(x);

	return x;
}

vector<double> LUDecomp(vector<vector<double>> A, vector<double> b) {
	vector<vector<double>> L = getZeroSquareMatrix(A.size()), U = getZeroSquareMatrix(A.size());
	LUFactorisation(A, L, U);
	return LUSolve(L, U, b);
}

void LUFactorisation(vector<vector<double>> A, vector<vector<double>>& L, vector<vector<double>>& U) {
	vector<vector<double>> temp = A;
	double mult;
	int n = A.size();

	for (int k = 0; k < (n - 1); k++) {
		for (int i = k + 1; i < n; i++) {
			if (fabs(temp[k][k]) < 1.e-07) {
				cout << "Pivot is zero" << endl;
				exit(1);
			}
			mult = temp[i][k] / temp[k][k];
			temp[i][k] = mult;
			for (int j = k + 1; j < n; j++) {
				temp[i][j] -= mult * temp[k][j];
				if (fabs(temp[i][i]) < 1.e-07) {
					cout << "Pivot is zero!" << endl;
					exit(1);
				}
			}
		}
	}

	for (int i = 0; i < n; i++) {
		L[i][i] = 1.0;
	}
	for (int i = 1; i < n; i++) {
		for (int j = 0; j < i; j++) {
			L[i][j] = temp[i][j];
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			U[i][j] = temp[i][j];
		}
	}
}

vector<double> LUSolve(vector<vector<double>> L, vector<vector<double>> U, vector<double> b) {
	vector<double> temp = b;
	int n = L.size();

	// Forward substitution
	for (int i = 1; i < n; i++) {
		for (int j = 0; j < i; j++) {
			temp[i] -= L[i][j] * temp[j];
		}
	}

	// Backward substitution
	for (int i = n - 2; i > 0; i--) {
		for (int j = i + 1; j < n; j++) {
			temp[i] -= U[i][j] * temp[j];
		}
		temp[i] /= U[i][i];
	}
	return temp;
}
