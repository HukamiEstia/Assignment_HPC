#include "Utils.h"

void printMatrix(vector<vector<double>> matrix) {
	for (int i = 0; i < matrix.size(); i++) {
		for (int j = 0; j < matrix[i].size(); j++) {
			cout << left << setw(8) << matrix[i][j] << " ";
		}
		cout << endl;
	}
}

void printVector(vector<double> vec) {
	for (int i = 0; i < vec.size(); i++) {
		cout << setprecision(4) <<  vec[i] << " ";
	}
	cout << endl;
}

vector<vector<double>> getZeroSquareMatrix(int size) {
	vector<vector<double>> matrix(size);
	for (int i = 0; i < size; i++) {
		matrix[i] = vector<double>(size);
	}
	return matrix;
}

vector<double> dotProduct(vector<double> vec, vector<vector<double>> mat) {
	vector<double> res(vec.size());
	for (int i = 0; i < vec.size(); i++) {
		for (int j = 0; j < vec.size(); j++) {
			res[i] += mat[i][j] * vec[j];
		}
	}
	return res;
}