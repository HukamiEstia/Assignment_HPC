#pragma once
#include <vector>
#include <iostream>
#include <iomanip>
using namespace std;

void printMatrix(vector<vector<double>> matrix);
vector<vector<double>> getZeroSquareMatrix(int size);
void printVector(vector<double> vec);
vector<double> dotProduct(vector<double> vec, vector<vector<double>> mat);
