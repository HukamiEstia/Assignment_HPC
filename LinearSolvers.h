#pragma once
#include <cmath>
#include <iostream>
#include <vector>

#include "Utils.h"
using namespace std;

vector<double> Thomas(vector<vector<double>> A, vector<double> b);
vector<double> LUDecomp(vector<vector<double>> A, vector<double> b);
