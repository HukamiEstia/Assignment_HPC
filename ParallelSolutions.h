#pragma once
#include "DataStorage.h"
#include "Parameters.h"
#include "LinearSolvers.h"
#include <cmath>
#include <iostream>


const double PI = atan(1.0) * 4;

void ParallelAnalytical(DataStorage& storage, const Parameters params);
void FTCS(DataStorage &storage, const Parameters params); 
void laasonenImplicit(DataStorage& storage, const Parameters params, vector<double>(*f)(vector<vector<double>>, vector<double>));
void crankNicholson(DataStorage &storage, const Parameters params, vector<double>(*f)(vector<vector<double>>, vector<double>));
