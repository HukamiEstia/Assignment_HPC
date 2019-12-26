#pragma once
#include "DataStorage.h"
#include "SetVariables.h"
#include "LinearSolvers.h"
#include <cmath>
#include <iostream>


const double PI = atan(1.0) * 4;

void duFortFrankel(DataStorage &storage, const SetVariables variables);
void richardson(DataStorage &storage, const SetVariables variables); 
void laasonenImplicit(DataStorage& storage, const SetVariables variables, vector<double>(*f)(vector<vector<double>>, vector<double>));
void crankNicholson(DataStorage &storage, const SetVariables variables, vector<double>(*f)(vector<vector<double>>, vector<double>));
void analytical(DataStorage& storage, const SetVariables variables);
