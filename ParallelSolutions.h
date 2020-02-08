#pragma once
#include "DataStorage.h"
#include "Parameters.h"
#include "ParallelSolvers.h"
#include <cmath>
#include <iostream>
#include "mpi.h"
#include "algorithm"
#include <assert.h>


const double PI = atan(1.0) * 4;

double ParallelAnalytical(DataStorage& storage, const Parameters params);
double FTCS(DataStorage &storage, const Parameters params); 
double laasonenImplicit(DataStorage& storage, const Parameters params);
double crankNicholson(DataStorage &storage, const Parameters params);
