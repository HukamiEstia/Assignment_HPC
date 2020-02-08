#pragma once
#include <cmath>
#include <iostream>
#include <vector>
#include "mpi.h"

#include "Utils.h"
using namespace std;

void ParallelThomas(int world_rank, int world_size,
                    int N, int *strides, int *offset,
                    double* a, double* b, double *c,
                    double* x, double* q);