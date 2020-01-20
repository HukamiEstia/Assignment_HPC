#include "mpi.h"
#include "ParallelSolutions.h"
#include "algorithm"

void ParallelAnalytical(DataStorage &storage, const Parameters params) {
	/*
	Compute the analytical solution for the specified discretization step and problem parameters
	/**/
	
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
	MPI_Status status;

	// Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	unsigned int numberOfT = (params.getDuration() / (double)storage.getDT()) + 1;
	unsigned int numberOfX = ((double)params.getL() / (double)storage.getDX()) + 1;
	
	unsigned int subRange = (numberOfX / world_size) + 1;
	unsigned int subStart = world_rank * subRange;
	unsigned int subEnd = std::min(subStart + subRange, numberOfX);

	double *currentTime = NULL;
	vector<vector<double>> solutionMat;
	if (world_rank == 0){
		solutionMat.push_back(params.getTime_0(numberOfX));
		currentTime = new double[subRange * world_size];
	}

	for (unsigned int n = 1; n < numberOfT; n++) {

		double subCurrent[subRange];
		for (unsigned int i = 0; i < subRange; i++) {
			double sum{ 0 };
			for (int k = 1; k < 200; k++) {
				std::cout << "point is: " << storage.getDX() * (subStart + i) << "\n";
				sum += exp(-params.getD()
						* powf((k * PI / params.getL()), 2.0)
						* storage.getDT() * n) 
					* (1 - (pow(-1, k))) 
					* sin(k * PI * storage.getDX() * (subStart + i) / params.getL())
					/ (k * PI);
			}
			subCurrent[i] = params.getT_sur() + 2 * (params.getT_in() - params.getT_sur())* sum;
		}
		
		MPI_Gather(&subCurrent, subRange, MPI_DOUBLE,
				   currentTime, subRange, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if (world_rank == 0){
			vector<double> temp[numberOfX];
			for (int i = 0; i < numberOfX; i++){
				temp[i] = currentTime[i];
			}
			/*
			solutionMat.push_back(temp);
			std::cout << "time: " << n << "\n";
			for (auto e : solutionMat[n]){
				std::cout << e << " ";
			}
			std::cout << "\n\n\n";
			*/
		}
	}

	//storage.setData(solutionMat);

    // Finalize the MPI environment.
    MPI_Finalize();
}