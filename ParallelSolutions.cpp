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
		//Initializing vector holding the solution at t=n+1 for each process
		double subCurrent[subRange];

		for (unsigned int i = 0; i < subRange; i++) {
			double sum{ 0 };
			for (int k = 1; k < 200; k++) {
				//std::cout << "point is: " << storage.getDX() * (subStart + i) << "\n";
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
			vector<double> temp;
			for (int i = 0; i < numberOfX; i++){
				temp.push_back(currentTime[i]);
			}
			solutionMat.push_back(temp);
		}
	}
	if (world_rank == 0){
		for (auto e : solutionMat[10]){
				std::cout << e << " ";
			}
		storage.setData(solutionMat);
	}
    // Finalize the MPI environment.
    MPI_Finalize();
}

void FTCS(DataStorage &storage, const Parameters params) {
	/*
	This method is forward time central space
	/**/
	double left, right, central;
	
	//set dT limit so that the scheme is stable
	double dT = pow(storage.getDX(), 2) / params.getD();
	// std::cout << dT << "\n";
	//check whether scheme is stable and inform the user
	if (storage.getDT() < dT){
		std::cout << "Scheme is stable \n";
	} else {
		std::cout << "Scheme is unstable \n";
	}

	// Initialize the MPI environment
    MPI_Init(NULL, NULL);
	MPI_Status status;

	// Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	//compute the number of time and space discretisation steps
	unsigned int numberOfT = params.getDuration() / storage.getDT() + 1;
	unsigned int numberOfX = params.getL() / storage.getDX() + 1;

	//define r = D.(dt/dx^2)
	double r = params.getD() * (storage.getDT() / pow(storage.getDX(), 2));

	//compute domain division
	unsigned int subRange = (numberOfX / world_size);
	unsigned int remainder = numberOfX - (subRange * world_size); 

	int displacements[world_size];
	int strides[world_size];
	for (int i = 0; i < world_size; i++){
		strides[i] = subRange;
	}
	for (int i = 0; i < remainder; i++){
		strides[i]++;
	}

	int offset{0};
	for (int i = 0; i < world_size; i++) {
		displacements[i] = offset;
		offset += strides[i];
	}

	//Initialiazing the buffer to scatter the initials conditions
	double *sendBuf = NULL;
	if (world_rank == 0){
		vector<double> initialC = params.getTime_0(numberOfX);
		sendBuf = new double[numberOfX];
		for (int  i = 0; i < numberOfX; i++){
			sendBuf[i] = initialC[i];
		}
	}
	
	//Initializing matrix holding the solution
	double subMat[numberOfT][strides[world_rank]];

	//Scatter initials condition to each process	
	MPI_Datatype rtype;
	MPI_Type_vector(strides[world_rank], 1, 1, MPI_DOUBLE, &rtype);
	MPI_Type_commit(&rtype);
	double* rptr = &subMat[0][0];
	MPI_Scatterv(sendBuf, strides, displacements, MPI_DOUBLE,
				subMat, 1, rtype, 0, MPI_COMM_WORLD);
	
	for (unsigned int n = 1; n < numberOfT; n++) {
		double *subBoundaryLeft = NULL;
		double *subBoundaryRight = NULL;

		//set boundary conditions or compute
		//the internal boundary points
		
		//boundary conditions
		if (world_rank == 0){
			subMat[n][0] = params.getT_sur();
		} else if (world_rank == world_size - 1){
			subMat[n][strides[world_rank] - 1] = params.getT_sur();
		}

		//send boundary points to neigbour processes and 
		//getting data from neighbour processes
		if (world_rank < world_size - 1){
			//sending right boundary point
			double sendRightBoundary = subMat[n - 1][strides[world_rank] - 1];
			MPI_Send(&sendRightBoundary , 1, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD);
		}
		if (world_rank > 0){
			//receiving left boundary point
			double subBoundaryLeft;
			MPI_Recv(&subBoundaryLeft, 1, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			//computing first point
			//T(n,i-1)
			left = subBoundaryLeft;
			//T(n,i)
			central = subMat[n - 1][0];
			//T(n,i+1)
			right = subMat[n - 1][1];
			//compute T(n+1,i) = T(n,i) + r * (T(n,i+1) - 2 * T(n,i) + T(n,i-1))
			subMat[n][0] = central+r*(left-2*central+right);
		}
		if (world_rank > 0){
			//sending left boundary point
			double sendLeftBoundary = subMat[n - 1][0];
			// if (n == 1){
				// std::cout << "process: " << world_rank << " send: " << sendLeftBoundary << "\n";			
			// }
			MPI_Send(&sendLeftBoundary , 1, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD);
		}
		if (world_rank < world_size - 1){
			// receiving right boundary point
			double subBoundaryRight;
			MPI_Recv(&subBoundaryRight, 1, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			
			//computing last point
			//T(n,i-1)
			left = subMat[n - 1][strides[world_rank] - 2];
			//T(n,i)
			central = subMat[n - 1][strides[world_rank] - 1];
			//T(n,i+1)
			right = subBoundaryRight;
			//compute T(n+1,i) = T(n,i) + r * (T(n,i+1) - 2 * T(n,i) + T(n,i-1))
			subMat[n][strides[world_rank] - 1] = central+r*(left-2*central+right);
		}
		
		//computing process' sub domain
		for (unsigned int i = 1; i < strides[world_rank] - 1; i++) {
			//T(n,i-1)
			left = subMat[n - 1][i - 1];
			//T(n,i)
			central = subMat[n - 1][i];
			//T(n,i+1)
			right = subMat[n - 1][i + 1];
			//compute T(n+1,i) = T(n,i) + r * (T(n,i+1) - 2 * T(n,i) + T(n,i-1))
			subMat[n][i] = central+r*(left-2*central+right);
		}
		
		// if (world_rank == 0){
		// 	for (int i = 0; i < strides[world_rank]; i++){
		// 		std::cout << subMat[n][i] << " ";
		// 	}
		// 	std::cout << "\n";
		// 	// std::cout << "process " << world_rank << " " << subMat[n][0] << " " << subMat[n][strides[world_rank] - 1] << "\n"; 
		// }
	}
	
	
	vector<vector<double>> solutionMat;
	double *recvBuf = NULL;
	if (world_rank == 0) {
		recvBuf = new double[numberOfX];
	}
	MPI_Datatype stype;
	MPI_Type_vector(strides[world_rank], 1, 1, MPI_DOUBLE, &stype);
	MPI_Type_commit(&stype);
	for (int n = 0; n < numberOfT; n++){
		double* sptr = &subMat[n][0];
		//gather all subdomains computed solution
		MPI_Gatherv(sptr, 1, stype, recvBuf, strides, displacements,
			        MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		//store it in solution matrix
		if (world_rank == 0){
			vector<double> temp;
			for (int i = 0; i < numberOfX; i++){
				temp.push_back(recvBuf[i]);
			}
			solutionMat.push_back(temp);
		}
	}

	//store results
	if (world_rank == 0){
		std::cout << solutionMat.size() << "\n";
		for (auto e : solutionMat[4]){
				std::cout << e << " ";
			}
		storage.setData(solutionMat);
	}
    
	// Finalize the MPI environment.
    MPI_Finalize();
}