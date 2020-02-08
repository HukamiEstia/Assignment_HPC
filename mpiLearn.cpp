#include <mpi.h>
#include <iostream>
#include <vector>
#include <algorithm>

int main() {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int* vec = NULL;
    if (world_rank == 0){
        vec = new int[8];
        for (int i = 0; i < 8; i++){
            vec[i] = i;
        }
    }
    
    unsigned int subRange = 8 / world_size; 
	unsigned int subStart = world_rank * subRange;
	unsigned int subEnd;

    int res[2];

    MPI_Scatter(vec, 2, MPI_INT, &res, 2, MPI_INT, 0, MPI_COMM_WORLD);

    for (int i = 0; i < 2; i++){
        res[i] += 1;
    }

    int* recVec = NULL;
    if (world_rank == 0){
        recVec = new int[8];
    }

    MPI_Gather(&res, 2, MPI_INT, recVec, 2, MPI_INT, 0, MPI_COMM_WORLD);

    if (world_rank == 0){
        for (int i=0; i < 8; i++){
            std::cout << "test\n";
            std::cout << recVec[i] << " ";
        }
    }

    // Finalize the MPI environment.
    MPI_Finalize();
}


    

    