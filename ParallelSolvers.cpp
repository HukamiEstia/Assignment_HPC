#include "ParallelSolvers.h"

void ParallelThomas(int world_rank, int world_size,
                    int N, int *strides, int *offset,
                    double* a, double* b, double *c,
                    double* x, double* q){
    int i, j, k, i_global;
    double S[2][2], T[2][2], s1tmp, s2tmp;
    double *l, *d, *y;
    MPI_Status status;

    l = new double[N];
    d = new double[N];
    y = new double[N];

    for (i=0; i<N; i++){
        l[i] = d[i] = y[i] = 0.0;
    }

    //Identity matrix
    S[0][0] = S[1][1] = 1.0;
    S[1][0] = S[0][1] = 0.0;

    //form local products of R_k matrices
    if (world_rank == 0 ) {
        //form R_0
        S[0][0] = a[offset[world_rank]];
        S[0][1] = 0.0;
        S[1][0] = 1.0;
        S[1][1] = 0.0;
        for (int i = 1; i < strides[world_rank]; i++){
            // std::cout << "stride: " << strides[world_rank] << "\n";
            //multiply current S by R_i
            s1tmp = a[i + offset[world_rank]] * S[0][0] -
                    b[i + offset[world_rank] - 1] *
                    c[i + offset[world_rank] - 1] * S[1][0];
            s2tmp = a[i + offset[world_rank]] * S[0][1] -
                    b[i + offset[world_rank] - 1] *
                    c[i + offset[world_rank] - 1] * S[1][1];
            S[1][0] = S[0][0];
            S[1][1] = S[0][1];
            S[0][0] = s1tmp;
            S[0][1] = s2tmp;
            // std::cout << "i: " << i << "\n";
            // std::cout << S[0][0] << " " << S[0][1] << "\n";
            // std::cout << S[1][0] << " " << S[1][1] << "\n";
        }
    } else{
        for (int i = 0; i < strides[world_rank]; i++){
            s1tmp = a[i + offset[world_rank]] * S[0][0] -
                    b[i + offset[world_rank] - 1] * c[i + offset[world_rank] - 1] * S[1][0];
            s2tmp = a[i + offset[world_rank]] * S[0][1] -
                    b[i + offset[world_rank] - 1] * c[i + offset[world_rank] - 1] * S[1][1];
            S[1][0] = S[0][0];
            S[1][1] = S[0][1];
            S[0][0] = s1tmp;
            S[0][1] = s2tmp;
        }
    }
    for (int i=0; i<=log2(world_size); i++){
        if (world_rank+pow(2, i) < world_size){
            MPI_Send(S, 4, MPI_DOUBLE, int(world_rank+pow(2,i)), 0,
                    MPI_COMM_WORLD);
        }
        if(world_rank-pow(2,i) >= 0){
            MPI_Recv(T, 4, MPI_DOUBLE, int(world_rank-pow(2,i)), 0,
                    MPI_COMM_WORLD, &status);
            s1tmp = S[0][0] * T[0][0] + S[0][1] * T[1][0];
            S[0][1] = S[0][0] * T[0][1] + S[0][1] * T[1][1];
            S[0][0] = s1tmp;
            s1tmp = S[1][0] * T[0][0] + S[1][1] * T[1][0];
            S[1][1] = S[1][0] * T[0][1] + S[1][1] * T[1][1];
            S[1][0] = s1tmp;
        }
    }

    d[offset[world_rank] + strides[world_rank] - 1] = (S[0][0] + S[0][1])/(S[1][0] + S[1][1]);

    if (world_rank == 0){
        MPI_Send(&d[offset[world_rank] + strides[world_rank] - 1], 1, MPI_DOUBLE,
                 1, 0, MPI_COMM_WORLD);
    } else{
        MPI_Recv(&d[offset[world_rank] - 1], 1, MPI_DOUBLE,
                 world_rank - 1, 0, MPI_COMM_WORLD, &status);
        if (world_rank != world_size - 1){
            MPI_Send(&d[offset[world_rank] + strides[world_rank] - 1], 1, MPI_DOUBLE,
                     world_rank + 1, 0, MPI_COMM_WORLD);
        }
    }

    if (world_rank == 0){
        l[0] = 0;
        d[0] = a[0];
        for (int i=1; i<strides[world_rank] - 1; i++){
            l[offset[world_rank] + i] = b[offset[world_rank] + i]/
                                        d[offset[world_rank] + i-1];
            d[offset[world_rank] + i] = a[offset[world_rank] + i] -
                                        l[offset[world_rank] + i] * c[offset[world_rank] + i];
        }
        l[offset[world_rank] + strides[world_rank] - 1] = b[offset[world_rank] + strides[world_rank] - 1]/
                                                          d[offset[world_rank] + strides[world_rank] - 2];
    } else{
        for (int i=0; i<strides[world_rank] - 1; i++){
            l[offset[world_rank] + i] = b[offset[world_rank] + i]/
                                        d[offset[world_rank] + i-1];
            d[offset[world_rank] + i] = a[offset[world_rank] + i] -
                                        l[offset[world_rank] + i] * 
                                        c[offset[world_rank] + i];
        }
        l[offset[world_rank] + strides[world_rank] - 1] = b[offset[world_rank] + strides[world_rank] - 1]/
                                                          d[offset[world_rank] + strides[world_rank] - 2];
    }
    

    if (world_rank > 0){
        d[offset[world_rank] - 1] = 0;
    }
    double *tmp = new double[N];
    for (int i=0; i<N; i++){
        tmp[i] = d[i];
    }
    MPI_Allreduce(tmp, d, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    for (int i=0; i<N; i++){
        tmp[i] = l[i];
    }
    MPI_Allreduce(tmp, l, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    
    if (world_rank == 0){
        y[0] = q[0];
        for (int i=1; i<N; i++){
            y[i] = q[i] - l[i]*y[i-1];
            //std::cout << "y: " << i << y[i] << "\n";
        }

        x[N-1] = y[N-1]/d[N-1];

        for (int i=N-2; i>=0; i--){
            x[i] = (y[i] - c[i]*x[i+1])/d[i];
        }
    }

    return;
}