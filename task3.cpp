#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <functional>
#include <vector>
#include <numeric>

void magic(int& val)
{
    val += 10;
}

int main(int argc, char *argv[]) {
    int rank, count;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &count);
    
    if (rank == 0)
    {
        MPI_Status status;

        int kek = -90;
        std::cout << "Rank: " << rank << " | Init Value: " << kek << std::endl;

        MPI_Send(&kek, 1, MPI_INT, rank+1, rank, MPI_COMM_WORLD);
        MPI_Recv(&kek, 1, MPI_INT, count-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        
        std::cout << "Rank: " << rank << " | Final Value: " << kek << std::endl;
    }
    else
    {
        MPI_Status status;

        int kekR;
        
        MPI_Recv( &kekR, 1, MPI_INT, rank-1, rank-1, MPI_COMM_WORLD, &status);

        
        magic(kekR);
        std::cout << "Rank: " << rank << " | Value: " << kekR << std::endl;

        int next_rank = (rank == count-1) ? 0 : rank + 1;
        MPI_Send(&kekR, 1, MPI_INT, (rank + 1) % count, rank, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}