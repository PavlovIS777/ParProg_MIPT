#include <mpi.h>
#include <stdio.h>
int main(int argc, char *argv[]){
    int rank, count;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &count);
    printf("Hello World, rank: %d, count: %d\n", rank, count);
    MPI_Finalize();
}