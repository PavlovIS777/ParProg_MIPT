#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <functional>
#include <vector>
#include <numeric>


int main(int argc, char *argv[]){    
    if (argc < 2)
        return -1;

    int N = -1;
    sscanf(argv[1], "%d", &N); 
    if (N <= 0)
        return -1;

    int rank, count;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &count);
    
    auto sumF{ [](int start, int end) -> double
    {
        double sum_part = 0;
        for (int i = start; i < end; i++)
        {
            sum_part += 1.0/i;
        }

        return sum_part;
    }};
    
    int start = 1+rank*N/count;
    int end = rank + 1 == count ? N+1 : start + N/count;

    if (rank == 0)
    {
        std::vector<double> received_sums(count - 1);
        std::vector<MPI_Request> recv_requests(count - 1);

        for (int i = 0; i < count - 1; ++i) {
            int sender_rank = i + 1;
            MPI_Irecv(
                &received_sums[i],   
                1,
                MPI_DOUBLE,
                sender_rank,
                sender_rank,
                MPI_COMM_WORLD,
                &recv_requests[i]
            );
        }

        auto sum_part{ sumF(start, end) };

        MPI_Waitall(
            count - 1,
            recv_requests.data(),
            MPI_STATUSES_IGNORE
        );

        auto sum{ std::accumulate(received_sums.begin(), received_sums.end(), sum_part) };
        std::cout << sum << std::endl;
    }
    else
    {
        auto sum_part{ sumF(start, end) };
        MPI_Request request;
        MPI_Isend(&sum_part, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, &request);
    }
    MPI_Finalize();
}