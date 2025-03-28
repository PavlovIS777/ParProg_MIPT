#include <mpi.h>
#include <stdio.h>
#include <algorithm>
#include <functional>
#include <vector>
#include <numeric>
#include <iostream>
#include <gmpxx.h>
#include <gmp.h>
#include <cmath>

void compute_partial_sum(mpf_t sum, int start, int end)
{
    mpf_t term, factorial;
    mpf_init(term);
    mpf_init(factorial);
    mpf_set_ui(factorial, 1);

    for (int k = 1; k <= start; ++k)
    {
        
        mpf_mul_ui(factorial, factorial, k);
    }

    mpf_ui_div(term, 1, factorial);
    mpf_set(sum, term);

    for (int k = start + 1; k <= end; ++k)
    {
        mpf_div_ui(term, term, k);
        mpf_add(sum, sum, term);
    }

    mpf_clear(term);
    mpf_clear(factorial);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int precision = (argc > 1) ? std::stoi(argv[1]) : 20;
    mpf_set_default_prec(1000);
    int terms = 1;

    double kek = 0;
    while (kek < precision)
    {
        kek += log10(terms);
        terms *= 2;
    }

    int start = rank * terms / size;
    int end = rank + 1 == size ? terms : start + terms / size;

    if (rank == 0)
    {
        std::vector<mpf_t> received_sums(size-1);
        std::for_each(received_sums.begin(), received_sums.end(), [precision](auto& data)
        {
            mpf_init(data);
        });
        std::vector<MPI_Request> recv_requests(4*(size - 1));
        for (int i = 0; i < size - 1; ++i)
        {
            int sender_rank = i + 1;

            MPI_Irecv(
                &received_sums[i]->_mp_exp,
                sizeof(mp_exp_t) / sizeof(char),
                MPI_CHAR,
                sender_rank,
                4*sender_rank,
                MPI_COMM_WORLD,
                &recv_requests[i*4]
            );

            MPI_Irecv(
                &received_sums[i]->_mp_size,
                sizeof(int) / sizeof(char),
                MPI_CHAR,
                sender_rank,
                4*sender_rank + 1,
                MPI_COMM_WORLD,
                &recv_requests[i*4+1]
            );

            MPI_Irecv(
                &received_sums[i]->_mp_prec,
                sizeof(int) / sizeof(char),
                MPI_CHAR,
                sender_rank,
                4*sender_rank + 2,
                MPI_COMM_WORLD,
                &recv_requests[i*4+2]
            );
            
            MPI_Status status;
            MPI_Wait(&recv_requests[i*4+1], &status);

            received_sums[i]->_mp_d = (mp_limb_t*) calloc(received_sums[i]->_mp_size, sizeof(mp_limb_t));

            MPI_Irecv(
                received_sums[i]->_mp_d,
                sizeof(mp_limb_t)*received_sums[i]->_mp_size / sizeof(char),
                MPI_CHAR,
                sender_rank,
                4*sender_rank + 3,
                MPI_COMM_WORLD,
                &recv_requests[i*4+3]
            );
        }

        mpf_t sum;
        mpf_init(sum);
        compute_partial_sum(sum, start, end);

        MPI_Waitall(
            4*(size - 1),
            recv_requests.data(),
            MPI_STATUSES_IGNORE
        );

        for (mpf_t &pSum : received_sums)
        {
            mpf_add(sum, sum, pSum);
        }
        gmp_printf("%.*Ff\n", precision, sum);
        mpf_clear(sum);
    }
    else
    {
        mpf_t part_sum;
        mpf_init(part_sum);
        compute_partial_sum(part_sum, start, end);
        std::vector<MPI_Request> requests(4);
        std::vector<MPI_Status> statuses(4);

        MPI_Isend(
            &part_sum->_mp_exp,
            sizeof(mp_exp_t) / sizeof(char),
            MPI_CHAR,
            0,
            rank * 4,
            MPI_COMM_WORLD,
            &requests[0]);
        
        int abs_sz = abs(part_sum->_mp_size);
        MPI_Isend(
            &abs_sz,
            sizeof(int) / sizeof(char),
            MPI_CHAR,
            0,
            rank * 4 + 1,
            MPI_COMM_WORLD,
            &requests[1]);

        MPI_Isend(
            &part_sum->_mp_prec,
            sizeof(int) / sizeof(char),
            MPI_CHAR,
            0,
            rank * 4 + 2,
            MPI_COMM_WORLD,
            &requests[2]);
        
        MPI_Waitall( 3, requests.data(), statuses.data());

        MPI_Isend(
            part_sum->_mp_d,
            sizeof(mp_limb_t)*abs(part_sum->_mp_size) / sizeof(char),
            MPI_CHAR,
            0,
            rank * 4 + 3,
            MPI_COMM_WORLD,
            &requests[3]);

        MPI_Waitall( 4, requests.data(), statuses.data());
        mpf_clear(part_sum);
    }

    MPI_Finalize();
    return 0;
}