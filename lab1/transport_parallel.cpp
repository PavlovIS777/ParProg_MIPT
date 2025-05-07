#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <mpi.h>
#include <chrono>

double phi(double x) { return std::cos(M_PI * x); }
double psi(double t) { return std::exp(-t); }
double f(double x, double t) { return x + t; }

void saveToFile(const std::vector<std::vector<double>> &data, const std::string &filename) {
    std::ofstream out(filename);
    for (const auto &row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            out << row[i] << (i + 1 < row.size() ? ' ' : '\n');
        }
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 7) {
        if (rank == 0)
            std::cerr << "Usage: ./a.out <K> <M> <X> <T> <a> <write=0|1>\n";
        MPI_Finalize();
        return 1;
    }

    int K = std::stoi(argv[2]);
    int M = std::stoi(argv[3]);
    double X = std::stod(argv[4]);
    double T = std::stod(argv[5]);
    double a = std::stod(argv[6]);
    bool write = argc >= 8 ? std::stoi(argv[7]) : false;

    double h = X / (M - 1);
    double tau = T / (K - 1);
    double frac = tau / h;
    auto start_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration;
    std::chrono::duration<double> calc_duration;

    int base = M / size;
    int extra = M % size;
    int local_n = base + (rank < extra ? 1 : 0);
    int offset = rank * base + std::min(rank, extra);

    std::vector<double> u_prev(local_n + 1);
    std::vector<double> u_curr(local_n + 1);
    std::vector<std::vector<double>> local_result;

    if (rank == 0)
        u_prev[0] = psi(0);

    for (int i = 1; i <= local_n; ++i)
        u_prev[i] = phi((offset + i - 1) * h);

    local_result.push_back(std::vector<double>(u_prev.begin() + 1, u_prev.end()));

    MPI_Status status;
    for (int k = 1; k <= K; ++k) {
        if (rank != 0)
            MPI_Recv(&u_prev[0], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
        if (rank != size - 1)
            MPI_Send(&u_prev[local_n], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);

        if (rank == 0)
            u_curr[1] = psi(k * tau);

        for (int i = (rank == 0 ? 2 : 1); i <= local_n; ++i) {
            double x = (offset + i - 1) * h;
            double t = k * tau;
            u_curr[i] = u_prev[i] - a * frac * (u_prev[i] - u_prev[i - 1]) + tau * f(x, t);
        }

        u_prev.swap(u_curr);
        local_result.push_back(std::vector<double>(u_prev.begin() + 1, u_prev.end()));
    }

    auto end_calc_time = std::chrono::high_resolution_clock::now();
    calc_duration = end_calc_time - start_time;

    // Сбор данных
    {
        std::vector<double> sendbuf;
        for (int k = 0; k <= K; ++k)
            sendbuf.insert(sendbuf.end(), local_result[k].begin(), local_result[k].end());

        int local_size = local_n * (K + 1);

        std::vector<double> recvbuf;
        std::vector<int> recvcounts, displs;

        if (rank == 0) {
            recvbuf.resize((M) * (K + 1));

            recvcounts.resize(size);
            displs.resize(size);

            for (int r = 0; r < size; ++r) {
                int local_r = base + (r < extra ? 1 : 0);
                recvcounts[r] = local_r * (K + 1);
                displs[r] = (r == 0 ? 0 : displs[r - 1] + recvcounts[r - 1]);
            }
        }

        MPI_Gatherv(sendbuf.data(), local_size, MPI_DOUBLE,
                    recvbuf.data(), recvcounts.data(), displs.data(), MPI_DOUBLE,
                    0, MPI_COMM_WORLD);

        if (rank == 0) {
            std::vector<std::vector<double>> full(K + 1, std::vector<double>(M));

            for (int r = 0, pos = 0; r < size; ++r) {
                int local_r = base + (r < extra ? 1 : 0);
                for (int k = 0; k <= K; ++k) {
                    for (int j = 0; j < local_r; ++j) {
                        full[k][r * base + std::min(r, extra) + j] = recvbuf[pos++];
                    }
                }
            }

            auto end_time = std::chrono::high_resolution_clock::now();
            duration = end_time - start_time;

            if (write)
                saveToFile(full, "out" + std::to_string(K) + std::to_string(M) + "_p.txt");
        }

        // std::cout << "Elapsed calculation time for rank: " << rank+1 << " - " << calc_duration.count() << " seconds\n";
        // MPI_Barrier(MPI_COMM_WORLD);
        if(rank == 0)
        {
            std::cout << "Elapsed time: " << duration.count() << " seconds\n";
        }
    }
    MPI_Finalize();
    return 0;
}
