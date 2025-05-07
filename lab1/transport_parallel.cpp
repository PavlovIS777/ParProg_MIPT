#include <iostream>
#include <string>
#include <cmath>
#include <chrono>
#include <vector>
#include <mpi.h>
#include <cstdio>
#include <fstream>
#include <numbers>


void saveToFile(const std::vector<std::vector<double>>& data, const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    for (const auto& row : data) {
        
        for (size_t i = 0; i < row.size(); ++i) {
            out << row[i];
            if (i + 1 < row.size()) out << ' ';
        }
        out << '\n';
    }
    out.close();
}

double phi(double x) { return std::cos(M_PI*x); }
double psi(double t) { return std::exp(-t); }
double f(double x, double t) { return x+t; }

int main(int argc, char *argv[]) {
    if (argc < 6) return 0;
    auto start = std::chrono::high_resolution_clock::now();
    int K = std::stoi(argv[1]);
    int M = std::stoi(argv[2]);
    double X = std::stod(argv[3]);
    double T = std::stod(argv[4]);
    double a = std::stod(argv[5]);

    bool write = false;
    if (argc >= 7)
        write = std::stoi(argv[6]);

    double h = X / (M - 1);
    double tau = T / (K - 1);
    double fraction = tau / h;

    int rank, count;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &count);

    int sizeOffset = M / count;
    int sampleSize = sizeOffset + ((rank + 1 == count) ? M % count : 0);
    int dataSize = sampleSize + (rank != 0 ? 1 : 0);

    int globalOffset = rank * sizeOffset;

    std::vector<std::vector<double>> result;
    std::vector<double> u_prev(dataSize), u_curr(dataSize);

    if (rank == 0)
        u_prev[0] = psi(0);
    else
        MPI_Recv(&u_prev[0], 1, MPI_DOUBLE, rank-1, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int m = 1; m < dataSize; ++m) {
        u_prev[m] = phi((globalOffset + m - 1) * h);
    }
    if (rank + 1 != count)
    {
        MPI_Send(&u_prev[dataSize-1], 1, MPI_DOUBLE, rank+1, rank+1, MPI_COMM_WORLD);
    }

    result.push_back(std::vector<double>(u_prev.begin() + (rank == 0 ? 0 : 1), u_prev.end()));

    for (int k = 1; k <= K; k++) {
        if (rank == 0)
            u_curr[0] = psi(k * tau);
        else
        {
            MPI_Recv(&u_curr[0], 1, MPI_DOUBLE, rank-1, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        for (int m = 1; m < dataSize; ++m) {
            u_curr[m] = u_prev[m] - a * fraction * (u_prev[m]-u_prev[m-1]) +  tau * f((globalOffset+m-1)*h, k*tau);
        }
        if (rank + 1 != count)
        {
            MPI_Send(&u_curr[dataSize-1], 1, MPI_DOUBLE, rank+1, rank+1, MPI_COMM_WORLD);
        }
        result.push_back(std::vector<double>(u_curr.begin() + (rank == 0 ? 0 : 1), u_curr.end()));
        u_prev = u_curr;
    }

    if (rank == 0) {
        std::vector<std::vector<double>> finalData(K + 1, std::vector<double>(M));
        std::vector<MPI_Request> reqs;

        for (int i = 1; i < count; i++) {
            int chunk = sizeOffset + (i+1 == count? M % count : 0);
            int offset = i * sizeOffset;
            for (int j = 0; j < K + 1; j++) {
                reqs.push_back({});
                MPI_Irecv(finalData[j].data() + offset, chunk, MPI_DOUBLE, i, i*(K+1)+j, MPI_COMM_WORLD, &reqs.back());
            }
        }

        for (int j = 0; j < K + 1; j++) {
            std::copy(result[j].begin(), result[j].end(), finalData[j].begin());
        }

        MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "Elapsed time: " << duration.count() << " seconds" << std::endl;

        if (write)
            saveToFile(finalData, "out" + std::to_string(K) + std::to_string(M) + "_p.txt");
    } else {
        std::vector<MPI_Request> reqs(K + 1);

        for (int i = 0; i < K + 1; i++) {
            MPI_Isend(result[i].data(), dataSize-1, MPI_DOUBLE, 0, rank*(K+1)+i, MPI_COMM_WORLD, &reqs[i]);
        }

        MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
