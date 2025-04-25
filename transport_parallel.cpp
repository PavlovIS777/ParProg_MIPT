#include <iostream>
#include <string>
#include <cmath>
#include <chrono>
#include <vector>
#include <mpi.h>
#include <cstdio>
#include <fstream>


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

double phi(double x) { return std::sin(10 * M_PI * x); }
double psi(double t) { return 0; }
double f(double x, double t) { return 0; }

int main(int argc, char *argv[]) {
    if (argc < 6) return 0;
    auto start = std::chrono::high_resolution_clock::now();
    int K = std::stoi(argv[1]);
    int M = std::stoi(argv[2]);
    double X = std::stod(argv[3]);
    double T = std::stod(argv[4]);
    double a = std::stod(argv[5]);

    double h = X / (M - 1);
    double tau = T / (K - 1);

    int rank, count;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &count);

    int sizeOffset = M / count;
    int sampleSize = sizeOffset + ((rank + 1 == count) ? M % count : 0);

    int globalOffset = rank * sizeOffset;
    int dataSize = sampleSize + (rank == 0 ? 1 : 2);

    std::vector<std::vector<double>> result;
    std::vector<double> u_prev(dataSize), u_curr(dataSize), u_next(dataSize);

    for (int m = 1; m < dataSize-1; m++) {
        u_prev[m] = phi((globalOffset + m) * h);
    }
    u_prev[0] = phi((globalOffset - 1) * h);
    u_prev[dataSize-1] = phi((globalOffset + sampleSize + 1) * h);

    if (rank == 0) {
        u_prev[0] = psi(0);
        u_curr[0] = psi(tau);
    }

    if (rank + 1 == count) {
        u_prev[dataSize - 1] = 3 * u_prev[dataSize - 2] - 3 * u_prev[dataSize - 3] + u_prev[dataSize - 4];
    }

    if (rank != 0) {
        MPI_Recv(&u_prev[0], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&u_curr[0], 1, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    for (int m = 1; m < dataSize-1; m++) {
        u_curr[m] = u_prev[m] - a * (tau / h) * (u_prev[m + 1] - u_prev[m - 1]) + 2 * tau * f((globalOffset + m - 1) * h, 0);
    }
    u_curr[dataSize-1] = 3 * u_curr[dataSize-2] - 3 * u_curr[dataSize-3] + u_curr[dataSize-4];

    if (rank + 1 != count) {
        MPI_Send(&u_prev[dataSize - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        MPI_Send(&u_curr[dataSize - 1], 1, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
    }

    result.push_back(std::vector<double>(u_prev.begin() + (rank == 0 ? 0 : 1), u_prev.end() - 1));
    result.push_back(std::vector<double>(u_curr.begin() + (rank == 0 ? 0 : 1), u_curr.end() - 1));

    for (int k = 2; k < K + 1; k++) {
        std::vector<MPI_Request> reqs;

        if (rank == 0)
            u_next[0] = psi(k * tau);
        if (rank + 1 == count)
            u_next[dataSize-1] = 3*u_next[dataSize-2] - 3*u_next[dataSize-3] + u_next[dataSize-4];
        
        if (rank != 0)
        {
            reqs.emplace_back();
            MPI_Irecv(&u_next[0], 1, MPI_DOUBLE, rank-1, 2, MPI_COMM_WORLD, &reqs.back());
        }
        if (rank + 1 != count)
        {
            reqs.emplace_back();
            MPI_Irecv(&u_next[dataSize-1], 1, MPI_DOUBLE, rank+1, 3, MPI_COMM_WORLD, &reqs.back());
        }
    
        for (int m = 1; m < dataSize-1; m++) {
            u_next[m] = u_prev[m] - a * (tau / h) * (u_curr[m + 1] - u_curr[m - 1]) + 2 * tau * f((globalOffset + m - 1) * h, k * tau);
        }

        if (rank + 1 != count )
        {
            reqs.emplace_back();
            MPI_Isend(&u_next[dataSize-2], 1, MPI_DOUBLE, rank+1, 2, MPI_COMM_WORLD, &reqs.back());
        }

        if (rank != 0)
        {
            reqs.emplace_back();
            MPI_Isend(&u_next[1], 1, MPI_DOUBLE, rank-1, 3, MPI_COMM_WORLD, &reqs.back());
        }

        MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUS_IGNORE);
        result.push_back(std::vector<double>(u_next.begin() + (rank == 0 ? 0 : 1), u_next.end() - 1));
        u_prev = u_curr;
        u_curr = u_next;
    }

    if (rank == 0) {
        std::vector<std::vector<double>> finalData(K + 1, std::vector<double>(M));
        std::vector<MPI_Request> reqs((count - 1) * (K + 1));

        for (int i = 1; i < count; i++) {
            int chunk = sizeOffset + (i+1 == count? M % count : 0);
            int offset = i * sizeOffset;
            for (int j = 0; j < K + 1; j++) {
                MPI_Irecv(finalData[j].data() + offset, chunk, MPI_DOUBLE, i, i*(K+1)+j, MPI_COMM_WORLD, &reqs[(i - 1) * (K + 1) + j]);
            }
        }

        for (int j = 0; j < K + 1; j++) {
            std::copy(result[j].begin(), result[j].end(), finalData[j].begin());
        }

        MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "Elapsed time: " << duration.count() << " seconds" << std::endl;

        saveToFile(finalData, "outParallel.txt");

    } else {
        std::vector<MPI_Request> reqs(K + 1);

        for (int i = 0; i < K + 1; i++) {
            MPI_Isend(result[i].data(), result[i].size(), MPI_DOUBLE, 0, rank*(K+1)+i, MPI_COMM_WORLD, &reqs[i]);
        }
        MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
