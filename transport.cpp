#include <iostream>
#include <string>
#include <cmath>
#include <chrono>
#include <vector>
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

int main(int argc, char** argv) {
    if (argc < 6) return 0;

    auto start = std::chrono::high_resolution_clock::now();

    int K = std::stoi(argv[1]);
    int M = std::stoi(argv[2]);
    double X = std::stod(argv[3]);
    double T = std::stod(argv[4]);
    double a = std::stod(argv[5]);

    double h = X / (M - 1);
    double tau = T / (K - 1);

    int dataSize = M + 2; // 2 призрачных слоя

    std::vector<std::vector<double>> result;
    std::vector<double> u_prev(dataSize), u_curr(dataSize), u_next(dataSize);

    // начальное условие
    for (int m = 0; m < M; ++m) {
        u_prev[m + 1] = phi(m * h);
    }

    // граничные условия на первом слое
    u_prev[0] = psi(0);
    u_prev[M + 1] = 3 * u_prev[M] - 3 * u_prev[M - 1] + u_prev[M - 2];

    // граничные условия на втором слое
    u_curr[0] = psi(tau);

    // шаг 1: t = tau
    for (int m = 1; m <= M; ++m) {
        u_curr[m] = u_prev[m] - a * (tau / h) * (u_prev[m + 1] - u_prev[m - 1]) + 2 * tau * f((m - 1) * h, 0);
    }
    u_curr[M + 1] = 3 * u_curr[M] - 3 * u_curr[M - 1] + u_curr[M - 2];

    result.push_back(std::vector<double>(u_prev.begin() + 1, u_prev.end() - 1));
    result.push_back(std::vector<double>(u_curr.begin() + 1, u_curr.end() - 1));

    // остальные временные слои
    for (int k = 2; k <= K; ++k) {
        u_next[0] = psi(k * tau);
        for (int m = 1; m <= M; ++m) {
            u_next[m] = u_prev[m] - a * (tau / h) * (u_curr[m + 1] - u_curr[m - 1]) + 2 * tau * f((m - 1) * h, k * tau);
        }
        u_next[M + 1] = 3 * u_next[M] - 3 * u_next[M - 1] + u_next[M - 2];

        result.push_back(std::vector<double>(u_next.begin() + 1, u_next.end() - 1));

        u_prev = u_curr;
        u_curr = u_next;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Elapsed time: " << duration.count() << " seconds" << std::endl;

    saveToFile(result, "out.txt");

    return 0;
}
