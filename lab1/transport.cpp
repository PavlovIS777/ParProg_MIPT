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

double phi(double x) { return std::cos(M_PI*x); }
double psi(double t) { return std::exp(-t); }
double f(double x, double t) { return x+t; }

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
    double fraction = tau/h;

    if (std::abs(a*fraction) >= 1 )
    {
        std::cout << "Схема не сходится, r = " << std::abs(a*fraction) << std::endl;
        return 0;
    }

    int dataSize = M;

    std::vector<std::vector<double>> result;
    std::vector<double> u_prev(dataSize), u_curr(dataSize);

    for (int m = 0; m < M; ++m) {
        u_prev[m] = phi(m * h);
    }
    u_prev[0] = psi(0);

    result.push_back(std::vector<double>(u_prev.begin(), u_prev.end()));

    // остальные временные слои
    for (int k = 1; k <= K; ++k) {
        u_curr[0] = psi(k * tau);
        for (int m = 1; m < M; ++m) {
            u_curr[m] = u_prev[m] - a*fraction * (u_prev[m]-u_prev[m-1]) +  tau * f(m*h, k*tau);
        }
        result.push_back(std::vector<double>(u_curr.begin(), u_curr.end()));
        u_prev = u_curr;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Elapsed time: " << duration.count() << " seconds" << std::endl;

    saveToFile(result, "out" + std::to_string(K) + std::to_string(M) + ".txt");

    return 0;
}
