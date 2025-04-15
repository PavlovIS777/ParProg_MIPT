#include <iostream>
#include <string>
#include <cmath>
#include <chrono>
#include <vector>


double phi(double x)
{
    return 0;
}

double psi(double t)
{
    return 0;
}
double f(double x, double t)
{
    return 100;
}

int main(int argc, char** argv) {    
    if (argc < 2) {
        return 0;
    }

    int K = std::stoi(argv[1]);
    int M = std::stoi(argv[2]);
    double X = std::stod(argv[3]);
    double T = std::stod(argv[4]);
    double a = std::stod(argv[5]);

    double h = X/M;
    double tau = T/K;
    double r = a*a * tau / (h*h);

    if (r > 0.5)
    {
        std::cout << "не устойчива";
        return 0;
    }
    
    std::vector<std::vector<double>> result;
    std::vector<double> u_prev(M+2);
    std::vector<double> u_curr(M+2);
    std::vector<double> u_next(M+2);

    // Начальный конфиг
    for (int m = 0; m < M+1; m++)
    {
        u_prev[m] = phi(m*h);
    }

    u_prev[0] = psi(tau*0);
    u_curr[0] = psi(tau*1);

    u_prev[M+1] = 3*u_prev[M] - 3*u_prev[M-1] + u_prev[M-2];
    result.push_back(std::vector<double>(u_prev.begin(), u_prev.end()-1));

    for (int m = 1; m < M+2; m++)
    {
        u_curr[m] = u_prev[m] - tau * (u_prev[m+1] + u_prev[m-1])/h + 2*tau * f(h*m, 0);
    }
    u_curr[M+1] = 3*u_curr[M] - 3*u_curr[M-1] + u_curr[M-2];

    result.push_back(std::vector<double>(u_curr.begin(), u_curr.end()-1));

    for (int k = 2; k < K+1; k++)
    {
        u_next[0] = psi(k*tau);
        for (int m = 1; m < M+2; m++)
        {
            u_next[m] = u_prev[m] - tau*((u_curr[m+1]-u_curr[m-1])/h - 2*f(m*h, k*tau));
        }
        u_next[M+1] = 3*u_next[M] - 3*u_next[M-1] + u_next[M-2];
        result.push_back(std::vector<double>(u_next.begin(), u_next.end()-1));

        u_prev = u_curr;
        u_curr = u_next;
    }

    for (auto& tLayer : result)
    {
        for (int i = 0; i < tLayer.size() - 1; i++)
        {
            std::cout << tLayer[i] << " "; 
        }
        std::cout << tLayer.back() << "\n";
    }

    return 0;
}