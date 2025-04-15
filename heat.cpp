#include <iostream>
#include <string>
#include <cmath>
#include <chrono>


double phi(x):
    return x;

double psi(t):
    return t;

double f(x, t):
    return x*t;

int main(int argc, char** argv) {    
    if (argc < 2) {
        return 0;
    }

    int K = std::stoi(argv[1]);
    int M = std::stoi(argv[2]);
    int X = std::stoi(argv[3]);
    int T = std::stoi(argv[4]);

    double h = X/M;
    double tau = T/K;

    std::vector<std::vector<double>> result(K+1);
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

    for (int m = 1; m < M+2; m++)
    {
        u_curr = u_prev[m] - tau * (u_prev[m+1] + u_prev[m-1]) / (2*h) + tau * f(h*m, 0) + tau*2*(u_prev[m+1]-2*u_prev[m] + u_prev[m-1])/(2*h**2)
    }
    u_curr[M+1] = 3*u_curr[M] - 3*u_curr[M-1] + u_curr[M-2];

    result.push_back(std::vector<double>({u_curr.begin(), u_curr.end()-1}));

    for (int k = 2; k < K+1; k++)
    {
        u_next[0] = psi(k*tau);
        for (int m = 1; m < M+2; m++)
        {
            u_next[m] = u_prev[m] - 2*tau*((u_curr[m+1]-u[m-1])/(2*h) - f(m*h, k*tau));
        }
        u_next[M+1] = 3*u_next[M] - 3*u_next[M-1] + u_next[M-2];
        result.push_back(std::vector<double>({u_next.begin(), u_next.end()-1}));

        u_prev = u_curr;
        u_curr = u_next;
    }

    for (auto& tLayer : result)
    {
        for (auto& x : tLayer)
        {
            std::cout << x << " "; 
        }

        std::cout << std::endl;
    }

    return 0;
}