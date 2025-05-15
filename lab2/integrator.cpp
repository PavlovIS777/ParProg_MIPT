#include <iostream>
#include <pthread.h>
#include <mpfr.h>
#include <vector>
#include <cmath>

struct ThreadData {
    mpfr_t a, h, sum, eps;
    long start, end;
    mpfr_prec_t prec;
    mpfr_rnd_t rnd;
};

void* integrate_thread(void* arg) {
    ThreadData* data = (ThreadData*)arg;

    mpfr_t xi, fx, temp;
    mpfr_inits2(data->prec, xi, fx, temp, nullptr);
    
    for (long i = data->start; i < data->end; ++i) {
        mpfr_set_si(temp, 2 * i + 1, data->rnd);       // temp = 2i + 1
        mpfr_mul(xi, data->h, temp, data->rnd);        // xi = h * (2i + 1)
        mpfr_div_si(xi, xi, 2, data->rnd);             // xi = h * (2i + 1) / 2
        mpfr_add(xi, data->a, xi, data->rnd);          // xi = a + h*(2i+1)/2

        mpfr_cos(fx, xi, data->rnd);                   // fx = cos(xi)
        mpfr_mul(temp, xi, xi, data->rnd);             // temp = xi^2
        mpfr_div(fx, fx, temp, data->rnd);             // fx = cos(xi)/xi^2

        mpfr_add(data->sum, data->sum, fx, data->rnd); // sum += fx
    }

    mpfr_clears(xi, fx, temp, nullptr);
    return nullptr;
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <num_threads> <eps> <a> <b>\n";
        return 1;
    }

    int num_threads = std::stoi(argv[1]);
    double eps_d = std::stod(argv[2]);
    double a_d = std::stod(argv[3]);
    double b_d = std::stod(argv[4]);

    if (a_d <= 0.0) {
        std::cerr << "Error: a must be > 0 to avoid division by zero\n";
        return 1;
    }

    mpfr_prec_t prec = 256;
    mpfr_rnd_t rnd = MPFR_RNDN;

    mpfr_t a, b, eps, prev_result, result, h, diff, temp;
    mpfr_inits2(prec, a, b, eps, prev_result, result, h, diff, temp, nullptr);

    mpfr_set_d(a, a_d, rnd);
    mpfr_set_d(b, b_d, rnd);
    mpfr_set_d(eps, eps_d, rnd);

    long n = 1000;
    bool converged = false;

    while (!converged) {
        mpfr_set_ui(result, 0, rnd);
        mpfr_sub(temp, b, a, rnd); // temp = b - a
        mpfr_div_si(h, temp, n, rnd); // h = (b - a) / n

        std::vector<pthread_t> threads(num_threads);
        std::vector<ThreadData> data(num_threads);

        long chunk = n / num_threads;

        for (int i = 0; i < num_threads; ++i) {
            data[i].start = i * chunk;
            data[i].end = (i == num_threads - 1) ? n : (i + 1) * chunk;
            data[i].prec = prec;
            data[i].rnd = rnd;

            mpfr_inits2(prec, data[i].a, data[i].h, data[i].sum, data[i].eps, nullptr);
            mpfr_set(data[i].a, a, rnd);
            mpfr_set(data[i].h, h, rnd);
            mpfr_set_ui(data[i].sum, 0, rnd);
            mpfr_set(data[i].eps, eps, rnd);

            pthread_create(&threads[i], nullptr, integrate_thread, &data[i]);
        }

        for (int i = 0; i < num_threads; ++i) {
            pthread_join(threads[i], nullptr);
            mpfr_add(result, result, data[i].sum, rnd);
            mpfr_clears(data[i].a, data[i].h, data[i].sum, data[i].eps, nullptr);
        }

        mpfr_mul(result, result, h, rnd); // result *= h

        mpfr_sub(diff, result, prev_result, rnd);
        mpfr_abs(diff, diff, rnd);
        if (mpfr_cmp(diff, eps) < 0)
            converged = true;
        else {
            mpfr_set(prev_result, result, rnd);
            n *= 2;
        }
    }

    std::cout << "Integral of cos(x)/x^2 from " << a_d << " to " << b_d << " ≈ ";
    mpfr_out_str(stdout, 10, 50, result, rnd);
    std::cout << std::endl;

    mpfr_clears(a, b, eps, prev_result, result, h, diff, temp, nullptr);
    return 0;
}
