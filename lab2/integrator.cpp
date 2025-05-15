#include <iostream>
#include <pthread.h>
#include <mpfr.h>
#include <vector>

struct ThreadData {
    mpfr_t a, h, sum;
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
        mpfr_div_si(xi, xi, 2, data->rnd);             // xi /= 2
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

    mpfr_t a, b, eps, h, temp, diff, prev, curr;
    mpfr_inits2(prec, a, b, eps, h, temp, diff, prev, curr, nullptr);
    mpfr_set_d(a, a_d, rnd);
    mpfr_set_d(b, b_d, rnd);
    mpfr_set_d(eps, eps_d, rnd);

    long n = 4;
    bool converged = false;
    long m = 10;

    // Этап 1: Найти подходящее n
    while (!converged) {
        mpfr_sub(temp, b, a, rnd);          // temp = b - a
        mpfr_div_si(h, temp, n, rnd);       // h = (b - a) / n

        mpfr_set_ui(prev, 0, rnd);
        mpfr_set_ui(curr, 0, rnd);

        // Интеграл на [a, a+h] с m отрезками
        mpfr_t h_m;
        mpfr_init2(h_m, prec);
        mpfr_div_si(h_m, h, m, rnd);

        // Первый прогон
        mpfr_set_ui(prev, 0, rnd);
        for (long i = 0; i < m; ++i) {
            mpfr_set_si(temp, 2 * i + 1, rnd);
            mpfr_mul(temp, temp, h_m, rnd);
            mpfr_div_si(temp, temp, 2, rnd);
            mpfr_add(temp, a, temp, rnd);

            mpfr_cos(curr, temp, rnd);
            mpfr_mul(temp, temp, temp, rnd);
            mpfr_div(curr, curr, temp, rnd);
            mpfr_add(prev, prev, curr, rnd);
        }
        mpfr_mul(prev, prev, h_m, rnd);

        // Второй прогон с удвоением m
        mpfr_set_ui(curr, 0, rnd);
        long m2 = 2 * m;
        mpfr_div_si(h_m, h, m2, rnd);

        for (long i = 0; i < m2; ++i) {
            mpfr_set_si(temp, 2 * i + 1, rnd);
            mpfr_mul(temp, temp, h_m, rnd);
            mpfr_div_si(temp, temp, 2, rnd);
            mpfr_add(temp, a, temp, rnd);

            mpfr_cos(diff, temp, rnd);
            mpfr_mul(temp, temp, temp, rnd);
            mpfr_div(diff, diff, temp, rnd);
            mpfr_add(curr, curr, diff, rnd);
        }
        mpfr_mul(curr, curr, h_m, rnd);

        mpfr_sub(diff, curr, prev, rnd);
        mpfr_abs(diff, diff, rnd);

        mpfr_div_si(temp, eps, n, rnd);
        if (mpfr_cmp(diff, temp) < 0) {
            converged = true;
        } else {
            n *= 2;
        }

        mpfr_clear(h_m);
    }

    // Этап 2: Параллельный расчет по всему отрезку с найденным n
    mpfr_sub(temp, b, a, rnd);
    mpfr_div_si(h, temp, n, rnd);

    std::vector<pthread_t> threads(num_threads);
    std::vector<ThreadData> data(num_threads);
    long chunk = n / num_threads;

    mpfr_t final_result;
    mpfr_init2(final_result, prec);
    mpfr_set_ui(final_result, 0, rnd);

    for (int i = 0; i < num_threads; ++i) {
        data[i].start = i * chunk;
        data[i].end = (i == num_threads - 1) ? n : (i + 1) * chunk;
        data[i].prec = prec;
        data[i].rnd = rnd;

        mpfr_inits2(prec, data[i].a, data[i].h, data[i].sum, nullptr);
        mpfr_set(data[i].a, a, rnd);
        mpfr_set(data[i].h, h, rnd);
        mpfr_set_ui(data[i].sum, 0, rnd);

        pthread_create(&threads[i], nullptr, integrate_thread, &data[i]);
    }

    for (int i = 0; i < num_threads; ++i) {
        pthread_join(threads[i], nullptr);
        mpfr_add(final_result, final_result, data[i].sum, rnd);
        mpfr_clears(data[i].a, data[i].h, data[i].sum, nullptr);
    }

    mpfr_mul(final_result, final_result, h, rnd);

    std::cout << "Integral of cos(x)/x^2 from " << a_d << " to " << b_d << " ≈ ";
    mpfr_out_str(stdout, 10, 50, final_result, rnd);
    std::cout << "\n";
    std::cout << "Used n = " << n << " subintervals.\n";

    mpfr_clears(a, b, eps, h, temp, diff, prev, curr, final_result, nullptr);
    return 0;
}
