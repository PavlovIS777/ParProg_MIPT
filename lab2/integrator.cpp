#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <pthread.h>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <mpfr.h>


struct ThreadData {
    mpf_class start;
    mpf_class end;
    mpf_class precision;
    mpf_class result;
};

void f_mpfr(mpf_class& result, const mpf_class& x, mpfr_rnd_t rnd = MPFR_RNDN) {
    if (sgn(x) == 0) {
        result = 0;
        return;
    }

    mpfr_t x_mpfr, cos_x, x_squared;
    mpfr_inits2(256, x_mpfr, cos_x, x_squared, (mpfr_ptr) 0);

    mpfr_set_f(x_mpfr, x.get_mpf_t(), rnd);         // x → mpfr
    mpfr_cos(cos_x, x_mpfr, rnd);                   // cos(x)
    mpfr_mul(x_squared, x_mpfr, x_mpfr, rnd);       // x^2
    mpfr_div(cos_x, cos_x, x_squared, rnd);         // cos(x)/x^2

    mpf_set(result.get_mpf_t(), cos_x);             // convert back to mpf_class

    mpfr_clears(x_mpfr, cos_x, x_squared, (mpfr_ptr) 0);
}


void compute_integral_segment(ThreadData* data) {
    mpf_class x = data->start;
    mpf_class step = data->precision;
    mpf_class sum = 0;

    while ( x < data->end) {
        mpf_class fx;
        f_mpfr(fx, x);
        sum += fx * step;
    }

    data->result = sum;
}

void* thread_function(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    compute_integral_segment(data);
    return nullptr;
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <num_threads> <precision> <start:end>\n";
        return 1;
    }

    int num_threads = std::stoi(argv[1]);
    mpf_class precision(argv[2]);

    char* range_str = argv[3];
    char* colon = strchr(range_str, ':');
    if (!colon) {
        std::cerr << "Invalid range format. Use start:end\n";
        return 1;
    }

    *colon = '\0';
    mpf_class range_start(range_str);
    mpf_class range_end(colon + 1);

    // Calculate interval length
    mpf_class total_length = range_end - range_start;
    mpf_class thread_length = total_length / num_threads;

    pthread_t* threads = new pthread_t[num_threads];
    ThreadData* thread_data = new ThreadData[num_threads];

    clock_t start_time = clock();

    for (int i = 0; i < num_threads; ++i) {
        thread_data[i].start = range_start + i * thread_length;
        thread_data[i].end = (i == num_threads - 1) ? range_end : (range_start + (i + 1) * thread_length);
        thread_data[i].precision = precision / num_threads;
        pthread_create(&threads[i], nullptr, thread_function, &thread_data[i]);
    }

    mpf_class total_result = 0;

    for (int i = 0; i < num_threads; ++i) {
        pthread_join(threads[i], nullptr);
        total_result += thread_data[i].result;
    }

    clock_t end_time = clock();
    double elapsed_secs = double(end_time - start_time) / CLOCKS_PER_SEC;

    std::cout.precision(20);
    std::cout << "Result: " << total_result << "\n";
    std::cout << "Time: " << elapsed_secs << " s\n";

    delete[] threads;
    delete[] thread_data;

    return 0;
}
