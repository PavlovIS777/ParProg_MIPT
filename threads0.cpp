#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>
#include <unistd.h>

#define MAX_THREADS 128

pthread_mutex_t mutex;
int shared_var = 0;

typedef struct {
    int thread_num;
    int total_threads;
} HelloArgs;

void* hello_thread(void* arg) {
    HelloArgs* args = (HelloArgs*) arg;
    printf("Hello World from thread %d of %d\n", args->thread_num, args->total_threads);
    return NULL;
}

void run_hello_world(int thread_count) {
    pthread_t threads[thread_count];
    HelloArgs args[thread_count];

    for (int i = 0; i < thread_count; ++i) {
        args[i].thread_num = i;
        args[i].total_threads = thread_count;
        pthread_create(&threads[i], NULL, hello_thread, &args[i]);
    }

    for (int i = 0; i < thread_count; ++i) {
        pthread_join(threads[i], NULL);
    }
}

typedef struct {
    int start;
    int end;
    double partial_sum;
} SumArgs;

void* sum_thread(void* arg) {
    SumArgs* args = (SumArgs*) arg;
    args->partial_sum = 0.0;
    for (int i = args->start; i <= args->end; ++i) {
        args->partial_sum += 1.0 / i;
    }
    return NULL;
}

void run_harmonic_sum(int N, int thread_count) {
    pthread_t threads[thread_count];
    SumArgs args[thread_count];

    int chunk = N / thread_count;
    int remainder = N % thread_count;
    int current = 1;

    for (int i = 0; i < thread_count; ++i) {
        args[i].start = current;
        args[i].end = current + chunk - 1 + (i < remainder ? 1 : 0);
        current = args[i].end + 1;
        pthread_create(&threads[i], NULL, sum_thread, &args[i]);
    }

    double total = 0.0;
    for (int i = 0; i < thread_count; ++i) {
        pthread_join(threads[i], NULL);
        total += args[i].partial_sum;
    }

    printf("Partial harmonic sum 1/k for k=1 to %d = %.10f\n", N, total);
}


typedef struct {
    int thread_num;
} SeqArgs;

void* seq_thread(void* arg) {
    SeqArgs* args = (SeqArgs*) arg;

    pthread_mutex_lock(&mutex);
    shared_var++;
    printf("Thread %d incremented shared_var to %d\n", args->thread_num, shared_var);
    pthread_mutex_unlock(&mutex);

    return NULL;
}

void run_sequential_access(int thread_count) {
    pthread_t threads[thread_count];
    SeqArgs args[thread_count];

    pthread_mutex_init(&mutex, NULL);
    shared_var = 0;

    for (int i = 0; i < thread_count; ++i) {
        args[i].thread_num = i;
        pthread_create(&threads[i], NULL, seq_thread, &args[i]);
    }

    for (int i = 0; i < thread_count; ++i) {
        pthread_join(threads[i], NULL);
    }

    pthread_mutex_destroy(&mutex);
}


int main(int argc, char* argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s [task_number 1|2|3] [thread_count] [N - only for task 2]\n", argv[0]);
        return 1;
    }

    int task = atoi(argv[1]);
    int thread_count = atoi(argv[2]);

    if (thread_count < 1 || thread_count > MAX_THREADS) {
        fprintf(stderr, "Invalid thread count\n");
        return 1;
    }

    switch (task) {
        case 1:
            run_hello_world(thread_count);
            break;
        case 2:
            if (argc < 4) {
                fprintf(stderr, "For task 2, provide N (number of terms)\n");
                return 1;
            }
            run_harmonic_sum(atoi(argv[3]), thread_count);
            break;
        case 3:
            run_sequential_access(thread_count);
            break;
        default:
            fprintf(stderr, "Invalid task number. Choose 1, 2 or 3.\n");
            return 1;
    }

    return 0;
}
