/**
 * autor: Angelos Motsios
 * Semester: 2nd Semester 2024-2025
 * Master of Computer, Network and Telecommunication Engineering
 * Department of Informatics, University of Athens
 */

#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <ctime>
#include <iomanip>
#include <memory>
#include <random>
#include <cmath>
#include <fstream>
#include <chrono>
#include <cstdint>
#include <limits>

std::mutex output_mutex;
unsigned long long NUM_ITER = 100000000ULL;

// #define USE_INT64
// #define USE_FLOAT
// #define USE_DOUBLE
// #define USE_MATMUL

#ifdef USE_INT64
using num_t = double; // For result storage
#define TEST_NAME "int64_calc"
#elif defined(USE_FLOAT)
using num_t = float;
#define TEST_NAME "float_calc"
#elif defined(USE_DOUBLE)
using num_t = double;
#define TEST_NAME "double_calc"
#elif defined(USE_MATMUL)
using num_t = double;
#define TEST_NAME "mat_mul"
#elif defined(USE_ROOT)
using num_t = double;
#define TEST_NAME "root_calc"
#else
using num_t = float;
#define TEST_NAME "error"
#error "Please define a test (USE_FLOAT, USE_DOUBLE, USE_INT64, USE_MATMUL)"
#endif

#ifdef USE_INT64
void calc(const std::vector<int64_t>& numbers, double* result) {
    const size_t N = numbers.size();

    int64_t acc = 1;
    for (size_t i = 0; i < NUM_ITER; ++i) {
        int64_t num = numbers[i % N];
        if (acc > std::numeric_limits<int64_t>::max() / num) {
            acc %= 1000000007LL;
        }
        acc *= num;
    }

    int64_t divisor = numbers[N / 2];
    if (divisor == 0) divisor = 2;

    int64_t final_result = acc / divisor;
    if (final_result == 0) final_result = 1;

    *result = static_cast<double>(final_result);
}
#elif defined(USE_FLOAT)
void calc(const std::vector<num_t>& numbers, num_t divisor, num_t* result) {
    num_t acc = 0.0;
    size_t N = numbers.size();
    for (size_t i = 0; i < NUM_ITER; ++i) {
        if (i % 2 == 0)
            acc += numbers[i % N];
        else
            acc -= numbers[i % N];
    }
    acc /= divisor;
    acc *= divisor;
    *result = std::sqrt(std::abs(acc));
}
#elif defined(USE_DOUBLE)
void calc(const std::vector<num_t>& numbers, num_t divisor, num_t* result) {
    num_t acc = 0.0;
    size_t N = numbers.size();
    for (size_t i = 0; i < NUM_ITER; ++i) {
        if (i % 2 == 0)
            acc += numbers[i % N];
        else
            acc -= numbers[i % N];
    }
    acc /= divisor;
    acc *= divisor;
    *result = std::sqrt(std::abs(acc));
}
#elif defined(USE_MATMUL)
void mat_mul(double* result,
             const std::vector<std::vector<double>>& A,
             const std::vector<std::vector<double>>& B) {
    size_t N = A.size();
    double sum = 0.0;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            double cell = 0.0;
            for (size_t k = 0; k < N; ++k) {
                cell += A[i][k] * B[k][j];
            }
            sum += cell;
        }
    }
    *result = sum;
}
#elif defined(USE_ROOT)
void calc(const std::vector<num_t>& numbers, num_t* result) {
    num_t acc = 0.0;
    size_t N = numbers.size();
    for (size_t i = 0; i < NUM_ITER; ++i) {
        num_t val = numbers[i % N];
        if (i % 2 == 0)
            acc += std::sqrt(std::abs(val)) + std::cbrt(std::abs(val + 1.0));
        else
            acc -= std::sqrt(std::abs(val)) + std::cbrt(std::abs(val + 1.0));
        acc += std::sqrt(std::abs(acc + 2.0));
        acc -= std::cbrt(std::abs(acc + 3.0));
    }
    *result = std::sqrt(std::abs(acc));
}
#endif

int main(int argc, char* argv[]) {

    long num_cores = 1;
    time_t now = time(0);
    tm *ltm = localtime(&now);

    if (argc > 1) {
        num_cores = std::strtol(argv[1], nullptr, 10);
    }

    if (num_cores <= 0) {
        std::cerr << "Invalid number of cores: " << num_cores << std::endl;
        return EXIT_FAILURE;
    }

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <num_cores> <output_file>\n";
        return EXIT_FAILURE;
    }

    std::vector<std::thread> threads;
    std::vector<std::unique_ptr<num_t>> results(num_cores);

    std::random_device rd;
    unsigned int seed = rd();

    auto exec_start = std::chrono::steady_clock::now();
    unsigned int N = 1000000;

    std::mt19937_64 gen(static_cast<num_t>(seed));

    #ifdef USE_INT64
    constexpr int64_t max_val = std::numeric_limits<int64_t>::max() / 10000000000.0;
    constexpr int64_t min_val = -max_val;
    std::uniform_int_distribution<int64_t> dist(min_val, max_val);
    std::vector<int64_t> numbers(N);
    for (size_t i = 0; i < N; ++i) {
        numbers[i] = dist(gen);
    }
#elif defined(USE_FLOAT) || defined(USE_DOUBLE)
    constexpr num_t max_val = std::numeric_limits<num_t>::max() / 1000000000.0;
    constexpr num_t min_val = -max_val;
    std::uniform_real_distribution<num_t> dist(min_val, max_val);
    std::vector<num_t> numbers(N);
    for (size_t i = 0; i < N; ++i) {
        numbers[i] = dist(gen);
    }
    num_t divisor = dist(gen);
    if (divisor == 0.0) divisor = 1.0;
#elif defined(USE_MATMUL)
    N = 1000;
    std::vector<std::vector<num_t>> A(N, std::vector<num_t>(N));
    std::vector<std::vector<num_t>> B(N, std::vector<num_t>(N));
    std::uniform_real_distribution<num_t> dist(-10000.0, 10000.0);

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            A[i][j] = dist(gen);
            B[i][j] = dist(gen);
        }
    }
#elif defined(USE_ROOT)
    NUM_ITER = 100000000;
    constexpr num_t max_val = std::numeric_limits<num_t>::max() / 100000000.0;
    constexpr num_t min_val = -max_val;
    std::uniform_real_distribution<num_t> dist(min_val, max_val);
    std::vector<num_t> numbers(N);
    for (size_t i = 0; i < N; ++i) {
        numbers[i] = dist(gen);
    }
    num_t divisor = dist(gen);
    if (divisor == 0.0) divisor = 1.0;
#endif

for (long i = 0; i < num_cores; ++i) {
        results[i] = std::make_unique<num_t>(0.0);
#ifdef USE_INT64
        threads.emplace_back(calc, std::cref(numbers), results[i].get());
#elif defined(USE_FLOAT) || defined(USE_DOUBLE)
        threads.emplace_back(calc, std::cref(numbers), divisor, results[i].get());
#elif defined(USE_MATMUL)
        threads.emplace_back(mat_mul, results[i].get(), std::cref(A), std::cref(B));
#elif defined(USE_ROOT)
        threads.emplace_back(calc, std::cref(numbers), results[i].get());
#endif
    }

    for (auto& t : threads) {
        t.join();
    }
    auto exec_end = std::chrono::steady_clock::now();
    double exec_time = std::chrono::duration<double>(exec_end - exec_start).count();

    bool all_equal = true;
    num_t first_value = *results[0];
    for (long i = 0; i < num_cores; ++i) {
        // std::cout << "Thread " << i << " returned: " << *results[i] << std::endl;
        if (*results[i] != first_value) {
            all_equal = false;
        }
    }

    std::string out_filename = argv[2];
    std::ofstream csv(out_filename, std::ios::app);

    csv << 1900+ltm->tm_year << "-"
        << std::setw(2) << std::setfill('0') << 1+ltm->tm_mon << "-"
        << std::setw(2) << std::setfill('0') << ltm->tm_mday << " "
        << std::setw(2) << std::setfill('0') << ltm->tm_hour << ":"
        << std::setw(2) << std::setfill('0') << ltm->tm_min << ":"
        << std::setw(2) << std::setfill('0') << ltm->tm_sec << ","
        << TEST_NAME << "," << exec_time << "," << seed << "," << (all_equal ? 1 : 0);

    for (long i = 0; i < num_cores; ++i) {
        csv << "," << *results[i];
    }

    csv << "\n";
    csv.close();

    if (all_equal) {
        // std::cout << "All threads returned the same value: " << first_value << std::endl;
    } else {
        std::cout << "Not all threads returned the same value." << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}