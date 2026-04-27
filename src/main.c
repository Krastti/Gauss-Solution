#include "gauss/gauss.h"
#include "logging/logger.h"
#include "lu/lu.h"
#include "matgen/matgen.h"

#include <direct.h>
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>

#ifndef DATA_DIR
#define DATA_DIR "data"
#endif

#ifndef REPORTS_DIR
#define REPORTS_DIR "data/reports"
#endif

#define REPORT_SINGLE_PATH REPORTS_DIR "/single_system_experiment.csv"
#define REPORT_MULTI_PATH REPORTS_DIR "/multiple_rhs_experiment.csv"
#define REPORT_HILBERT_PATH REPORTS_DIR "/hilbert_accuracy_experiment.csv"

#define ARRAY_COUNT(array) (sizeof(array) / sizeof((array)[0]))

static void close_file(FILE** file) {
    if (*file) {
        fclose(*file);
        *file = NULL;
    }
}

static void close_reports(FILE** single, FILE** multi, FILE** hilbert) {
    close_file(single);
    close_file(multi);
    close_file(hilbert);
}

static void ensure_report_file(const char* path) {
    FILE* file = fopen(path, "a");
    if (!file) {
        LOG_ERROR("ensure_report_file: не удалось создать файл отчёта: %s", path);
        return;
    }

    LOG_INFO("ensure_report_file: файл отчёта существует или создан: %s", path);
    fclose(file);
}

static void ensure_report_directories(void) {
    LOG_INFO("ensure_report_directories: подготовка папок %s и %s", DATA_DIR, REPORTS_DIR);
    _mkdir(DATA_DIR);
    _mkdir(REPORTS_DIR);

    ensure_report_file(REPORT_SINGLE_PATH);
    ensure_report_file(REPORT_MULTI_PATH);
    ensure_report_file(REPORT_HILBERT_PATH);
}

static void fill_ones(size_t n, double* x) {
    LOG_DEBUG("fill_ones: заполнение вектора единицами, n=%zu", n);
    for (size_t i = 0; i < n; i++)
        x[i] = 1.0;
}

static void write_timing_row(FILE* out, size_t n, const char* method, int status,
                             double total_ms, double decompose_ms, double solve_ms,
                             const double* A, const double* b, const double* x) {
    const double residual = status == ALGEBRA_OK ? gauss_residual(n, A, b, x) : -1.0;

    fprintf(out, "%zu,%s,%d,%.6f,%.6f,%.6f,%.12e\n",
            n, method, status, total_ms, decompose_ms, solve_ms, residual);
}

static void write_accuracy_row(FILE* out, size_t n, const char* method, int status,
                               const double* A, const double* b,
                               const double* x, const double* exact) {
    const double relative_error = status == ALGEBRA_OK ? matgen_relative_error(n, x, exact) : -1.0;
    const double residual = status == ALGEBRA_OK ? gauss_residual(n, A, b, x) : -1.0;

    fprintf(out, "%zu,%s,%d,%.12e,%.12e\n",
            n, method, status, relative_error, residual);
}

static void write_gauss_solution(FILE* out, size_t n, const char* method_name,
                                 int method, const double* A,
                                 const double* b, double* x) {
    double elapsed = 0.0;
    const int status = gauss_solve(n, A, b, x, method, &elapsed);

    write_timing_row(out, n, method_name, status, elapsed, 0.0, elapsed, A, b, x);
}

static void write_lu_solution(FILE* out, size_t n, const double* A,
                              const double* b, double* x) {
    double decompose_ms = 0.0;
    double solve_ms = 0.0;
    double total_ms = 0.0;
    const int status = lu_decompose_solve(n, A, b, x, &decompose_ms, &solve_ms, &total_ms);

    write_timing_row(out, n, "lu", status, total_ms, decompose_ms, solve_ms, A, b, x);
}

static void write_gauss_accuracy(FILE* out, size_t n, const char* method_name,
                                 int method, const double* A, const double* b,
                                 double* x, const double* exact) {
    const int status = gauss_solve(n, A, b, x, method, NULL);

    write_accuracy_row(out, n, method_name, status, A, b, x, exact);
}

static void write_lu_accuracy(FILE* out, size_t n, const double* A,
                              const double* b, double* x, const double* exact) {
    const int status = lu_decompose_solve(n, A, b, x, NULL, NULL, NULL);

    write_accuracy_row(out, n, "lu", status, A, b, x, exact);
}

static void write_single_system_experiment(FILE* out) {
    LOG_INFO("%s", "write_single_system_experiment: начало эксперимента для одной системы");
    const size_t sizes[] = {100, 200, 500, 1000};

    fprintf(out, "n,method,status,total_ms,decompose_ms,solve_ms,residual\n");

    for (size_t i = 0; i < ARRAY_COUNT(sizes); i++) {
        const size_t n = sizes[i];
        LOG_INFO("write_single_system_experiment: обработка n=%zu", n);
        double* A = matgen_random(n, -1.0, 1.0, 1000u + (unsigned int)n);
        double* b = matgen_random_vector(n, -1.0, 1.0, 2000u + (unsigned int)n);
        double* x = calloc(n, sizeof(double));

        if (!A || !b || !x) {
            LOG_ERROR("write_single_system_experiment: ошибка выделения памяти для n=%zu", n);
            fprintf(out, "%zu,allocation,alloc_error,0,0,0,0\n", n);
            gauss_free_matrix(A);
            free(b);
            free(x);
            continue;
        }

        write_gauss_solution(out, n, "gauss_classic", GAUSS_METHOD_CLASSIC, A, b, x);
        write_gauss_solution(out, n, "gauss_pivot", GAUSS_METHOD_PIVOT, A, b, x);
        write_lu_solution(out, n, A, b, x);

        gauss_free_matrix(A);
        free(b);
        free(x);
        LOG_INFO("write_single_system_experiment: n=%zu завершён", n);
    }
    LOG_INFO("%s", "write_single_system_experiment: эксперимент завершён");
}

static void write_multiple_rhs_experiment(FILE* out) {
    LOG_INFO("%s", "write_multiple_rhs_experiment: начало эксперимента с несколькими правыми частями");
    const size_t n = 500;
    const size_t rhs_counts[] = {1, 10, 100};

    fprintf(out, "n,k,method,status,total_ms,decompose_ms,solve_ms\n");

    double* A = matgen_random(n, -1.0, 1.0, 3000u);
    double* L = gauss_alloc_matrix(n);
    double* U = gauss_alloc_matrix(n);
    double* b = NULL;
    double* x = calloc(n, sizeof(double));

    if (!A || !L || !U || !x) {
        LOG_ERROR("%s", "write_multiple_rhs_experiment: ошибка выделения памяти");
        fprintf(out, "%zu,0,allocation,alloc_error,0,0,0\n", n);
        gauss_free_matrix(A);
        gauss_free_matrix(L);
        gauss_free_matrix(U);
        free(x);
        return;
    }

    double decompose_ms = 0.0;
    int lu_status = lu_decompose(n, A, L, U, &decompose_ms);

    for (size_t i = 0; i < ARRAY_COUNT(rhs_counts); i++) {
        const size_t k = rhs_counts[i];
        LOG_INFO("write_multiple_rhs_experiment: обработка k=%zu", k);
        double gauss_total_ms = 0.0;
        double lu_solve_total_ms = 0.0;
        int gauss_status = ALGEBRA_OK;
        int solve_status = lu_status;

        for (size_t rhs = 0; rhs < k; rhs++) {
            b = matgen_random_vector(n, -1.0, 1.0, 4000u + (unsigned int)(rhs + k));
            if (!b) {
                LOG_ERROR("write_multiple_rhs_experiment: ошибка генерации правой части #%zu для k=%zu", rhs, k);
                gauss_status = ALGEBRA_ALLOC_ERR;
                solve_status = ALGEBRA_ALLOC_ERR;
                break;
            }

            double elapsed = 0.0;
            gauss_status = gauss_solve(n, A, b, x, GAUSS_METHOD_PIVOT, &elapsed);
            gauss_total_ms += elapsed;

            if (solve_status == ALGEBRA_OK) {
                elapsed = 0.0;
                solve_status = lu_solve(n, L, U, b, x, &elapsed);
                lu_solve_total_ms += elapsed;
            }

            free(b);
            b = NULL;
        }

        fprintf(out, "%zu,%zu,gauss_pivot,%d,%.6f,0,%.6f\n",
                n, k, gauss_status, gauss_total_ms, gauss_total_ms);
        fprintf(out, "%zu,%zu,lu,%d,%.6f,%.6f,%.6f\n",
                n, k, solve_status, decompose_ms + lu_solve_total_ms,
                decompose_ms, lu_solve_total_ms);
        LOG_INFO("write_multiple_rhs_experiment: k=%zu завершён", k);
    }

    gauss_free_matrix(A);
    gauss_free_matrix(L);
    gauss_free_matrix(U);
    free(b);
    free(x);
    LOG_INFO("%s", "write_multiple_rhs_experiment: эксперимент завершён");
}

static void write_hilbert_experiment(FILE* out) {
    LOG_INFO("%s", "write_hilbert_experiment: начало эксперимента на матрицах Гильберта");
    const size_t sizes[] = {5, 10, 15};

    fprintf(out, "n,method,status,relative_error,residual\n");

    for (size_t i = 0; i < ARRAY_COUNT(sizes); i++) {
        const size_t n = sizes[i];
        LOG_INFO("write_hilbert_experiment: обработка n=%zu", n);
        double* H = matgen_hilbert(n);
        double* exact = malloc(n * sizeof(double));
        double* x = calloc(n, sizeof(double));
        if (!H || !exact || !x) {
            LOG_ERROR("write_hilbert_experiment: ошибка выделения памяти для n=%zu", n);
            fprintf(out, "%zu,allocation,alloc_error,-1,-1\n", n);
            gauss_free_matrix(H);
            free(exact);
            free(x);
            continue;
        }

        fill_ones(n, exact);
        double* b = matgen_rhs_from_exact(n, H, exact);
        if (!b) {
            LOG_ERROR("write_hilbert_experiment: ошибка вычисления правой части для n=%zu", n);
            fprintf(out, "%zu,allocation,alloc_error,-1,-1\n", n);
            gauss_free_matrix(H);
            free(exact);
            free(x);
            continue;
        }

        write_gauss_accuracy(out, n, "gauss_classic", GAUSS_METHOD_CLASSIC, H, b, x, exact);
        write_gauss_accuracy(out, n, "gauss_pivot", GAUSS_METHOD_PIVOT, H, b, x, exact);
        write_lu_accuracy(out, n, H, b, x, exact);

        gauss_free_matrix(H);
        free(exact);
        free(b);
        free(x);
        LOG_INFO("write_hilbert_experiment: n=%zu завершён", n);
    }
    LOG_INFO("%s", "write_hilbert_experiment: эксперимент завершён");
}

int main(void) {
    SetConsoleOutputCP(CP_UTF8);

    logger_set_level(LOG_ERROR);
    LOG_INFO("%s", "main: запуск программы лабораторной работы");
    ensure_report_directories();

    LOG_INFO("%s", "main: открытие CSV-файлов отчётов");
    FILE* single = fopen(REPORT_SINGLE_PATH, "w");
    FILE* multi = fopen(REPORT_MULTI_PATH, "w");
    FILE* hilbert = fopen(REPORT_HILBERT_PATH, "w");

    if (!single || !multi || !hilbert) {
        LOG_ERROR("%s", "Не удалось открыть файлы отчётов в data/reports!");
        close_reports(&single, &multi, &hilbert);
        return 1;
    }

    write_single_system_experiment(single);
    write_multiple_rhs_experiment(multi);
    write_hilbert_experiment(hilbert);

    close_reports(&single, &multi, &hilbert);

    LOG_INFO("%s", "main: все отчёты успешно записаны");
    printf("Эксперименты завершены. CSV-файлы сохранены в data/reports.\n");
    return 0;
}
