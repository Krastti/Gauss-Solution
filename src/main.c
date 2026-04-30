#include "experiments/experiments.h"
#include "gauss/gauss.h"
#include "logging/logger.h"
#include "lu/lu.h"

#include <stdio.h>
#include <stdlib.h>

#ifdef _WIN32
#include <windows.h>
#endif

#define MAX_MATRIX_SIZE 10

int main(void) {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif
    logger_set_level(LOG_ERROR);

    int choice = -1;

    while (choice != 0) {
        printf("Меню для запуска разных СЛАУ\n");
        printf("1. Решить методом Гаусса без выбора главного элемента\n");
        printf("2. Решить методом Гаусса с частичным выбором главного элемента\n");
        printf("3. Решить через LU-разложение\n");
        printf("4. Запустить эксперименты и сохранить CSV\n");
        printf("5. Изменить уровень логирования\n");
        printf("0. Выход\n");
        printf("Ваш выбор: ");

        if (scanf("%d", &choice) != 1) {
            printf("Ошибка ввода.\n");
            return 1;
        }

        if (choice == 0) {
            printf("Работа завершена.\n");
        } else if (choice == 4) {
            if (experiments_run_all() == 0) {
                printf("Эксперименты завершены. CSV-файлы сохранены в data/reports.\n\n");
            } else {
                printf("Не удалось записать CSV-файлы.\n\n");
            }
        } else if (choice == 5) {
            int level_choice = -1;
            printf("Уровни логирования:\n");
            printf("0. Debug\n");
            printf("1. Info\n");
            printf("2. Warning\n");
            printf("3. Error\n");
            printf("Установить уровень: ");

            if (scanf("%d", &level_choice) != 1) {
              printf("Ошибка ввода.\n");
              return 1;
            }

            if (level_choice == logger_current_level()) {
              printf("Выбранный уровень уже установлен!\n");
            } else if (level_choice == 1) {
              logger_set_level(level_choice);
            } else if (level_choice == 2) {
              logger_set_level(level_choice);
            } else if (level_choice == 3) {
              logger_set_level(level_choice);
            } else if (level_choice == 4) {
              logger_set_level(level_choice);
            }
        } else if (choice >= 1 && choice <= 3) {
            size_t n;

            printf("Введите размер системы n: ");
            if (scanf("%zu", &n) != 1 || n == 0 || n > MAX_MATRIX_SIZE) {
                printf("Некорректный размер системы. Размер должен быть от 1 до %d.\n", MAX_MATRIX_SIZE);
                return 1;
            }

            double* A = gauss_alloc_matrix(n);
            double* b = calloc(n, sizeof(double));
            double* x = calloc(n, sizeof(double));

            if (A == NULL || b == NULL || x == NULL) {
                printf("Не удалось выделить память.\n");
                gauss_free_matrix(A);
                free(b);
                free(x);
                return 1;
            }

            printf("Введите коэффициенты матрицы A:\n");
            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    printf("A[%zu][%zu] = ", i + 1, j + 1);
                    scanf("%lf", gauss_at_ptr(A, n, i, j));
                }
            }

            printf("Введите правую часть b:\n");
            for (size_t i = 0; i < n; i++) {
                printf("b[%zu] = ", i + 1);
                scanf("%lf", &b[i]);
            }

            double total_ms = 0.0;
            double decompose_ms = 0.0;
            double solve_ms = 0.0;
            int gauss_status = GAUSS_OK;
            int lu_status = LU_OK;

            if (choice == 1) {
                printf("\nМетод: Гаусс без выбора главного элемента\n");
                gauss_status = gauss_solve(n, A, b, x, GAUSS_METHOD_CLASSIC, &total_ms);
            } else if (choice == 2) {
                printf("\nМетод: Гаусс с частичным выбором главного элемента\n");
                gauss_status = gauss_solve(n, A, b, x, GAUSS_METHOD_PIVOT, &total_ms);
            } else {
                printf("\nМетод: LU-разложение\n");
                lu_status = lu_decompose_solve(n, A, b, x, &decompose_ms, &solve_ms, &total_ms);
            }

            if ((choice == 3 && lu_status == LU_OK) || (choice != 3 && gauss_status == GAUSS_OK)) {
                gauss_print_vector(n, x, "Решение x");
                printf("Невязка ||Ax - b|| = %.12e\n", gauss_residual(n, A, b, x));
                printf("Общее время: %.6f мс\n", total_ms);

                if (choice == 3) {
                    printf("LU-разложение: %.6f мс\n", decompose_ms);
                    printf("Решение после разложения: %.6f мс\n", solve_ms);
                }
            } else {
                printf("Матрица вырождена или близка к вырожденной.\n");
            }
            printf("\n");

            gauss_free_matrix(A);
            free(b);
            free(x);
        } else {
            printf("Неизвестный пункт меню.\n\n");
        }
    }

    return 0;
}
