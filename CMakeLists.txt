cmake_minimum_required(VERSION 4.0)
project(GaussSolution C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(ALGEBRA1_SOURCES
        src/gauss/gauss.c
        src/logging/logger.c
        src/timer/timer.c
        src/matgen/matgen.c
        src/lu/lu.c
        src/experiments/experiments.c
)

add_executable(GaussSolution
        src/main.c
        ${ALGEBRA1_SOURCES}
)

target_compile_definitions(GaussSolution PRIVATE
        DATA_DIR="${CMAKE_SOURCE_DIR}/data"
        REPORTS_DIR="${CMAKE_SOURCE_DIR}/data/reports"
)
