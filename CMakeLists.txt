cmake_minimum_required(VERSION 3.16)
project(Algebra1 C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(ALGEBRA1_SOURCES
        src/gauss/gauss.c
        src/logging/logger.c
        src/timer/timer.c
        src/matgen/matgen.c
        src/lu/lu.c
)

add_executable(Algebra1
        src/main.c
        ${ALGEBRA1_SOURCES}
)

target_compile_definitions(Algebra1 PRIVATE
        DATA_DIR="${CMAKE_SOURCE_DIR}/data"
        REPORTS_DIR="${CMAKE_SOURCE_DIR}/data/reports"
)

option(ALGEBRA1_BUILD_TESTS "Build Algebra1 unit tests" ON)

if(ALGEBRA1_BUILD_TESTS)
    enable_testing()

    add_executable(Algebra1Tests
            test/test_algebra1.c
            ${ALGEBRA1_SOURCES}
    )

    add_test(NAME Algebra1Tests COMMAND Algebra1Tests)
endif()
