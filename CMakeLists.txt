cmake_minimum_required(VERSION 3.10)
project(ecg_dealination C)

set(CMAKE_C_STANDARD 99)

set(SOURCES
    src/main.c
    src/ecg_processing.c
    src/csv_reader.c
    src/json_writer.c
)

include_directories(include)

add_executable(ecg_dealination ${SOURCES})

target_link_libraries(ecg_dealination m)
