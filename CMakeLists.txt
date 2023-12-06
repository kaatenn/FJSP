cmake_minimum_required(VERSION 3.26)
project(fjsp)

set(CMAKE_CXX_STANDARD 17)

file(GLOB nlohmann_header "${PROJECT_SOURCE_DIR}/nlohmann/*.hpp")
file(GLOB nlohmann_source "${PROJECT_SOURCE_DIR}/nlohmann/*.cpp")

list(APPEND NOLHMANN_FILES ${nlohmann_header} ${nlohmann_source})

add_executable(fjsp main.cpp
        Job.cpp
        Job.h
        ${NOLHMANN_FILES}
        Chromosome.cpp
        Chromosome.h
)

target_include_directories(fjsp PRIVATE ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/nlohmann)
