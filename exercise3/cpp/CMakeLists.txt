cmake_minimum_required(VERSION 3.24)
project(Huckel)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#OpenMP
find_package(OpenMP REQUIRED)
if(NOT OpenMP_FOUND)
    message(FATAL_ERROR "OpenMP not found. Please install OpenMP and try again.")
endif()

# Set OpenBLAS
if(APPLE)
    set(OpenBLAS_ROOT /opt/homebrew/opt/openblas)
    find_package(OpenBLAS REQUIRED)
if (NOT OpenBLAS_FOUND)
    message(FATAL_ERROR "OpenBLAS not found. Please install OpenBLAS and try again.")
endif()
else()
    find_package(OpenBLAS REQUIRED)
if (NOT OpenBLAS_FOUND)
    message(FATAL_ERROR "OpenBLAS not found. Please install OpenBLAS and try again.")
endif()
endif()

message(STATUS "Found OpenBLAS: ${OpenBLAS_LIBRARIES}")
include_directories(${OpenBLAS_INCLUDE_DIRS})

# Add executable
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE ${OpenBLAS_LIBRARIES} OpenMP::OpenMP_CXX)