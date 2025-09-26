cmake_minimum_required(VERSION 3.24)
project(J_Matrix)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
target_link_libraries(${PROJECT_NAME} PRIVATE ${OpenBLAS_LIBRARIES})