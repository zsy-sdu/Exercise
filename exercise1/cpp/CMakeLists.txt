cmake_minimum_required(VERSION 3.24)
project(mol_class)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add source files
add_library(molecular_lib molecular.cpp molecular.h)

# Add executable
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE molecular_lib)