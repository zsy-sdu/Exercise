cmake_minimum_required(VERSION 3.28)
project(GMEV VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

message(STATUS "CMake version: ${CMAKE_VERSION}")
message(STATUS "Configuring project: ${PROJECT_NAME} version ${PROJECT_VERSION} with C++ Standard ${CMAKE_CXX_STANDARD}")
message(STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")

if(CMAKE_BUILD_TYPE STREQUAL Debug)
    add_compile_definitions(_DEBUG_)
elseif(CMAKE_BUILD_TYPE STREQUAL Release)
    add_compile_definitions(_RELEASE_)
else()
    add_compile_definitions(_DEBUG_)
endif()
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")



# Options
option(BUILD_EXE "Build executable" ON)
option(BUILD_LIB "Build Python library" ON)
#option(BUILD_CLIB "Build C++ library" OFF)
set(BLAS_LAPACK "OpenBLAS" CACHE STRING "Choose BLAS/LAPACK library: MKL, OpenBLAS, Accelerate")
set_property(CACHE BLAS_LAPACK PROPERTY STRINGS MKL OpenBLAS Accelerate)


# Set OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(OMP_FLAG ${OpenMP_CXX_FLAGS})
    set(OMP_LIB_NAME ${OpenMP_CXX_LIBRARIES})
else()
    message(FATAL_ERROR "OpenMP not found. Please install OpenMP or set CMAKE_CXX_COMPILER to a compiler that supports OpenMP.")
endif()

# Set BLAS/LAPACK
if(${BLAS_LAPACK} STREQUAL "OpenBLAS")
    add_compile_definitions(_USE_OPENBLAS_)
    if(UNIX AND NOT APPLE)
        find_package(BLAS REQUIRED)
        find_package(LAPACK REQUIRED)
        set(BLAS_LIBS ${BLAS_LIBRARIES})
        set(LAPACK_LIBS ${LAPACK_LIBRARIES})
    elseif(APPLE)
        # Directly specify Homebrew OpenBLAS
        set(OpenBLAS_ROOT "/opt/homebrew/opt/openblas")
        find_package(OpenBLAS REQUIRED)
        include_directories(${OpenBLAS_INCLUDE_DIRS})
        set(BLAS_LIBS ${OpenBLAS_LIBRARIES})
        set(LAPACK_LIBS ${OpenBLAS_LIBRARIES})
    else()
        message(FATAL_ERROR "OpenBLAS is not supported on this platform.")
    endif()
else()
    if(APPLE)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        message(FATAL_ERROR "GNU compiler is not supported with Accelerate on macOS. Please use Clang or AppleClang.")
    endif()
        find_package(BLAS REQUIRED)
        find_package(LAPACK REQUIRED)
        set(BLAS_LIBS ${BLAS_LIBRARIES})
        set(LAPACK_LIBS ${LAPACK_LIBRARIES})
        add_compile_definitions(_USE_ACCELERATE_) 
        add_compile_definitions(ACCELERATE_NEW_LAPACK)
        add_compile_definitions(ACCELERATE_LAPACK_ILP64)    
    else()
        message(FATAL_ERROR "Accelerate is only supported on macOS. Please select MKL or OpenBLAS.")
    endif()
endif()

# add libcint
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/libcint)
include_directories(${cint_BINARY_DIR}/include)
set(CINTLIB "${cint_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}cint${CMAKE_SHARED_LIBRARY_SUFFIX}")
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/libcint/include)

message(STATUS "OMP_LIB_NAME = ${OMP_LIB_NAME}")
message(STATUS "BLAS_LIBS = ${BLAS_LIBS}")
message(STATUS "LAPACK_LIBS = ${LAPACK_LIBS}")

# ---------------- CORE LIB ----------------
add_library(gmev_core_lib INTERFACE)
target_include_directories(gmev_core_lib INTERFACE src/core)
target_link_libraries(gmev_core_lib INTERFACE ${OMP_LIB_NAME} ${BLAS_LIBS} ${LAPACK_LIBS} ${CINTLIB})

# ---------------- EXE ----------------
if(BUILD_EXE)
    add_executable(gemv_exe src/main.cpp)
    target_link_libraries(gemv_exe PRIVATE gmev_core_lib)
endif()

# ---------------- PYTHON LIB ----------------
if(BUILD_LIB)
    find_package(pybind11 REQUIRED)
    pybind11_add_module(gmev_core src/pybind.cpp)
    target_link_libraries(gmev_core PRIVATE gmev_core_lib)
endif()

