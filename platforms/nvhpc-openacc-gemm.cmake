set(CMAKE_Fortran_COMPILER      "mpif90")
set(CMAKE_C_COMPILER            "mpicc")
set(OPENMP_FLAGS                "-Mnoopenmp")

set(General_Fortran_FLAGS       "-Wall -fstrict-aliasing")
set(General_C_FLAGS             "-Wall -alias=ansi")

# OpenACC + cuBLAS (Fortran-only batched-GEMM path for pseudo-pt, see
# salmon-gpu-optimization-ideas idea 2), NOT the hand-written CUDA kernels
# -- same -cuda/-cudalib flags as nvhpc-openacc-cuda.cmake (needed for the
# `use cublas` Fortran module and cublasZgemmStridedBatched), but
# -DUSE_GEMM instead of -DUSE_CUDA, and USE_CUDA left OFF below so
# src/common/CMakeLists.txt doesn't add zpseudo.cu/stencil_current.cu to
# the build at all.
set(OpenACC_FLAGS               "-acc=strict -gpu=cc80,cc90,cc100,cc120,managed,ptxinfo -cudalib=cublas,cusolver -cuda -Minfo=accel -DUSE_OPENACC -DUSE_GEMM")

set(CMAKE_Fortran_FLAGS_DEBUG   "-O2 -g -traceback ${General_Fortran_FLAGS} ${OpenACC_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG         "-O2 -g -traceback ${General_C_FLAGS} ${OpenACC_FLAGS}")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 ${General_Fortran_FLAGS} ${OpenACC_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE       "-O3 ${General_C_FLAGS} ${OpenACC_FLAGS}")

set(USE_MPI_DEFAULT             ON)
set(USE_OPENACC                 ON)
set(USE_CUDA                    OFF)

########
# CMake Platform-specific variables
########
set(CMAKE_SYSTEM_NAME "Linux" CACHE STRING "Compiling for x86_64 + Nvidia GPU")
set(CMAKE_SYSTEM_PROCESSOR "openacc")
