cmake_minimum_required(VERSION 3.14)
project(fastsos)

set(CMAKE_CXX_STANDARD 14)

# MOSEK

set (MOSEK_DIR /Users/joe/solvers/mosek/9.0)

include_directories(${MOSEK_DIR}/tools/platform/osx64x86/h)
link_directories(${MOSEK_DIR}/tools/platform/osx64x86/bin)
set (MOSEK_SRCS build_mosek.cpp build_mosek.h)

# SCS

set (SCS_DIR /Users/joe/solvers/scs)

include_directories(${SCS_DIR}/include)
include_directories(${SCS_DIR}/linsys)
include_directories(${SCS_DIR}/test)
include_directories(${SCS_DIR}/test/problems)
link_directories(${SCS_DIR}/out)
set (SCS_SRCS build_scs.cpp build_scs.h)

# SOSADMM

set (LAPACK_DIR /usr/lib)
include_directories(${LAPACK_DIR})
include_directories(sosadmm_src/h_files)
set (SOSADMM_SRCS sosadmm_src/h_files/projection.h sosadmm_src/cpp_files/projection.cpp sosadmm_src/cpp_files/sdp.cpp sosadmm_src/h_files/sdp.h sosadmm_src/cpp_files/ordering.cpp sosadmm_src/h_files/ordering.h sosadmm_src/cpp_files/admm.cpp sosadmm_src/h_files/admm.h sosadmm_src/cpp_files/vecop.cpp sosadmm_src/h_files/vecop.h sosadmm_src/speedTests/vecopTests.h sosadmm_src/cpp_files/h_ops.cpp sosadmm_src/h_files/h_ops.h)

# Build list of sources according to what's available

set (SRCS main.cpp manip_poly.cpp manip_poly.h sos.cpp sos.h tests.cpp tests.h timing.cpp timing.h)

add_executable(fastsos ${SRCS} ${MOSEK_SRCS} ${SCS_SRCS} ${SOSADMM_SRCS})

target_link_libraries(fastsos mosek64.9.0 fusion64.9.0 libscsdir.dylib liblapack.dylib)