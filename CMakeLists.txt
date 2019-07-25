cmake_minimum_required(VERSION 3.14)
project(fastsos)

set(CMAKE_CXX_STANDARD 14)

include_directories(/usr/local/include)
include_directories(/Users/joe/solvers/mosek/9.0/tools/platform/osx64x86/h)
link_directories(/usr/local/lib)
link_directories(/Users/joe/solvers/mosek/9.0/tools/platform/osx64x86/bin)

include_directories(/Users/joe/solvers/scs/include)
include_directories(/Users/joe/solvers/scs/linsys)
include_directories(/Users/joe/solvers/scs/test)
include_directories(/Users/joe/solvers/scs/test/problems)
link_directories(/Users/joe/solvers/scs/out)

add_executable(fastsos main.cpp polystring.cpp polystring.h sos.cpp sos.h tests.cpp tests.h build_mosek.cpp build_mosek.h build_scs.cpp build_scs.h)

target_link_libraries(fastsos mosek64.9.0 fusion64.9.0 libscsindir.dylib)