cmake_minimum_required(VERSION 3.14)
project(git_cppsos)

set(CMAKE_CXX_STANDARD 14)

include_directories(/usr/local/include)
include_directories(/Users/joe/solvers/mosek/9.0/tools/platform/osx64x86/h)
link_directories(/usr/local/lib)
link_directories(/Users/joe/solvers/mosek/9.0/tools/platform/osx64x86/bin)
add_executable(cppsos main.cpp polystring.cpp polystring.h sos.cpp sos.h tests.cpp tests.h)

target_link_libraries(cppsos mosek64.9.0 fusion64.9.0)