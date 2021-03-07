cmake_minimum_required(VERSION 3.10)
project(triplelens)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(src triplelens)
add_subdirectory(tests)
# Including googletest without googlemock
#add_subdirectory(tests/googletest/googletest)

add_subdirectory(tests/googletest)
