cmake_minimum_required(VERSION 3.0)
project(ReDyMo)
include(CTest)
include(GoogleTest)

#######
# Add the include folder
#######
include_directories(include)

#######
# Compile componenets in libraries
#######
add_library(chromosome src/chromosome.cpp)

#######
# Testing
#######
if(BUILD_TESTING)
    enable_testing()

    add_executable(test_chromosome test/test_chromosome.cpp)

    target_link_libraries(test_chromosome chromosome gtest)

    gtest_discover_tests(test_chromosome)
endif()
