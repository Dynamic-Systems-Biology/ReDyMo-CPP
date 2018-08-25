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
add_library(deps src/chromosome.cpp src/data_manager.cpp src/fork_manager.cpp src/genome.cpp src/genomic_location.cpp src/replication_fork.cpp src/util.cpp)

#######
# Testing
#######
if(BUILD_TESTING)
    enable_testing()

    add_executable(test_chromosome test/test_chromosome.cpp)
    add_executable(test_genome test/test_genome.cpp)

    target_link_libraries(test_chromosome deps gtest)
    target_link_libraries(test_genome deps gtest)

    gtest_discover_tests(test_chromosome)
    gtest_discover_tests(test_genome)
endif()
