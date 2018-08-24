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
add_library(genome src/genome.cpp)
add_library(genomic_location src/genomic_location.cpp)
add_library(util src/util.cpp)
add_library(replication_fork src/replication_fork.cpp)
add_library(fork_manager src/fork_manager.cpp)
add_library(deps src/chromosome.cpp src/genome.cpp src/genomic_location.cpp src/util.cpp src/replication_fork.cpp)

#######
# Testing
#######
if(BUILD_TESTING)
    enable_testing()

    add_executable(test_chromosome test/test_chromosome.cpp)
    add_executable(test_genome test/test_genome.cpp)

    target_link_libraries(test_chromosome chromosome util gtest)
    target_link_libraries(test_genome chromosome genome genomic_location util gtest)

    gtest_discover_tests(test_chromosome)
    gtest_discover_tests(test_genome)
endif()
