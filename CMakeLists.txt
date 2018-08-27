cmake_minimum_required(VERSION 3.0)
set(CMAKE_BUILD_TYPE Debug)
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
    add_executable(test_genomic_location test/test_genomic_location.cpp)
    add_executable(test_replication_fork test/test_replication_fork.cpp)

    target_link_libraries(test_chromosome deps gtest)
    target_link_libraries(test_genome deps gtest)
    target_link_libraries(test_genomic_location deps gtest)
    target_link_libraries(test_replication_fork deps gtest)

    gtest_discover_tests(test_chromosome)
    gtest_discover_tests(test_genome)
    gtest_discover_tests(test_genomic_location)
    gtest_discover_tests(test_replication_fork)
endif()
