cmake_minimum_required(VERSION 3.0)
project(ReDyMo)
include(CTest)
include(GoogleTest)

SET(CUSTOM_BUILD_TYPE "" CACHE STRING "Sets custom comp flags")

if(CUSTOM_BUILD_TYPE STREQUAL Debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp -g -pg")

elseif(CUSTOM_BUILD_TYPE STREQUAL Performance)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp -Ofast -march=native -frename-registers")

elseif(CUSTOM_BUILD_TYPE STREQUAL Development)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp -g")

else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp -O3 -frename-registers")

endif()

#######
# OpenCL
#######
find_package(OpenCL REQUIRED)
include_directories(${OpenCL_INCLUDE_DIRS})
link_directories(${OpenCL_LIBRARY})

#######
# Rapid YAML
#######
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/thirdparty/rapidyaml ryml)
include_directories(${CMAKE_CURRENT_LIST_DIR}/thirdparty/rapidyaml/ext/c4core/src)
include_directories(${CMAKE_CURRENT_LIST_DIR}/thirdparty/rapidyaml/src)

#######
# Add the include folder
#######
include_directories(include)

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/thirdparty/googletest)
include_directories(${gtest_SOURCE_DIR}/include)
include_directories(${gmock_SOURCE_DIR}/include)

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/thirdparty/SQLiteCpp)
include_directories(
  ${CMAKE_CURRENT_LIST_DIR}/thirdparty/SQLiteCpp/include
)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/thirdparty/zstd/build/cmake)
include_directories(${CMAKE_CURRENT_LIST_DIR}/thirdparty/zstd/lib)

#######
# Compile componenets in libraries
#######
add_library(
    deps
    src/chromosome.cpp
    src/data_manager.cpp
    src/fork_manager.cpp
    src/genome.cpp
    src/genomic_location.cpp
    src/replication_fork.cpp
    src/util.cpp
    src/s_phase.cpp
    
    src/configuration.cpp
    src/evolution_data_provider.cpp
    src/evolution.cpp
#    src/gpu_s_phase.cpp
#    src/gpu_chromosome.cpp
)
add_executable(simulator src/main.cpp)
target_link_libraries(simulator deps gtest SQLiteCpp libzstd_static sqlite3 pthread dl ryml ${OpenCL_LIBRARY})

#######
# Testing
#######
if(BUILD_TESTING)
    enable_testing()

    add_executable(test_chromosome test/test_chromosome.cpp)
    add_executable(test_genome test/test_genome.cpp)
    add_executable(test_genomic_location test/test_genomic_location.cpp)
    add_executable(test_replication_fork test/test_replication_fork.cpp)
    add_executable(test_fork_manager test/test_fork_manager.cpp)
    add_executable(test_data_manager test/test_data_manager.cpp)

    target_link_libraries(test_chromosome deps gtest gmock)
    target_link_libraries(test_genome deps libzstd_static gtest)
    target_link_libraries(test_genomic_location deps libzstd_static gtest)
    target_link_libraries(test_replication_fork deps libzstd_static gtest)
    target_link_libraries(test_fork_manager deps libzstd_static gtest)
    target_link_libraries(test_data_manager deps gtest SQLiteCpp libzstd_static sqlite3 pthread dl)

    gtest_discover_tests(test_chromosome)
    gtest_discover_tests(test_genome)
    gtest_discover_tests(test_genomic_location)
    gtest_discover_tests(test_replication_fork)
    gtest_discover_tests(test_fork_manager)
    gtest_discover_tests(test_data_manager)
endif()
