cmake_minimum_required(VERSION 3.0)
project(ReDyMo)

#######
# Make sure libs are available locally
#######
execute_process(COMMAND git submodule update --init)
# We have to run these commands twice because of submodules that have
# submodules. They will be accounted for only after the first command.
execute_process(COMMAND git submodule foreach --recursive git submodule update --init)
execute_process(COMMAND git submodule foreach --recursive git submodule update --init)

#######
# Declare our custom options
#######
SET(CUSTOM_BUILD_TYPE "Release" CACHE STRING "Sets custom comp flags")
OPTION(BUILD_GPGPU "Wether to build the OpenCL dependant GPU code." OFF)
OPTION(BUILD_TESTING "Wether to build the tests." ON)
OPTION(COVERAGE "Wether to configure test coverage report." OFF)

#######
# Set flags for our custom build types
#######
if(CUSTOM_BUILD_TYPE STREQUAL Debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp -g -pg")
elseif(CUSTOM_BUILD_TYPE STREQUAL Performance)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp -Ofast -march=native -frename-registers")
    set(COVERAGE OFF)
elseif(CUSTOM_BUILD_TYPE STREQUAL Development)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp -g")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -fopenmp -O3")
endif()

#######
# OpenCL
#######
if (BUILD_GPGPU)
    find_package(OpenCL REQUIRED)
    include_directories(${OpenCL_INCLUDE_DIRS})
    link_directories(${OpenCL_LIBRARY})
endif()

#######
# Rapid YAML
#######
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/thirdparty/rapidyaml ryml)
include_directories(${CMAKE_CURRENT_LIST_DIR}/thirdparty/rapidyaml/ext/c4core/src)
include_directories(${CMAKE_CURRENT_LIST_DIR}/thirdparty/rapidyaml/src)

#######
# SQLite Cpp
#######
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/thirdparty/SQLiteCpp)
include_directories(
  ${CMAKE_CURRENT_LIST_DIR}/thirdparty/SQLiteCpp/include
)

#######
# Zstd
#######
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/thirdparty/zstd/build/cmake)
include_directories(${CMAKE_CURRENT_LIST_DIR}/thirdparty/zstd/lib)

#######
# Add the ReDyMo include folder
#######
include_directories(include)

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
)

if (BUILD_GPGPU)
    add_library(
        gpudeps
        src/gpu_s_phase.cpp
        src/gpu_chromosome.cpp
    )
endif()

message("")
message("================")
message("Current CXX flags: ${CMAKE_CXX_FLAGS}")
message("================")
message("")

add_executable(simulator src/main.cpp)
target_link_libraries(simulator deps SQLiteCpp libzstd_static sqlite3 pthread dl ryml)

if (BUILD_GPGPU)
    target_link_libraries(simulator gpudeps ${OpenCL_LIBRARY})
endif()

#######
# Testing
#######
if(BUILD_TESTING)
    include(CTest)
    include(GoogleTest)

    #######
    # Google Test
    #######
    add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/thirdparty/googletest)
    include_directories(${gtest_SOURCE_DIR}/include)
    include_directories(${gmock_SOURCE_DIR}/include)

    enable_testing()

    add_executable(test_chromosome test/test_chromosome.cpp)
    add_executable(test_genome test/test_genome.cpp)
    add_executable(test_genomic_location test/test_genomic_location.cpp)
    add_executable(test_replication_fork test/test_replication_fork.cpp)
    add_executable(test_fork_manager test/test_fork_manager.cpp)
    add_executable(test_data_manager test/test_data_manager.cpp)
    add_executable(test_util test/test_util.cpp)

    target_link_libraries(simulator deps gtest)

    target_link_libraries(test_chromosome deps gtest gmock gcov)
    target_link_libraries(test_genome deps libzstd_static gtest gcov)
    target_link_libraries(test_genomic_location deps libzstd_static gtest gcov)
    target_link_libraries(test_replication_fork deps libzstd_static gtest gcov)
    target_link_libraries(test_fork_manager deps libzstd_static gtest gcov)
    target_link_libraries(test_data_manager deps gtest SQLiteCpp libzstd_static sqlite3 pthread dl gcov)
    target_link_libraries(test_util deps gtest gmock gcov)

    gtest_discover_tests(test_chromosome)
    gtest_discover_tests(test_genome)
    gtest_discover_tests(test_genomic_location)
    gtest_discover_tests(test_replication_fork)
    gtest_discover_tests(test_fork_manager)
    gtest_discover_tests(test_data_manager)
    gtest_discover_tests(test_util)

    #######
    # Coverage report
    #######
    if(COVERAGE AND CMAKE_COMPILER_IS_GNUCXX)
        LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/thirdparty/CMakeModules")

        include(CodeCoverage)
        APPEND_COVERAGE_COMPILER_FLAGS()
        set(COVERAGE_GCOVR_EXCLUDES "thirdparty/*" "include/*" "test/*")

        # The ctest_no_fail.sh is used because the coverage tests won't
        # continue if the ctest returns non zero. Ctest will returns 0 only if
        # all tests pass.
        setup_target_for_coverage_gcovr_html(
            NAME coverage
            EXECUTABLE ${CMAKE_CURRENT_LIST_DIR}/script/ctest_no_fail.sh
            EXCLUDE "thirdparty/*"
            DEPENDENCIES test_chromosome test_genome test_genomic_location test_replication_fork test_fork_manager test_data_manager test_util
        )
        setup_target_for_coverage_gcovr_xml(
            NAME coverage_xml
            EXECUTABLE ${CMAKE_CURRENT_LIST_DIR}/script/ctest_no_fail.sh
            EXCLUDE "thirdparty/*"
            DEPENDENCIES test_chromosome test_genome test_genomic_location test_replication_fork test_fork_manager test_data_manager test_util
        )
    endif()
endif()
