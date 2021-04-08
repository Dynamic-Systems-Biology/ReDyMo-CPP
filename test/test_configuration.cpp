
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <iostream>
#include <memory>

#include "../include/configuration.hpp"

class ConfigurationTest : public ::testing::Test
{
    // Runs before each test
    void SetUp() {}

    // Runs after each test
    void TearDown() {}
};

TEST_F(ConfigurationTest, ValidCmdOptions)
{
    std::vector<char *> argv_mock = {
        "program_name",
        "--gpu",
        "false",
        "--cells",
        "2",
        "--organism",
        "dummy",
        "--resources",
        "2",
        "--speed",
        "3",
        "--dormant",
        "--summary",
        "--seed",
        "4",
        "--name",
        "abc",
        "--config",
        "config_file.yaml",
        "--timeout",
        "5",
        "--period",
        "6",
        "--constitutive",
        "7",
        "--data-dir",
        "data_dir/",
        "--probability",
        "8",
        "--output",
        "out_folder/",
        "--threads",
        "9",
    };
    cl_configuration_data expected;
    expected.mode                                          = "";
    expected.cells                                         = 2;
    expected.organism                                      = "dummy";
    expected.resources                                     = 2;
    expected.speed                                         = 3;
    expected.timeout                                       = 5;
    expected.dormant                                       = true;
    expected.seed                                          = 4;
    expected.name                                          = "abc";
    expected.period                                        = 6;
    expected.constitutive                                  = 7;
    expected.data_dir                                      = "data_dir/";
    expected.probability                                   = 8;
    expected.output                                        = "out_folder/";
    expected.threads                                       = 9;
    expected.evolution.population                          = 0;
    expected.evolution.generations                         = 0;
    expected.evolution.survivors                           = 0;
    expected.evolution.mutations.probability_landscape.add = 0;
    expected.evolution.mutations.probability_landscape.del = 0;
    expected.evolution.mutations.probability_landscape.change_mean.prob = 0;
    expected.evolution.mutations.probability_landscape.change_mean.std  = 0;
    expected.evolution.mutations.probability_landscape.change_std.prob  = 0;
    expected.evolution.mutations.probability_landscape.change_std.std   = 0;
    expected.evolution.mutations.probability_landscape.change_std.max   = 0;
    expected.evolution.mutations.genes.move.prob                        = 0;
    expected.evolution.mutations.genes.move.std                         = 0;
    expected.evolution.mutations.genes.swap.prob                        = 0;
    expected.evolution.fitness.min_sphase                               = 0;
    expected.evolution.fitness.match_mfaseq                             = 0;
    expected.evolution.fitness.max_coll_all                             = 0;
    expected.evolution.fitness.min_coll_all                             = 0;
    expected.evolution.fitness.max_coll.weight                          = 0;
    expected.evolution.fitness.max_coll.gene                            = "";
    expected.evolution.fitness.min_coll.weight                          = 0;
    expected.evolution.fitness.min_coll.gene                            = "";

    cl_configuration_data result;
    result = Configuration(argv_mock.size(), argv_mock.data()).arguments();
    ASSERT_EQ(expected, result);
}

TEST_F(ConfigurationTest, InvalidCmdOptions) { ASSERT_TRUE(true); }

TEST_F(ConfigurationTest, ValidConfigFile) { ASSERT_TRUE(true); }

TEST_F(ConfigurationTest, InvalidConfigFile) { ASSERT_TRUE(true); }

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
