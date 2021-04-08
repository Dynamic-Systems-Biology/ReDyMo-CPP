
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <iostream>
#include <memory>

#include "../include/configuration.hpp"

class ConfigurationTest : public ::testing::Test
{
    // Runs before each test
    void SetUp()
    {
        // Reset getopt's global variable so it can be used again in the next
        // test.
        optind = 1;
    }
};

std::vector<char *> argv_mock()
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

    std::vector<char *> argv_non_const;
    for (char *a : argv_mock)
    {
        char *b = (char *)malloc(strlen(a) * sizeof(char));
        for (int i = 0; i < strlen(a); i++)
            b[i] = a[i];
        b[strlen(a)] = 0;
        argv_non_const.push_back(b);
    }
    return argv_non_const;
}

TEST_F(ConfigurationTest, ValidCmdOptions)
{
    cl_configuration_data expected;
    expected.mode                                          = "basic";
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
    result = Configuration(argv_mock().size(), argv_mock().data()).arguments();
    ASSERT_EQ(expected, result);
}

TEST_F(ConfigurationTest, HelpCmdOption)
{
    std::vector<char *> argv_mock = {
        "program_name",
        "--h",
    };
    ASSERT_THROW(Configuration(argv_mock.size(), argv_mock.data()),
                 std::invalid_argument);
}

TEST_F(ConfigurationTest, InvalidCmdOption)
{
    std::vector<char *> argv_mock = {
        "program_name",
        "--z",
    };
    ASSERT_THROW(Configuration(argv_mock.size(), argv_mock.data()),
                 std::invalid_argument);
}

TEST_F(ConfigurationTest, MandatoryCmdOption)
{
    std::vector<char *> argv_mock = {
        "program_name",
    };
    ASSERT_THROW(Configuration(argv_mock.size(), argv_mock.data()),
                 std::invalid_argument);

    optind = 1;
    argv_mock.push_back("--cells");
    argv_mock.push_back("2");
    ASSERT_THROW(Configuration(argv_mock.size(), argv_mock.data()),
                 std::invalid_argument);

    optind = 1;
    argv_mock.push_back("--organism");
    argv_mock.push_back("my-organism");
    ASSERT_THROW(Configuration(argv_mock.size(), argv_mock.data()),
                 std::invalid_argument);

    optind = 1;
    argv_mock.push_back("--resources");
    argv_mock.push_back("5");
    ASSERT_THROW(Configuration(argv_mock.size(), argv_mock.data()),
                 std::invalid_argument);

    optind = 1;
    argv_mock.push_back("--timeout");
    argv_mock.push_back("10");
    ASSERT_NO_THROW(Configuration(argv_mock.size(), argv_mock.data()));
}

TEST_F(ConfigurationTest, ValidBasicConfigFile)
{
    cl_configuration_data expected;
    expected.mode                   = "basic";
    expected.cells                  = 100;
    expected.organism               = "TcruziCLBrenerEsmeraldo-like";
    expected.resources              = 50;
    expected.speed                  = 65;
    expected.timeout                = 1000000;
    expected.dormant                = true;
    expected.name                   = "abc";
    expected.period                 = 100;
    expected.constitutive           = 0;
    expected.data_dir               = "data_dir/";
    expected.probability            = 8;
    expected.output                 = "out_folder/";
    expected.threads                = 9;
    std::vector<char *> argv_config = {"program_name", "-C",
                                       "../test/config/config.yaml"};
    cl_configuration_data result =
        Configuration(argv_config.size(), argv_config.data()).arguments();

    // Seed is not set by config file
    expected.seed = result.seed;

    ASSERT_EQ(expected, result);
}

TEST_F(ConfigurationTest, InvalidBasicConfigFileOption)
{
    std::vector<char *> argv_config = {
        "program_name", "-C", "../test/config/invalid_option_config.yaml"};
    ASSERT_THROW(Configuration(argv_config.size(), argv_config.data()),
                 std::invalid_argument);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
