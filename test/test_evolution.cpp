#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <omp.h>

#include "../include/evolution.hpp"
class MockEvolutionManager : public EvolutionManager
{
  public:
    MOCK_METHOD(void, simulate, (), (override));
    MOCK_METHOD(void, reproduce, (), (override));

    MockEvolutionManager(Configuration config, int seed)
        : EvolutionManager(config, seed)
    {
    }
    int get_current_generation() { return this->current_generation; }
};

class PublicEvolutionManager : public EvolutionManager
{
  public:
    PublicEvolutionManager(Configuration config, int seed)
        : EvolutionManager(config, seed)
    {
    }

    cl_configuration_data get_arguments() { return this->arguments; }
    std::vector<std::shared_ptr<EvolutionDataProvider>> get_data_providers()
    {
        return this->data_providers;
    }
    std::vector<std::vector<simulation_stats>> get_population()
    {
        return this->population;
    }
    void reproduce() {
        EvolutionManager::reproduce();
    }
};

class EvolutionTest : public ::testing::Test
{
  protected: // Runs before each test
    // Allow the test class to set args directly
    friend class Configuration;
    void SetUp() { omp_set_num_threads(1); }
};

TEST_F(EvolutionTest, CalculateFitness)
{
    instance_metrics metrics;
    cl_evolution_data config;
    double fitness;

    metrics.collisions          = 2000;
    config.fitness.max_coll_all = 0.5;
    config.fitness.min_coll_all = 0.5;

    fitness = calculate_fitness(metrics, config);
    ASSERT_NEAR(fitness, 0.5, 0.001);

    metrics.collisions          = 50000;
    config.fitness.max_coll_all = 1.0;
    config.fitness.min_coll_all = 0.0;

    fitness = calculate_fitness(metrics, config);
    ASSERT_NEAR(fitness, 1.0, 0.001);
}

TEST_F(EvolutionTest, TestConstructor)
{
    std::vector<char *> argv_mock = {"program_name",
                                     "--cells",
                                     "2",
                                     "--organism",
                                     "dummy",
                                     "--resources",
                                     "2",
                                     "--timeout",
                                     "5",
                                     "-C",
                                     "../test/config/config_evolution.yaml"};

    // Reset getopt global variable
    optind               = 1;
    Configuration config = Configuration(argv_mock.size(), argv_mock.data());
    PublicEvolutionManager evo(config, 0);
    ASSERT_EQ(evo.get_arguments(), config.arguments());
    ASSERT_FALSE(evo.get_data_providers().empty());
    ASSERT_FALSE(evo.get_population().empty());
}

TEST_F(EvolutionTest, TestGeneration)
{
    std::vector<char *> argv_mock = {
        "program_name", "--cells", "2",         "--organism", "dummy",
        "--resources",  "2",       "--timeout", "5"};

    // Reset getopt global variable
    optind               = 1;
    Configuration config = Configuration(argv_mock.size(), argv_mock.data());
    MockEvolutionManager evo(config, 0);
    int last_gen = evo.get_current_generation();
    EXPECT_CALL(evo, simulate());
    EXPECT_CALL(evo, reproduce());
    evo.generation();
    ASSERT_EQ(evo.get_current_generation(), last_gen + 1);
}

TEST_F(EvolutionTest, TestReproduce)
{

    std::vector<char *> argv_mock = {"program_name",
                                     "--cells",
                                     "2",
                                     "--organism",
                                     "dummy",
                                     "--resources",
                                     "2",
                                     "--timeout",
                                     "5",
                                     "-C",
                                     "../test/config/config_evolution.yaml"};

    // Reset getopt global variable
    optind               = 1;
    Configuration config = Configuration(argv_mock.size(), argv_mock.data());
    PublicEvolutionManager evo(config, 8);
    evo.reproduce();
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
