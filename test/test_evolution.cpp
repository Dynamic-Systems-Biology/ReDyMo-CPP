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

    cl_configuration_data get_arguments() { return this->arguments; }
    int get_current_generation() { return this->current_generation; }
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
    std::vector<char *> argv_mock = {
        "program_name", "--cells", "2",         "--organism",
        "--resources",  "2",       "--timeout", "5"};

    // Reset getopt global variable
    optind               = 1;
    Configuration config = Configuration(argv_mock.size(), argv_mock.data());
    MockEvolutionManager evo(config, 0);
    ASSERT_EQ(evo.get_arguments(), config.arguments());
}

TEST_F(EvolutionTest, TestGeneration) {

}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
