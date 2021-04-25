#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "../include/s_phase.hpp"
#include <fstream>
#include <vector>

class SPhaseTest : public ::testing::Test
{

  protected:
    std::shared_ptr<ReplicationFork> fork;
    std::vector<std::shared_ptr<Chromosome>> chrms;
    unsigned long long seed = 0;

  protected:
    SPhaseTest() {}

    void SetUp()
    {
        // Reset getopt's global variable so it can be used again in the next
        // test.
        optind = 1;
    }

    void TearDown() {}
};

class PublicSPhase : public SPhase
{
  public:
    PublicSPhase(int origins_range, int n_resources, int replication_speed,
                 int timeout, int transcription_period, bool has_dormant,
                 std::shared_ptr<DataProvider> data, std::string organism,
                 std::string name, std::string output_folder = "output",
                 unsigned long long seed = 0)
        : SPhase(origins_range, n_resources, replication_speed, timeout,
                 transcription_period, has_dormant, data, organism, name,
                 output_folder, seed){};

    std::shared_ptr<Genome> get_genome() { return SPhase::genome; }
    std::shared_ptr<ForkManager> get_fork_manager()
    {
        return SPhase::fork_manager;
    }
};

std::vector<char *>
unconst_string_vector(std::vector<const char *> const_strings)
{
    std::vector<char *> strings_non_const;
    for (const char *a : const_strings)
    {
        char *b = (char *)malloc(strlen(a) * sizeof(char));
        for (unsigned long long i = 0; i < strlen(a); i++)
            b[i] = a[i];
        b[strlen(a)] = 0;
        strings_non_const.push_back(b);
    }
    return strings_non_const;
}

TEST_F(SPhaseTest, NoCollisionsWhenNoTranscription)
{
    std::vector<const char *> const_argv_mock = {"program_name",
                                                 "--gpu",
                                                 "false",
                                                 "--cells",
                                                 "1",
                                                 "--organism",
                                                 "dummy",
                                                 "--resources",
                                                 "2",
                                                 "--speed",
                                                 "1",
                                                 "--name",
                                                 "test",
                                                 "--timeout",
                                                 "1000000",
                                                 "--data-dir",
                                                 "../data/",
                                                 "--output",
                                                 "test_out_folder/",
                                                 "--threads",
                                                 "1",
                                                 "--seed",
                                                 "1"};

    std::vector<char *> argv_mock = unconst_string_vector(const_argv_mock);
    Configuration config(argv_mock.size(), argv_mock.data());

    cl_configuration_data arg_values = config.arguments();

    EXPECT_EQ(arg_values.seed, 1);

    srand(arg_values.seed);

    std::vector<std::pair<int, s_phase_checkpoints_t>> checkpoint_times;

    std::shared_ptr<DataManager> data = std::make_shared<DataManager>(
        arg_values.organism, arg_values.data_dir + "/database.sqlite",
        arg_values.data_dir + "/MFA-Seq_" + arg_values.organism + "/",
        arg_values.probability);

    unsigned long long seed = arg_values.seed;

    for (long long unsigned int i = 0; i < arg_values.cells; i++)
    {
        // Run all simulations with the same parameters, except for
        // seed, otherwise it would be exactly the same simulation
        // every time.
        PublicSPhase *s_phase = new PublicSPhase(
            arg_values.constitutive, arg_values.resources, arg_values.speed,
            arg_values.timeout, arg_values.period, arg_values.dormant, data,
            arg_values.organism, arg_values.name, arg_values.output, i ^ seed);

        s_phase->simulate(i);

        // check chromosome states

        std::vector<std::shared_ptr<Chromosome>> chrms =
            s_phase->get_genome()->chromosomes;

        std::ifstream chrm_0_times;
        chrm_0_times.open("../test/expected_outputs/dummy01_times.out");
        std::string replication_times(
            (std::istreambuf_iterator<char>(chrm_0_times)),
            std::istreambuf_iterator<char>());
        ASSERT_EQ(replication_times, chrms[0]->to_string());

        printf("metrics attach: %d\ndetach_normal: %d\ndetach_collision: %d\n",
               s_phase->get_fork_manager()->metric_times_attached,
               s_phase->get_fork_manager()->metric_times_detached_normal,
               s_phase->get_fork_manager()->metric_times_detached_collision);
        printf("iod: %d\n",
               s_phase->get_genome()->average_interorigin_distance());

        // check fork states

        std::cout << s_phase->get_fork_manager()->metric_times_attached
                  << std::endl;
        checkpoint_times.push_back(
            std::pair<int, s_phase_checkpoints_t>(i, s_phase->getTimes()));

        // check for interorigin distance states

        delete s_phase;
    }
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
