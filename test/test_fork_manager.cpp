#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "../include/chromosome.hpp"
#include "../include/fork_manager.hpp"
#include "../include/util.hpp"

class TestingProvider : public DataProvider
{
  private:
    int size;
    std::vector<double> prob_landscape;
    std::vector<transcription_region_t> transcription_regions;
    std::vector<constitutive_origin_t> cons_origins;

  public:
    TestingProvider(uint size) : size(size)
    {
        prob_landscape.resize(size, (double)1 / (size + 1));

        transcription_region_t reg;
        reg.start = 1000;
        reg.end   = 2600;

        transcription_region_t reg2;
        reg2.start = 260;
        reg2.end   = 100;

        transcription_regions.push_back(reg);
        transcription_regions.push_back(reg2);

        constitutive_origin_t origin;
        origin.base = 0;
        cons_origins.push_back(origin);
    }

    const std::vector<std::string> &get_codes()
    {
        std::vector<std::string> codes;
        return codes;
    }

    int get_length(std::string code) { return size; }

    const std::vector<double> &get_probability_landscape(std::string code)
    {
        return prob_landscape;
    }

    const std::shared_ptr<std::vector<transcription_region_t>>
    get_transcription_regions(std::string code)
    {
        return std::make_shared<std::vector<transcription_region_t>>(
            transcription_regions);
    }

    const std::shared_ptr<std::vector<constitutive_origin_t>>
    get_constitutive_origins(std::string code)
    {
        return std::make_shared<std::vector<constitutive_origin_t>>(
            cons_origins);
    }
};

class ForkManagerTest : public ::testing::Test
{
  protected:
    std::shared_ptr<ForkManager> manager;
    std::shared_ptr<Genome> gen;

  protected:
    ForkManagerTest() {}

    void SetUp()
    {
        std::vector<std::shared_ptr<Chromosome>> chrms(1);
        chrms[0] = create_chromosome(3000, "2");
        gen      = std::make_shared<Genome>(chrms);
        manager  = std::make_shared<ForkManager>(3, gen, 15);
    }

    void TearDown() {}

    std::shared_ptr<Chromosome> create_chromosome(uint size      = 3000,
                                                  std::string id = "1")
    {
        std::shared_ptr<TestingProvider> provider(new TestingProvider(size));

        return std::make_shared<Chromosome>(id, provider);
    }
};

TEST_F(ForkManagerTest, CheckConflicts)
{
    GenomicLocation loc(1400, gen->chromosomes[0]);
    manager->attach_forks(loc, 10);
    ASSERT_EQ(
        manager->check_replication_transcription_conflicts(1400, 1000, true),
        1);
    ASSERT_TRUE(manager->replication_forks[0]->is_attached());
    ASSERT_FALSE(manager->replication_forks[1]->is_attached());
}

TEST_F(ForkManagerTest, CheckConflictsReversed)
{
    GenomicLocation loc(220, gen->chromosomes[0]);
    manager->attach_forks(loc, 10);
    ASSERT_EQ(
        manager->check_replication_transcription_conflicts(140, 100, true),
        1);
    ASSERT_FALSE(manager->replication_forks[0]->is_attached());
    ASSERT_TRUE(manager->replication_forks[1]->is_attached());
}

TEST_F(ForkManagerTest, AdvanceAttachedForks)
{
    GenomicLocation loc(1800, gen->chromosomes[0]);
    manager->attach_forks(loc, 10);
    manager->replication_forks[2]->set_just_detached(true);
    manager->advance_attached_forks(12);
    ASSERT_FALSE(manager->replication_forks[2]->get_just_detached());
    ASSERT_TRUE(manager->replication_forks[0]->is_attached());
    ASSERT_TRUE(manager->replication_forks[1]->is_attached());
    ASSERT_EQ(manager->replication_forks[0]->get_base(), 1815);
    ASSERT_EQ(manager->replication_forks[1]->get_base(), 1785);
}

TEST_F(ForkManagerTest, AttachForks)
{
    GenomicLocation loc(1800, gen->chromosomes[0]);
    manager->attach_forks(loc, 10);
    ASSERT_TRUE(manager->replication_forks[0]->is_attached());
    ASSERT_TRUE(manager->replication_forks[1]->is_attached());
    ASSERT_FALSE(manager->replication_forks[2]->is_attached());
    ASSERT_EQ(manager->replication_forks[0]->get_base(), 1800);
    ASSERT_EQ(manager->replication_forks[1]->get_base(), 1800);
    ASSERT_EQ(gen->chromosomes[0]->get_n_replicated_bases(), 1);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
