#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "../include/genome.hpp"
#include "../include/genomic_location.hpp"
#include "../include/replication_fork.hpp"
#include "../include/util.hpp"

class ReplicationForkTest : public ::testing::Test
{

  protected:
    std::shared_ptr<ReplicationFork> fork;
    std::vector<std::shared_ptr<Chromosome>> chrms;

  protected:
    ReplicationForkTest() {}

    void SetUp()
    {
        for (int i = 0; i < 100; i++)
            chrms.push_back(create_chromosome(300, std::to_string(i)));
        std::shared_ptr<Genome> gen = std::make_shared<Genome>(chrms);
        fork = std::make_shared<ReplicationFork>(gen, 40);
    }

    void TearDown() {}

    std::shared_ptr<Chromosome> create_chromosome(uint size      = 300,
                                                  std::string id = "1")
    {
        uint test_size = size;
        std::vector<float> prob_landscape;
        std::vector<transcription_region_t> transcription_regions;
        std::vector<constitutive_origin_t> cons_origins;

        prob_landscape.resize(test_size, (float)1 / (test_size + 1));

        transcription_region_t reg;
        reg.start = 0;
        reg.end   = 10;

        transcription_regions.resize(test_size / 4, reg);

        constitutive_origin_t origin;
        origin.base = 70;
        cons_origins.resize(1, origin);

        return std::make_shared<Chromosome>(
            id, test_size, prob_landscape, transcription_regions, cons_origins);
    }
};

TEST_F(ReplicationForkTest, AlreadyAttached)
{
    GenomicLocation loc(2, chrms[1]);
    fork->attach(loc, 1, 2);
    ASSERT_ANY_THROW(fork->attach(loc, 1, 2));
}

TEST_F(ReplicationForkTest, AttachAndGetters)
{
    GenomicLocation loc(2, chrms[1]);
    fork->attach(loc, 1, 2);
    ASSERT_EQ(fork->get_base(), 2);
    ASSERT_EQ(fork->get_direction(), 1);
    ASSERT_TRUE(*fork->get_chromosome() == *chrms[1]);
    ASSERT_FALSE(fork->get_just_detached());
}

TEST_F(ReplicationForkTest, Detach)
{
    GenomicLocation loc(2, chrms[1]);
    fork->attach(loc, 1, 2);
    fork->detach();
    ASSERT_EQ(-1, fork->get_base());
    ASSERT_EQ(0, fork->get_direction());
    ASSERT_EQ(nullptr, fork->get_chromosome());
}

TEST_F(ReplicationForkTest, Advance)
{
    GenomicLocation loc(2, chrms[1]);
    fork->attach(loc, 1, 2);
    ASSERT_TRUE(fork->advance(3));

    for (int i = 0; i < 40; i++)
        ASSERT_TRUE(chrms[1]->base_is_replicated(2 + i));
}

TEST_F(ReplicationForkTest, IsAttached)
{
    ASSERT_FALSE(fork->is_attached());
    GenomicLocation loc(2, chrms[1]);
    fork->attach(loc, 1, 2);
    ASSERT_TRUE(fork->is_attached());
    fork->detach();
    ASSERT_FALSE(fork->is_attached());
}

TEST_F(ReplicationForkTest, JustDetached)
{
    ASSERT_FALSE(fork->get_just_detached());

    GenomicLocation loc(2, chrms[1]);
    fork->attach(loc, 1, 2);

    ASSERT_FALSE(fork->get_just_detached());

    fork->detach();

    ASSERT_FALSE(fork->get_just_detached());
    GenomicLocation loc2(298, chrms[1]);
    fork->attach(loc2, 1, 4);
    fork->advance(5);

    ASSERT_TRUE(fork->get_just_detached());
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}