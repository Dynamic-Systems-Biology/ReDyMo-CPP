#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "../include/chromosome.hpp"
#include "../include/genomic_location.hpp"
#include "../include/util.hpp"

class GenomicLocationTest : public ::testing::Test
{
  protected:
    GenomicLocation *gen_loc;

  protected:
    GenomicLocationTest() {}
    void SetUp()
    {
        Chromosome *chrm = create_chromosome();
        uint base        = rand() % chrm->size();
        gen_loc          = new GenomicLocation(base, *chrm);
    }

    Chromosome *create_chromosome(uint size = 300, uint id = 1)
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

        return new Chromosome(id, test_size, prob_landscape,
                              transcription_regions, cons_origins);
    }
    void TearDown()
    {
        delete &gen_loc->chromosome;
        delete gen_loc;
    }
};

TEST_F(GenomicLocationTest, OutOfRangeBase)
{
    ASSERT_THROW(GenomicLocation(302, gen_loc->chromosome),
                 std::invalid_argument);
    ASSERT_NO_THROW(GenomicLocation(10, gen_loc->chromosome));
}

TEST_F(GenomicLocationTest, IsReplicated)
{
    ASSERT_FALSE(gen_loc->is_replicated());
    gen_loc->chromosome.replicate(gen_loc->base, gen_loc->base, 2);
    ASSERT_TRUE(gen_loc->is_replicated());
}

TEST_F(GenomicLocationTest, WillActivate)
{
    float sum           = 0;
    Chromosome chrm     = *create_chromosome(1, 2);
    GenomicLocation loc = GenomicLocation(0, chrm);

    for (int i = 0; i < 1000; i++)
        sum += loc.will_activate(false, 1) ? 1 : 0;
    ASSERT_NEAR(sum, (float)1 / (1 + 1) * 1000, 6);
}

TEST_F(GenomicLocationTest, GetConstitutiveOrigin)
{
    ASSERT_TRUE(gen_loc->chromosome.constitutive_origins[0].base ==
                gen_loc->get_constitutive_origin(301)->base);
}

TEST_F(GenomicLocationTest, PutFiredConstitutiveOrigin)
{
    ASSERT_TRUE(gen_loc->chromosome.fired_constitutive_origins.empty());
    gen_loc->put_fired_constitutive_origin(
        gen_loc->chromosome.constitutive_origins[0]);
    ASSERT_FALSE(gen_loc->chromosome.fired_constitutive_origins.empty());
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}