#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "../include/chromosome.hpp"
#include "../include/genome.hpp"
#include "../include/util.hpp"

// class MockChromosome
// {
//   public:
//     uint n_fired_origins;

//   public:
//     MOCK_METHOD0(size, uint());
//     MOCK_METHOD0(is_replicated, bool());
//     MOCK_METHOD1(base_is_replicated, bool(int bas));
//     MOCK_METHOD0(n_constitutive_origins, uint());
// };

class GenomeTest : public ::testing::Test
{
  protected:
    Genome *gen;

  protected:
    GenomeTest() {}

    void SetUp()
    {
        std::vector<Chromosome *> chrms;
        for (int i = 0; i < 200; i++)
        {
            chrms.push_back(create_chromosome(300, std::to_string(i)));
        }
        gen = new Genome(chrms);
    }

    Chromosome *create_chromosome(uint size = 300, std::string id = "1")
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
        origin.base = 0;
        cons_origins.resize(test_size / 3, origin);

        return new Chromosome(id, test_size, prob_landscape,
                              transcription_regions, cons_origins);
    }

    void TearDown()
    {
        for (int i = 0; i < 200; i++)
        {
            delete gen->chromosomes[i];
        }
        delete gen;
    }
};

TEST_F(GenomeTest, Size) { ASSERT_EQ(200 * 300, gen->size()); }

TEST_F(GenomeTest, RandomGenomicLocation)
{
    for (int i = 0; i < 50; i++)
    {
        GenomicLocation loc = *gen->random_genomic_location();
        bool found          = false;
        for (auto chrm : gen->chromosomes)
        {
            if (*chrm == loc.chromosome) found = true;
        }
        ASSERT_TRUE(found);
    }
}

TEST_F(GenomeTest, RandomUnreplicatedGenLoc)
{

    for (int i = 0; i < 50; i++)
    {
        GenomicLocation loc = *gen->random_genomic_location();
        bool found          = false;
        for (auto chrm : gen->chromosomes)
        {
            if (*chrm == loc.chromosome) found = true;
        }
        ASSERT_TRUE(found);
        ASSERT_FALSE(loc.chromosome.is_replicated());
    }
}

TEST_F(GenomeTest, IsReplicated)
{
    ASSERT_FALSE(gen->is_replicated());
    for (auto chrm : gen->chromosomes)
    {
        chrm->replicate(0, 300, 1);
        ASSERT_TRUE(chrm->is_replicated());
    }
    ASSERT_TRUE(gen->is_replicated());
}

TEST_F(GenomeTest, AverageInterOriginDistance)
{
    ASSERT_EQ(300, gen->average_interorigin_distance());
    for (auto chrm : gen->chromosomes)
        chrm->add_fired_origin();
    ASSERT_EQ(150, gen->average_interorigin_distance());
    for (auto chrm : gen->chromosomes)
        chrm->add_fired_origin();
    ASSERT_EQ(100, gen->average_interorigin_distance());
    for (auto chrm : gen->chromosomes)
        chrm->add_fired_origin();
    ASSERT_EQ(75, gen->average_interorigin_distance());
}

TEST_F(GenomeTest, NConstitutiveOrigins)
{
    ASSERT_EQ(200 * 300 / 3, gen->n_constitutive_origins());
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}