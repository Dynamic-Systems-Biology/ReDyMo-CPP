#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "../include/chromosome.hpp"
#include "../include/genome.hpp"
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
        reg.start = 0;
        reg.end   = 10;

        transcription_regions.resize(size / 4, reg);

        constitutive_origin_t origin;
        origin.base = 0;
        cons_origins.resize(size / 3, origin);
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

class GenomeTest : public ::testing::Test
{
  protected:
    std::shared_ptr<Genome> gen;

  protected:
    GenomeTest() {}

    void SetUp()
    {
        std::vector<std::shared_ptr<Chromosome>> chrms;
        for (int i = 0; i < 200; i++)
        {
            chrms.push_back(create_chromosome(300, std::to_string(i)));
        }
        gen = std::make_shared<Genome>(chrms);
    }

    std::shared_ptr<Chromosome> create_chromosome(uint size      = 300,
                                                  std::string id = "1")
    {
        std::shared_ptr<TestingProvider> provider(new TestingProvider(size));

        return std::make_shared<Chromosome>(id, provider);
    }

    void TearDown() {}
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
            if (chrm == loc.chromosome) found = true;
        }
        ASSERT_TRUE(found);
    }
}

TEST_F(GenomeTest, RandomUnreplicatedGenomicLocation)
{
    std::vector<std::shared_ptr<Chromosome>> chrms;
    chrms.push_back(create_chromosome(300, "1"));
    gen = std::make_shared<Genome>(chrms);

    gen->chromosomes[0]->replicate(1, 299, 1);

    GenomicLocation loc = *gen->random_unreplicated_genomic_location();
    bool found          = false;
    for (auto chrm : gen->chromosomes)
    {
        if (chrm == loc.chromosome) found = true;
    }
    ASSERT_EQ(loc.base, 0);
    ASSERT_TRUE(found);

    gen->chromosomes[0]->replicate(0, 0, 2);
    ASSERT_THROW(gen->random_unreplicated_genomic_location(), std::runtime_error);
}

TEST_F(GenomeTest, RandomUnreplicatedGenLoc)
{

    for (int i = 0; i < 50; i++)
    {
        GenomicLocation loc = *gen->random_genomic_location();
        bool found          = false;
        for (auto chrm : gen->chromosomes)
        {
            if (chrm == loc.chromosome) found = true;
        }
        ASSERT_TRUE(found);
        ASSERT_FALSE(loc.chromosome->is_replicated());
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
