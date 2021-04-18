#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "../include/chromosome.hpp"
#include "../include/genomic_location.hpp"
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
        origin.base = 70;
        cons_origins.resize(1, origin);
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

class GenomicLocationTest : public ::testing::Test
{
  protected:
    std::shared_ptr<GenomicLocation> gen_loc;
    unsigned long long seed = 0;

  protected:
    GenomicLocationTest() {}
    void SetUp()
    {
        std::shared_ptr<Chromosome> chrm = create_chromosome();
        uint base                        = rand() % chrm->size();
        gen_loc = std::make_shared<GenomicLocation>(base, chrm, seed);
    }

    std::shared_ptr<Chromosome> create_chromosome(uint size      = 300,
                                                  std::string id = "1")
    {
        std::shared_ptr<TestingProvider> provider(new TestingProvider(size));
        return std::make_shared<Chromosome>(id, provider);
    }
    void TearDown() {}
};

TEST_F(GenomicLocationTest, OutOfRangeBase)
{
    ASSERT_THROW(GenomicLocation(302, gen_loc->chromosome, seed),
                 std::invalid_argument);
    ASSERT_NO_THROW(GenomicLocation(10, gen_loc->chromosome, seed));
}

TEST_F(GenomicLocationTest, IsReplicated)
{
    ASSERT_FALSE(gen_loc->is_replicated());
    gen_loc->chromosome->replicate(gen_loc->base, gen_loc->base, 2);
    ASSERT_TRUE(gen_loc->is_replicated());
}

TEST_F(GenomicLocationTest, WillActivateNonConstitutive)
{
    double sum                       = 0;
    std::shared_ptr<Chromosome> chrm = create_chromosome(1, "2");
    GenomicLocation loc              = GenomicLocation(0, chrm, seed);

    for (int i = 0; i < 1000; i++)
        sum += loc.will_activate(false, 1) ? 1 : 0;
    ASSERT_NEAR(sum, ((double)1 / (1 + 1)) * 1000, 100);
}

TEST_F(GenomicLocationTest, WillActivateConstitutive)
{
    std::shared_ptr<Chromosome> chrm = create_chromosome(100, "3");
    GenomicLocation loc              = GenomicLocation(0, chrm, seed);

    ASSERT_FALSE(loc.will_activate(true, 138));
    ASSERT_TRUE(loc.will_activate(true, 140));
    ASSERT_FALSE(loc.will_activate(true, 139));
}

TEST_F(GenomicLocationTest, WillActivateConstitutiveAllFired)
{
    std::shared_ptr<Chromosome> chrm = create_chromosome(100, "3");
    GenomicLocation loc              = GenomicLocation(0, chrm, seed);
    ASSERT_TRUE(
        loc.put_fired_constitutive_origin((*chrm->constitutive_origins)[0]));
    ASSERT_FALSE(chrm->fired_constitutive_origins->empty());
    ASSERT_FALSE(loc.will_activate(true, 138));
    ASSERT_FALSE(loc.will_activate(true, 140));
    ASSERT_FALSE(loc.will_activate(true, 139));
}

TEST_F(GenomicLocationTest, GetConstitutiveOrigin)
{
    ASSERT_EQ((*gen_loc->chromosome->constitutive_origins)[0].base,
              gen_loc->get_constitutive_origin(600).base);
}

TEST_F(GenomicLocationTest, GetConstitutiveOriginAllFired)
{
    gen_loc->put_fired_constitutive_origin(
        (*gen_loc->chromosome->constitutive_origins)[0]);
    ASSERT_EQ(gen_loc->get_constitutive_origin(600).base, -1);
}

TEST_F(GenomicLocationTest, PutFiredConstitutiveOrigin)
{
    ASSERT_TRUE(gen_loc->chromosome->fired_constitutive_origins->empty());
    ASSERT_TRUE(gen_loc->put_fired_constitutive_origin(
        (*gen_loc->chromosome->constitutive_origins)[0]));
    ASSERT_FALSE(gen_loc->chromosome->fired_constitutive_origins->empty());
}

TEST_F(GenomicLocationTest, PutFiredConstitutiveOriginAlreadyFired)
{
    ASSERT_TRUE(gen_loc->chromosome->fired_constitutive_origins->empty());
    ASSERT_TRUE(gen_loc->put_fired_constitutive_origin(
        (*gen_loc->chromosome->constitutive_origins)[0]));
    ASSERT_FALSE(gen_loc->chromosome->fired_constitutive_origins->empty());
    ASSERT_FALSE(gen_loc->put_fired_constitutive_origin(
        (*gen_loc->chromosome->constitutive_origins)[0]));
}

TEST_F(GenomicLocationTest, PutFiredInvalidConstitutiveOrigin)
{
    constitutive_origin_t invalid_origin;
    invalid_origin.base = 71;
    ASSERT_TRUE(gen_loc->chromosome->fired_constitutive_origins->empty());
    ASSERT_FALSE(gen_loc->put_fired_constitutive_origin(invalid_origin));
    ASSERT_TRUE(gen_loc->chromosome->fired_constitutive_origins->empty());
}

TEST_F(GenomicLocationTest, SumOperator)
{
    std::shared_ptr<Chromosome> chrm = create_chromosome(100, "3");
    GenomicLocation loc(0, chrm, seed);

    EXPECT_EQ((loc + (-1)).base, 0);
    EXPECT_EQ((loc + 100).base, 99);
    loc += -1;
    EXPECT_EQ(loc.base, 0);
    loc += 100;
    ASSERT_EQ(loc.base, 99);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
