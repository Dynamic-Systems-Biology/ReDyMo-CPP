#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <iostream>
#include <memory>

#include "../include/chromosome.hpp"

class ChromosomeTest : public ::testing::Test
{
  protected:
    Chromosome chrm;

  protected:
    ChromosomeTest() {}

    void SetUp() { chrm = *create_chromosome(); }

    std::shared_ptr<Chromosome> create_chromosome(uint size = 300)
    {
        uint test_size = size;
        std::vector<float> prob_landscape;
        std::vector<transcription_region_t> transcription_regions;
        std::vector<constitutive_origin_t> cons_origins;

        prob_landscape.resize(test_size, (float)1 / (test_size + 1));

        transcription_region_t reg;
        reg.start = 0;
        reg.end   = 10;

        transcription_regions.resize(test_size, reg);

        constitutive_origin_t origin;
        origin.base = 0;
        cons_origins.resize(test_size, origin);

        return std::shared_ptr<Chromosome>(
            new Chromosome("1", test_size, prob_landscape,
                           transcription_regions, cons_origins));
    }

    void TearDown() {}
};

/*! Tests if throws exception on invalid length.
 */
TEST_F(ChromosomeTest, ZeroOrNegativeLength)
{
    ASSERT_THROW(create_chromosome(0), std::invalid_argument);
}

/*! Tests if length returns the right value.
 */
TEST_F(ChromosomeTest, Size)
{
    auto chrm_300 = create_chromosome(300);
    ASSERT_EQ(300, chrm_300->size());
    auto chrm_500 = create_chromosome(500);
    ASSERT_EQ(500, chrm_500->size());
}

/*! Tests if an exception is thrown in case of non existant base.
 */
TEST_F(ChromosomeTest, OutOfRangeBaseIsReplicated)
{
    ASSERT_THROW(chrm.base_is_replicated(400), std::out_of_range);
    ASSERT_THROW(chrm.base_is_replicated(301), std::out_of_range);
}

/*! Tests if a newly created Chromosome has no replicated bases and throws no
 * exception on valid bases.
 */
TEST_F(ChromosomeTest, BaseIsReplicated)
{
    ASSERT_NO_THROW(chrm.base_is_replicated(0));
    ASSERT_NO_THROW(chrm.base_is_replicated(299));

    for (int i = 0; i < 300; i++)
        ASSERT_FALSE(chrm.base_is_replicated(i));

    ASSERT_NO_THROW(chrm.replicate(4, 6, 1));

    ASSERT_TRUE(chrm.base_is_replicated(4));
    ASSERT_TRUE(chrm.base_is_replicated(5));
}

/*! Tests if the method throws exception when the base is invalid.
 */
TEST_F(ChromosomeTest, OutOfRangeActivationProbability)
{
    ASSERT_THROW(chrm.activation_probability(301), std::out_of_range);
    ASSERT_THROW(chrm.activation_probability(400), std::out_of_range);
}

/*! Tests if the method does not throw exceptions when the base is valid.
 */
TEST_F(ChromosomeTest, InRangeActivationProbability)
{
    ASSERT_NO_THROW(chrm.activation_probability(0));
    ASSERT_NO_THROW(chrm.activation_probability(10));
}

/*! Tests if the method returns the right activation probability.
 */
TEST_F(ChromosomeTest, ActivationProbability)
{
    auto chrm_300 = create_chromosome(300);
    ASSERT_EQ((float)1 / 301, chrm_300->activation_probability(0));
    auto chrm_500 = create_chromosome(500);
    ASSERT_EQ((float)1 / 501, chrm_500->activation_probability(4));
}

/*! Tests if the method sets correctly(like a gaussian) the activation
 * probability.
 */
TEST_F(ChromosomeTest, SetDormantActivationProbability)
{
    int size       = 300000;
    auto chrm_300k = create_chromosome(size);
    chrm_300k->set_dormant_activation_probability(size / 2);
    for (int i = 130000; i < 170000; i++)
    {
        ASSERT_GT((float)1 / size + 1, chrm_300k->activation_probability(i));
    }
}

/*! Tests if the method throws exception when the start base is invalid.
 */
TEST_F(ChromosomeTest, OutOfRangeReplicate)
{
    ASSERT_THROW(chrm.replicate(-3, 6, 30), std::out_of_range);
    ASSERT_THROW(chrm.replicate(301, 150, 30), std::out_of_range);
}

/*! Tests if the method does not throw exceptions when the start base is valid.
 */
TEST_F(ChromosomeTest, InRangeReplicate)
{
    ASSERT_NO_THROW(chrm.replicate(3, 1500, 30));
    ASSERT_NO_THROW(chrm.replicate(3, 15, 30));
}

TEST_F(ChromosomeTest, Replicate)
{
    ASSERT_EQ(chrm.replicate(0, 30, 1), true);
    for (int base = 0; base < 30; base++)
        ASSERT_EQ(true, chrm.base_is_replicated(base));
    ASSERT_EQ(chrm.get_n_replicated_bases(), 31);
}

TEST_F(ChromosomeTest, IsReplicated)
{
    ASSERT_FALSE(chrm.is_replicated());
    chrm.replicate(0, 300, 1);
    ASSERT_TRUE(chrm.is_replicated());
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}