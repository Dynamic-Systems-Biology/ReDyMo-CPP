#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "../include/chromosome.hpp"
#include "../include/data_manager.hpp"
#include "../include/util.hpp"

class DataManagerTest : public ::testing::Test
{
  protected:
    std::shared_ptr<DataManager> data;

  protected:
    DataManagerTest() {}

    void SetUp()
    {
        data = std::make_shared<DataManager>("../data/database.sqlite",
                                             "../data/MFA-Seq_dummy/");
    }

    void TearDown() {}
};

TEST_F(DataManagerTest, GenerateProbLandscape)
{
    std::vector<double> expected = {
        0.00010000000000021103, 0.3334000000000006, 1.0,
        0.6667000000000001,     0.6667000000000001, 0.3334000000000006,
        0.1667500000000004};
    std::vector<double> result = data->generate_prob_landscape("dummy_01", 7);
    for (int i = 0; i < (int)result.size(); i++)
        ASSERT_NEAR(result[i], expected[i], 1e10);
}

TEST_F(DataManagerTest, GetTranscriptionRegions)
{
    std::vector<transcription_region_t> result =
        data->get_transcription_regions("dummy_01");
    ASSERT_EQ(result[0].start, 25);
    ASSERT_EQ(result[0].end, 80);
}

TEST_F(DataManagerTest, GetConstitutiveOrigin)
{
    std::vector<constitutive_origin_t> result =
        data->get_constitutive_origins("dummy_01");
    ASSERT_EQ(result[0].base, 1234);
}

TEST_F(DataManagerTest, GetChromosomeData)
{
    std::vector<std::shared_ptr<Chromosome>> result =
        data->get_chromosome_data("dummy_01");
    std::vector<double> probabilities = {
        0.00010000000000021103, 0.3334000000000006, 1.0,
        0.6667000000000001,     0.6667000000000001, 0.3334000000000006,
        0.1667500000000004};
    for (auto chrm : result)
    {
        ASSERT_EQ(chrm->get_transcription_regions()[0].start, 25);
        ASSERT_EQ(chrm->get_transcription_regions()[0].end, 80);
        ASSERT_EQ(chrm->constitutive_origins[0].base, 1234);
        for (int i = 0; i < (int)probabilities.size(); i++)
            ASSERT_EQ(chrm->activation_probability(i), probabilities[i]);
    }
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}