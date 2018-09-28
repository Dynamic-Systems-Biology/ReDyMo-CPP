#include "data_manager.hpp"
#include "chromosome.hpp"
#include <SQLiteCpp/SQLiteCpp.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

DataManager::DataManager(std::string database_path,
                         std::string mfa_seq_data_path)
    : database_path(database_path), mfa_seq_data_path(mfa_seq_data_path)
{
}

std::vector<std::shared_ptr<Chromosome>>
DataManager::get_chromosome_data(std::string organism)
{
    std::vector<std::shared_ptr<Chromosome>> chrms;
    SQLite::Database db(database_path, SQLite::OPEN_READONLY);
    SQLite::Statement query(db, "select * from Chromosome where organism = ?");
    query.bind(1, organism);

    std::string code;
    int length;

    while (query.executeStep())
    {
        code   = query.getColumn("code").getString();
        length = query.getColumn("length").getInt();

        std::vector<float> prob = generate_prob_landscape(code, length);
        std::vector<transcription_region_t> regions =
            get_transcription_regions(code);
        std::vector<constitutive_origin_t> origins =
            get_constitutive_origins(code);

        std::shared_ptr<Chromosome> chrm =
            std::make_shared<Chromosome>(code, length, prob, regions, origins);
        chrms.push_back(chrm);
    }
    return chrms;
}

std::vector<transcription_region_t>
DataManager::get_transcription_regions(std::string chromosome_code)
{
    SQLite::Database db(database_path, SQLite::OPEN_READONLY);
    SQLite::Statement query(
        db, "SELECT * FROM TranscriptionRegion WHERE chromosome_code = ?");
    query.bind(1, chromosome_code);

    std::vector<transcription_region_t> regions;
    while (query.executeStep())
    {
        transcription_region_t region;
        region.start = query.getColumn("start").getInt();
        region.end   = query.getColumn("end").getInt();
        regions.push_back(region);
    }
    return regions;
}

std::vector<constitutive_origin_t>
DataManager::get_constitutive_origins(std::string chromosome_code)
{
    SQLite::Database db(database_path, SQLite::OPEN_READONLY);
    SQLite::Statement query(
        db, "SELECT * FROM ReplicationOrigin WHERE chromosome_code = ?");
    query.bind(1, chromosome_code);

    std::vector<constitutive_origin_t> origins;
    while (query.executeStep())
    {
        constitutive_origin_t origin;
        origin.base = query.getColumn("position").getInt();
        origins.push_back(origin);
    }
    return origins;
}

std::vector<float> DataManager::generate_prob_landscape(std::string code,
                                                        uint length)
{
    std::ifstream mfa_file;
    std::vector<float> scores;
    std::vector<float> probabilities(length);
    float curr_score;

    mfa_file.open(mfa_seq_data_path + code + ".txt");

    while (!mfa_file.eof())
    {
        mfa_file >> curr_score;
        scores.push_back(curr_score);
    }
    mfa_file.close();

    int step = (int)(std::ceil(length / scores.size()));

    float a =
        (1 - pow(10, -4)) / (*std::max_element(scores.begin(), scores.end()) -
                             *std::min_element(scores.begin(), scores.end()));

    float b = 1 - (*std::max_element(scores.begin(), scores.end()) * a);

    std::ofstream probs_file;
    probs_file.open(mfa_seq_data_path + code + "_probability.txt");
    for (int i = 0; i < (int)scores.size(); i++)
    {
        float prob = a * scores[i] + b;
        probs_file << prob;
        for (int j = i * step; j < (i + 1) * step; j++)
        {
            probabilities[j] = prob;
            if (j == (int)length - 1) return probabilities;
        }
    }
    probs_file.close();
    return probabilities;
}
