#include "data_manager.hpp"
#include "chromosome.hpp"
#include <SQLiteCpp/SQLiteCpp.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#undef __UNIFORM_LANDSCAPE__

DataManager::DataManager(std::string organism, std::string database_path,
                         std::string mfa_seq_data_path, double p)
    : database_path(database_path), mfa_seq_data_path(mfa_seq_data_path),
      uniform(p)
{
    try {
        printf("[INFO] Loading database data for organism \"%s\" \n",
               organism.c_str());
        SQLite::Database db(database_path, SQLite::OPEN_READONLY);
        SQLite::Statement query(db, "select * from Chromosome where organism = ?");
        query.bind(1, organism);

        while (query.executeStep())
        {
            std::string code = query.getColumn("code").getString();
            printf("[INFO] Loading chromosome %s\n", code.c_str());
            int length       = query.getColumn("length").getInt();

            lengths.insert(std::pair<std::string, int>(code, length));

            codes.push_back(code);

            try
            {
                generate_prob_landscape(code, length);
                generate_constitutive_origins(code);
                generate_transcription_regions(code);
            }
            catch (std::out_of_range &e)
            {
                std::cerr << e.what() << std::endl;
                exit(-1);
            }
        }
    }
    catch (std::exception &e) {
        std::cout << "An error occurred while loading database data: " << e.what()
                  << std::endl
                  << std::fflush(nullptr);
        exit(1);
    }

    if (codes.empty())
    {
        std::cout << "[ERROR] No chromosomes found for organism " << organism << std::endl
                  << std::fflush(nullptr);
        exit(1);
    }

}

DataManager::~DataManager()
{
    std::cout << "[LIFE] Data Manager deleted!" << std::endl << std::flush;
}

const std::vector<std::string> &DataManager::get_codes() { return codes; }

void DataManager::generate_transcription_regions(std::string code)
{
    try
    {
        SQLite::Database db(database_path, SQLite::OPEN_READONLY);
        SQLite::Statement query(
            db, "SELECT * FROM TranscriptionRegion WHERE chromosome_code = ?");
        query.bind(1, code);

        transcription_regions.insert(
            std::pair<std::string, std::vector<transcription_region_t>>(
                code, std::vector<transcription_region_t>()));

        auto &regions = transcription_regions.at(code);

        while (query.executeStep())
        {
            transcription_region_t region;
            region.start = query.getColumn("start").getInt();
            region.end   = query.getColumn("end").getInt();

            regions.push_back(region);
        }
    }
    catch (int e)
    {
        std::cout << "An error ocurred while loading database data. Error " << e
                  << std::endl
                  << std::fflush(nullptr);
        exit(1);
    }
}

void DataManager::generate_constitutive_origins(std::string code)
{
    try
    {
        SQLite::Database db(database_path, SQLite::OPEN_READONLY);
        SQLite::Statement query(
            db, "SELECT * FROM ReplicationOrigin WHERE chromosome_code = ?");
        query.bind(1, code);

        constitutive_origins.insert(
            std::pair<std::string, std::vector<constitutive_origin_t>>(
                code, std::vector<constitutive_origin_t>()));

        auto &origins = constitutive_origins.at(code);

        while (query.executeStep())
        {
            constitutive_origin_t origin;
            origin.base = query.getColumn("position").getInt();
            origins.push_back(origin);
        }
    }
    catch (int e)
    {
        std::cout << "An error ocurred while loading database data. Error " << e
                  << std::endl
                  << std::fflush(nullptr);
        exit(1);
    }
}

void DataManager::generate_prob_landscape(std::string code, uint length)
{
    std::ifstream mfa_file;
    std::vector<double> scores;
    probability_landscape.insert(std::pair<std::string, std::vector<double>>(
        code, std::vector<double>(length, 0.0)));

    auto &landscape = probability_landscape.at(code);
    double curr_score;

    try
    {
        mfa_file.open(mfa_seq_data_path + code + ".txt");
        if (!mfa_file.is_open()) { throw 1; }
    }
    catch (int e)
    {
        std::cout << "An error occurred while loading MFA_Seq["
                  << mfa_seq_data_path + code + ".txt"
                  << "] data. Error " << e << std::endl
                  << std::fflush(nullptr);
        exit(1);
    }

    while (!mfa_file.eof())
    {
        mfa_file >> curr_score;
        if (curr_score != EOF) scores.push_back(curr_score);
        curr_score = EOF;
    }
    mfa_file.close();

    int step = (int)(std::ceil(length / (float)(scores.size())));

    double a =
        (1 - pow(10, -4)) / (*std::max_element(scores.begin(), scores.end()) -
                             *std::min_element(scores.begin(), scores.end()));

    double b = 1 - (*std::max_element(scores.begin(), scores.end()) * a);

    #ifdef __UNIFORM_LANDSCAPE__
        double sum  = 0;
        double mean = 0;
    #endif

    std::ofstream probs_file;
    probs_file.open(mfa_seq_data_path + code + "_probability.txt");
    for (int i = 0; i < (int)scores.size(); i++)
    {
        double prob = a * scores[i] + b;

        if (uniform) prob = uniform;

        probs_file << std::fixed << std::setprecision(17) << prob << std::endl;
        for (int j = i * step; j < (i + 1) * step; j++)
        {
            landscape[j] = prob;
            #ifdef __UNIFORM_LANDSCAPE__
                sum += prob;
            #endif
            if (j == (int)length - 1)
            {
                probs_file.close();
                #ifdef __UNIFORM_LANDSCAPE__
                    mean = sum / landscape.size();
                    for (int k = 0; k < landscape.size(); k++)
                    {
                        landscape[k] = mean;
                    }
                #endif
                return;
            }
        }
    }
}

int DataManager::get_length(std::string code)
{
    std::lock_guard<std::mutex> guard(lengths_mutex);
    try
    {
        return lengths.at(code);
    }
    catch (std::out_of_range &e)
    {
        std::cerr << code << ": " << e.what() << std::endl;
        exit(-1);
    }
}

const std::vector<double> &
DataManager::get_probability_landscape(std::string code)
{
    std::lock_guard<std::mutex> guard(prob_landscape_mutex);

    try
    {
        return probability_landscape.at(code);
    }
    catch (std::out_of_range &e)
    {
        auto codes = get_codes();
        for (auto code = codes.begin(); code != codes.end(); code++)
            std::cerr << *code << std::endl << std::flush;

        std::cerr << code << "[P]: " << e.what()
                  << " - size: " << probability_landscape.size() << std::endl;
        exit(-1);
    }
}

const std::shared_ptr<std::vector<transcription_region_t>>
DataManager::get_transcription_regions(std::string code)
{
    std::lock_guard<std::mutex> guard(transcription_regions_mutex);
    try
    {
        return std::make_shared<std::vector<transcription_region_t>>(
            transcription_regions.at(code));
    }
    catch (std::out_of_range &e)
    {
        auto codes = get_codes();
        for (auto code = codes.begin(); code != codes.end(); code++)
            std::cerr << *code << std::endl << std::flush;

        std::cerr << code << "[T]: " << e.what() << std::endl;
        exit(-1);
    }
}

const std::shared_ptr<std::vector<constitutive_origin_t>>
DataManager::get_constitutive_origins(std::string code)
{
    std::lock_guard<std::mutex> guard(constitutive_origins_mutex);
    try
    {
        return std::make_shared<std::vector<constitutive_origin_t>>(
            constitutive_origins.at(code));
    }
    catch (std::out_of_range &e)
    {

        auto codes = get_codes();
        for (auto code = codes.begin(); code != codes.end(); code++)
            std::cerr << *code << std::endl << std::flush;

        std::cerr << code << "[C]: " << test << std::endl;
        exit(-1);
    }
}
