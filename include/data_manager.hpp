#ifndef __DATA_MANAGER_HPP__
#define __DATA_MANAGER_HPP__

#include "chromosome.hpp"
#include "data_provider.hpp"
#include "util.hpp"
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

/*! The DataManager class is a data provider that loads data from the sqlite
 * database and MFA-Seq files to provide information to the Genomes.
 */
class DataManager : public DataProvider
{
  private:
    long long unsigned test = 0;
    std::string database_path;
    std::string mfa_seq_data_path;
    double uniform;

    std::vector<std::string> codes;
    std::unordered_map<std::string, int> lengths;

    void generate_prob_landscape(std::string code, uint length);

    void generate_transcription_regions(std::string code);

    void generate_constitutive_origins(std::string code);

  protected:
    std::mutex lengths_mutex;
    std::mutex prob_landscape_mutex;
    std::mutex transcription_regions_mutex;
    std::mutex constitutive_origins_mutex;

    std::unordered_map<std::string, std::vector<double>> probability_landscape;
    std::unordered_map<std::string, std::vector<constitutive_origin_t>>
        constitutive_origins;
    std::unordered_map<std::string, std::vector<transcription_region_t>>
        transcription_regions;

  public:
    DataManager(std::string organism, std::string database_path,
                std::string mfa_seq_data_path, double p = 0);
    ~DataManager();

    const std::vector<std::string> &get_codes();
    int get_length(std::string code);
    const std::vector<double> &get_probability_landscape(std::string code);
    const std::shared_ptr<std::vector<transcription_region_t>>
    get_transcription_regions(std::string code);
    const std::shared_ptr<std::vector<constitutive_origin_t>>
    get_constitutive_origins(std::string code);
};

#endif
