#ifndef __DATA_MANAGER_HPP__
#define __DATA_MANAGER_HPP__

#include "chromosome.hpp"
#include "util.hpp"
#include <memory>
#include <string>
#include <vector>

/*! The DataManager class handles all input and loading of data.
 * It stores Chromosome data, activation probabilities of replication origins
 * and information about Transcription Regions.
 */
class DataManager
{
  private:
    std::string database_path;
    std::string mfa_seq_data_path;

  public:
    DataManager(std::string database_path, std::string mfa_seq_data_path);

    std::vector<std::shared_ptr<Chromosome>>
    get_chromosome_data(std::string organism);

    std::vector<float> generate_prob_landscape(uint code, uint length);

    std::vector<std::shared_ptr<Chromosome>>
    assemble_chromosomes(std::string organism);

    std::vector<transcription_region_t>
    get_transcription_regions(uint chromosome_code);

    std::vector<constitutive_origin_t>
    get_constitutive_origins(uint chromosome_code);
};

#endif