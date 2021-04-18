#ifndef __EVOLUTION_DATA_PROVIDER_HPP__
#define __EVOLUTION_DATA_PROVIDER_HPP__

#include "chromosome.hpp"
#include "configuration.hpp"
#include "data_manager.hpp"
#include "genome.hpp"
#include "genomic_location.hpp"
#include "util.hpp"
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

class GenomicLocation;
class Genome;

struct BellCurve
{
    struct
    {
        std::string chromosome;
        long base;
    } location;

    double sigma;
};

/*! The EvolutionDataProvider class is a data provider that loads data from the
 * sqlite database and MFA-Seq files to provide information to the Genomes, and
 * contains extra functionality for simulating mutations and replication.
 */
class EvolutionDataProvider : public DataManager
{
  private:
    std::mt19937 rand_generator;
    std::vector<BellCurve> bell_curves;
    std::unordered_map<std::string, std::vector<double>> original_landscape;
    bool dead = false;

  public:
    EvolutionDataProvider(std::string organism, std::string database_path,
                          std::string mfa_seq_data_path,
                          unsigned long long seed, double p = 0);

    void clone(EvolutionDataProvider &provider);
    void mutate(cl_configuration_data config);

    void snapshot(std::string folder);
    void die();
    bool isdead();
};

#endif
