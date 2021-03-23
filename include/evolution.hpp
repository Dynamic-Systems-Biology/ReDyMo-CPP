/*! File Evolution.hpp
 *  Contains the Evolution class.
 */
#ifndef __EVOLUTION_HPP__
#define __EVOLUTION_HPP__

#include "configuration.hpp"
#include "evolution_data_provider.hpp"
#include "s_phase.hpp"
#include "util.hpp"
#include <memory>
#include <string>
#include <vector>

class EvolutionManager
{
  private:
    int seed;
    std::mt19937 rand_generator;
    Configuration &configuration;

    std::vector<std::vector<simulation_stats>> population;
    std::vector<std::shared_ptr<EvolutionDataProvider>> data_providers;

    cl_configuration_data arguments;

    int current_generation = 0;

    void simulate();
    void reproduce();

  public:
    EvolutionManager(Configuration &configuration, int seed);
    ~EvolutionManager();

    void run_all();
    void generation();

    void snapshot(std::string folder);
};

#endif