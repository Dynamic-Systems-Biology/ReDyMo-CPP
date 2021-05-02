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
  protected:
    // This allows us to test private methods and attributes.
    friend class EvolutionTest;

    unsigned long long seed;
    std::mt19937 rand_generator;
    Configuration &configuration;

    std::vector<std::vector<simulation_stats>> population;
    std::vector<std::shared_ptr<EvolutionDataProvider>> data_providers;

    cl_configuration_data arguments;

    int current_generation = 0;

    virtual void simulate();
    virtual void reproduce();

  public:
    EvolutionManager(Configuration &configuration, unsigned long long seed);
    ~EvolutionManager();

    void run_all();
    virtual void generation();

    virtual void snapshot(std::string folder);
};

typedef struct
{
    double collisions = 0;
    double time       = 0;
} instance_metrics;

double calculate_fitness(instance_metrics metrics,
                         cl_evolution_data config_data);

#endif
