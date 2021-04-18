#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <omp.h>
#include <set>
#include <stdexcept>

#include "evolution.hpp"

double calculate_fitness(instance_metrics metrics,
                         cl_evolution_data config_data)
{
    double max_coll_all = 1 - 1 / metrics.collisions;
    double min_coll_all = 1 / metrics.collisions;

    double result = max_coll_all * config_data.fitness.max_coll_all +
                    min_coll_all * config_data.fitness.min_coll_all;

    return result;
}

EvolutionManager::EvolutionManager(Configuration &configuration,
                                   unsigned long long seed)
    : configuration(configuration), seed(seed)
{
    rand_generator.seed(seed);

    arguments = configuration.arguments();

    // Read Configuration
    for (int i = 0; i < arguments.evolution.population; i++)
        data_providers.push_back(std::make_shared<EvolutionDataProvider>(
            arguments.organism, arguments.data_dir + "/database.sqlite",
            arguments.data_dir + "/MFA-Seq_" + arguments.organism + "/", seed,
            arguments.probability));

    // Create first generation
    for (int i = 0; i < arguments.evolution.population; i++)
    {
        std::vector<simulation_stats> stats;
        for (int j = 0; j < arguments.cells; j++)
            stats.push_back(simulation_stats());
        population.push_back(stats);
    }
}

EvolutionManager::~EvolutionManager()
{
    std::cout << "[LIFE] Evolution Manager deleted!" << std::endl << std::flush;
}

void EvolutionManager::generation()
{
    current_generation++;

    std::cout << "[INFO] "
              << "Simulating generation " << current_generation << std::endl
              << std::flush;
    simulate();

    std::cout << "[INFO] "
              << "Reproducing generation " << current_generation << std::endl
              << std::flush;
    reproduce();
}

void EvolutionManager::reproduce()
{
    // Prepare metrics
    std::vector<instance_metrics> population_metrics;
    for (int i = 0; i < arguments.evolution.population; i++)
    {
        instance_metrics metrics;

        for (int c = 0; c < arguments.cells; c++)
        {
            simulation_stats stats = population[i][c];

            metrics.collisions += stats.collisions;
            metrics.time += metrics.time;
        }

        metrics.collisions /= arguments.cells;
        metrics.time /= arguments.cells;

        population_metrics.push_back(metrics);
    }

    // Calculate fitness
    std::vector<double> fitness;
    std::vector<double> inv_fitness;

    double max_fitness = -INFINITY;

    for (int i = 0; i < arguments.evolution.population; i++)
    {
        auto organism_fitness =
            calculate_fitness(population_metrics[i], arguments.evolution);
        fitness.push_back(organism_fitness);

        if (max_fitness < organism_fitness) max_fitness = organism_fitness;
    }

    for (int i = 0; i < arguments.evolution.population; i++)
        inv_fitness.push_back(max_fitness -
                              fitness[i]); // Best organism is never killed

    std::discrete_distribution<int> reproduction_roulette(fitness.begin(),
                                                          fitness.end());
    std::discrete_distribution<int> killing_roulette(inv_fitness.begin(),
                                                     inv_fitness.end());

    // Kill least successful
    int to_kill =
        arguments.evolution.population - arguments.evolution.survivors;

    int killed = 0;

    while (killed < to_kill)
    {
        auto kill_index = killing_roulette(rand_generator);
        if (!data_providers[kill_index]
                 ->isdead()) // Check if organism is not dead
        {
            data_providers[kill_index]->die();
            killed++;
        }
    }

    // Replicate most successful

    std::set<int> dead_organisms;

    for (int i = 0; i < arguments.evolution.population; i++)
    {
        if (data_providers[i]->isdead()) dead_organisms.insert(i);
    }

    for (auto organism_ptr = dead_organisms.begin();
         organism_ptr != dead_organisms.end(); organism_ptr++)
    {
        auto organism = *organism_ptr;

        do
        {
            auto index = reproduction_roulette(rand_generator);

            data_providers[organism]->clone(*data_providers[index]);
        } while (data_providers[organism]->isdead());
    }

    // Mutate
    for (auto provider = data_providers.begin();
         provider != data_providers.end(); provider++)
    {
        (*provider)->mutate(arguments);
    }
}

void EvolutionManager::simulate()
{
#pragma omp parallel for
    for (int i = 0; i < arguments.cells * arguments.evolution.population; i++)
    {
        int instance = i % arguments.evolution.population;
        int cell     = i / arguments.cells;

        SPhase s_phase(configuration, data_providers[cell],
                       i ^ seed + current_generation);
        s_phase.simulate(i);

        population[cell][instance] = s_phase.get_stats();
    }

    std::cout << "[INFO] Finished simulating" << std::endl << std::flush;
}

void EvolutionManager::run_all()
{
    for (int i = 0; i < arguments.evolution.generations; i++)
        generation();

    snapshot(arguments.output + "/final_snapshot");

    std::cout << "[INFO] Finished evolution simulation" << std::endl
              << std::flush;
}

void EvolutionManager::snapshot(std::string folder)
{
    for (int i = 0; i < data_providers.size(); i++)
    {
        data_providers[i]->snapshot(folder + std::string("/snapshot-genome-") +
                                    std::to_string(i));
    }
}
