#include "evolution_data_provider.hpp"
#include "chromosome.hpp"
#include <SQLiteCpp/SQLiteCpp.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <random>
#include <string>
#include <vector>

EvolutionDataProvider::EvolutionDataProvider(std::string organism,
                                             std::string database_path,
                                             std::string mfa_seq_data_path,
                                             double p)
    : DataManager(organism, database_path, mfa_seq_data_path, p)
{
    rand_generator.seed((long long)this);
    // std::cout << "[DEBUG] Evolution Data Provider created!" << std::endl <<
    // std::flush;

    original_landscape = probability_landscape;
}

void EvolutionDataProvider::clone(EvolutionDataProvider &provider)
{
    dead                  = false;
    probability_landscape = provider.probability_landscape;
    constitutive_origins  = provider.constitutive_origins;
    transcription_regions = provider.transcription_regions;
}

void EvolutionDataProvider::mutate(cl_configuration_data config)
{
    std::lock_guard<std::mutex> prob_guard(prob_landscape_mutex);
    std::lock_guard<std::mutex> transcription_guard(
        transcription_regions_mutex);
    std::lock_guard<std::mutex> origins_guard(constitutive_origins_mutex);

    auto codes = get_codes();

    std::uniform_real_distribution<double> uniform1(0, 1);

    // Single Bell Curve Mutations
    for (auto bellptr = bell_curves.begin(); bellptr != bell_curves.end();
         bellptr++)
    {
        auto &curve = *bellptr;

        if (config.evolution.mutations.probability_landscape.change_mean.prob >
            uniform1(rand_generator))
        {
            std::normal_distribution<double> normal(
                0, config.evolution.mutations.probability_landscape.change_mean
                       .std);

            curve.location.base += normal(rand_generator);
            auto length = get_length(curve.location.chromosome);

            if (curve.location.base < 0)
                curve.location.base = 0;
            else if (curve.location.base >= length)
                curve.location.base = length;
        }

        if (config.evolution.mutations.probability_landscape.change_std.prob >
            uniform1(rand_generator))
        {
            std::normal_distribution<double> normal(
                0, config.evolution.mutations.probability_landscape.change_std
                       .std);

            curve.sigma += normal(rand_generator);
        }
    }

    // Bell Curve Global Mutations
    for (auto code = codes.begin(); code != codes.end(); code++)
    {

        if (config.evolution.mutations.probability_landscape.add >
            uniform1(rand_generator))
        {
            auto curve  = BellCurve();
            auto length = get_length(*code);
            std::uniform_int_distribution<> distribution(length);

            curve.sigma               = 1;
            curve.location.chromosome = *code;
            curve.location.base       = distribution(rand_generator);

            bell_curves.push_back(curve);
        }

        if (config.evolution.mutations.probability_landscape.del >
            uniform1(rand_generator))
        {
            std::uniform_int_distribution<> distribution(bell_curves.size());

            bell_curves.erase(bell_curves.begin() +
                              distribution(rand_generator));
        }
    }

    // Generate modified landscape
    for (auto code = codes.begin(); code != codes.end(); code++)
    {
        auto length = get_length(*code);

        std::vector<double> &landscape      = probability_landscape[*code];
        std::vector<double> &orig_landscape = original_landscape[*code];

        for (int i = 0; i < length; i++)
        {
            landscape[i] = orig_landscape[i];
        }

        // Sum all bell curves
        for (auto curve = bell_curves.begin(); curve != bell_curves.end();
             curve++)
        {
            if (!code->compare(curve->location.chromosome))
            {
                for (int i = 0; i < length; i++)
                {
                    // Normal equation
                    landscape[i] +=
                        1. / (curve->sigma * 2.50662827463) *
                        exp(-0.5 *
                            pow(((i - curve->location.base) / curve->sigma),
                                2));
                }
            }
        }

        // Normalizing probabilities
        double max = 0;
        for (int i = 0; i < length; i++)
        {
            if (landscape[i] > max) max = landscape[i];
        }

        for (int i = 0; i < length; i++)
            landscape[i] /= max;
    }

    // TODO: Gene replacement and moving
}
void EvolutionDataProvider::snapshot(std::string folder)
{
    // TODO: Snapshot saving
}

void EvolutionDataProvider::die() { dead = true; }

bool EvolutionDataProvider::isdead() { return dead; }