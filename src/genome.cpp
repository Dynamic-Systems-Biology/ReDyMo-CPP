#include "genome.hpp"
#include "genomic_location.hpp"
#include <cstdlib>
#include <exception>
#include <iostream>
#include <vector>

Genome::Genome(std::vector<std::shared_ptr<Chromosome>> &chromosomes,
               unsigned long long seed)
    : seed(seed)
{
    this->rand_generator = std::mt19937(seed);
    initialize(chromosomes);
}

void Genome::initialize(std::vector<std::shared_ptr<Chromosome>> &chromosomes)
{
    std::vector<int> chromosome_sizes;

    for (auto chromosome = 0; chromosome < chromosomes.size(); chromosome++)
    {
        this->chromosomes.push_back(chromosomes[chromosome]);
        chromosome_sizes.push_back(chromosomes[chromosome]->size());
    }

    this->chromosome_distribution = std::discrete_distribution<int>(
        chromosome_sizes.begin(), chromosome_sizes.end());
}

uint Genome::size()
{
    uint tot_length = 0;
    for (auto chromosome : chromosomes)
        tot_length += chromosome->size();
    return tot_length;
}

std::shared_ptr<GenomicLocation> Genome::random_genomic_location()
{
    uint rand_chromosome = chromosome_distribution(rand_generator);

    std::uniform_int_distribution<int>::param_type bases_dist(
        0, chromosomes[rand_chromosome]->size() - 1);

    // Change the distribution parameters to match the numer of possible bases
    base_distribution.param(bases_dist);

    uint rand_base = base_distribution(rand_generator);
    return std::make_shared<GenomicLocation>(
        rand_base, chromosomes[rand_chromosome], this->seed);
}

// Actually never used, still here for eventual future use
std::shared_ptr<GenomicLocation> Genome::random_unreplicated_genomic_location()
{
    if (this->is_replicated())
        throw std::runtime_error("There are no unreplicated bases available.");

    uint rand_chromosome, rand_base;

    do
        rand_chromosome = rand() % chromosomes.size();
    while (chromosomes[rand_chromosome]->is_replicated());

    do
        rand_base = rand() % chromosomes[rand_chromosome]->size();
    while (chromosomes[rand_chromosome]->base_is_replicated(rand_base));

    return std::make_shared<GenomicLocation>(
        rand_base, chromosomes[rand_chromosome], this->seed);
}

bool Genome::is_replicated()
{
    bool replicated = true;
    for (auto chromosome : chromosomes)
        if (!chromosome->is_replicated()) replicated = false;

    return replicated;
}

double Genome::average_interorigin_distance()
{
    uint n_interorigin_spaces = 0;
    for (auto chromosome : chromosomes)
    {
        n_interorigin_spaces += chromosome->get_n_fired_origins() + 1;
    }
    if (n_interorigin_spaces == 0) return 0;
    return (double)this->size() / n_interorigin_spaces;
}

uint Genome::n_constitutive_origins()
{
    uint n_origins = 0;
    for (auto chromosome : chromosomes)
    {
        n_origins += chromosome->n_constitutive_origins();
    }
    return n_origins;
}
