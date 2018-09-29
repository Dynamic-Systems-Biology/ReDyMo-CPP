#include "genome.hpp"
#include "genomic_location.hpp"
#include <cstdlib>
#include <exception>
#include <iostream>
#include <vector>

Genome::Genome() {}

Genome::Genome(std::vector<std::shared_ptr<Chromosome>> &chromosomes)
{
    initialize(chromosomes);
}

void Genome::initialize(std::vector<std::shared_ptr<Chromosome>> &chromosomes)
{
    this->chromosomes = chromosomes;
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
    uint rand_chromosome = rand() % chromosomes.size();
    uint rand_base       = rand() % chromosomes[rand_chromosome]->size();
    return std::shared_ptr<GenomicLocation>(std::make_shared<GenomicLocation>(
        rand_base, chromosomes[rand_chromosome]));
}

std::shared_ptr<GenomicLocation> Genome::random_unreplicated_genomic_location()
{
    if (this->is_replicated())
        throw "There are no unreplicated bases available.";

    uint rand_chromosome, rand_base;

    do
        rand_chromosome = rand() % chromosomes.size();
    while (chromosomes[rand_chromosome]->is_replicated());

    do
        rand_base = rand() % chromosomes[rand_chromosome]->size();
    while (chromosomes[rand_chromosome]->base_is_replicated(rand_base));

    return std::shared_ptr<GenomicLocation>(std::make_shared<GenomicLocation>(
        rand_base, chromosomes[rand_chromosome]));
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
        n_interorigin_spaces += chromosome->get_n_fired_origins() + 1;

    return (double)this->size() / n_interorigin_spaces;
}

uint Genome::n_constitutive_origins()
{
    uint n_origins = 0;
    for (auto chromosome : chromosomes)
        n_origins += chromosome->n_constitutive_origins();
    return n_origins;
}