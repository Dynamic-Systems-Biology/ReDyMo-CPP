#include "genome.hpp"
#include "genomic_location.hpp"
#include <cstdlib>
#include <exception>
#include <iostream>
#include <vector>

Genome::Genome() {}

Genome::Genome(std::vector<Chromosome> &chromosomes)
{
    initialize(chromosomes);
}

void Genome::initialize(std::vector<Chromosome> &chromosomes)
{
    this->chromosomes = chromosomes;
}

unsigned int Genome::size()
{
    unsigned int tot_length = 0;
    for (auto chromosome : chromosomes)
        tot_length += chromosome.size();
    return tot_length;
}

GenomicLocation Genome::random_genomic_location()
{
    unsigned int rand_chromosome = rand() % chromosomes.size();
    unsigned int rand_base       = rand() % chromosomes[rand_chromosome].size();
    try
    {
        GenomicLocation loc(rand_base, chromosomes[rand_chromosome]);
        return loc;
    }
    catch (std::exception e)
    {
        std::cout << e.what << std::endl;
    }
}

GenomicLocation Genome::random_unreplicated_genomic_location()
{
    if (this->is_replicated())
        throw "There are no unreplicated bases available.";

    unsigned int rand_chromosome, rand_base;

    do
        rand_chromosome = rand() % chromosomes.size();
    while (chromosomes[rand_chromosome].is_replicated());

    do
        rand_base = rand() % chromosomes[rand_chromosome].size();
    while (chromosomes[rand_chromosome].base_is_replicated(rand_base));

    try
    {
        GenomicLocation loc(rand_base, chromosomes[rand_chromosome]);
        return loc;
    }
    catch (std::exception e)
    {
        std::cout << e.what << std::endl;
    }
}

bool Genome::is_replicated()
{
    bool replicated = true;
    for (auto chromosome : chromosomes)
        if (!chromosome.is_replicated()) replicated == false;

    return replicated;
}

float Genome::average_interorigin_distance()
{
    unsigned int n_interorigin_spaces = 0;
    for (auto chromosome : chromosomes)
        n_interorigin_spaces += chromosome.n_fired_origins + 1;

    return (float)this->size() / n_interorigin_spaces;
}

unsigned int Genome::n_constitutive_origins()
{
    unsigned int n_origins = 0;
    for (auto chromosome : chromosomes)
        n_origins += chromosome.n_constitutive_origins();
    return n_origins;
}