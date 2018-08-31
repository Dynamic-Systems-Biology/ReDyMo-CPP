#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include "chromosome.hpp"
#include "genomic_location.hpp"

GenomicLocation::GenomicLocation(uint base, Chromosome &chromosome)
    : chromosome(chromosome)
{
    if (base >= chromosome.size())
        throw std::invalid_argument("Base is not inside given chromosome.");

    this->base = base;
}

bool GenomicLocation::is_replicated()
{
    return this->chromosome.base_is_replicated(this->base);
}

bool GenomicLocation::will_activate(bool use_constitutive_origin,
                                    uint origins_range)
{
    if (!use_constitutive_origin)
    {
        float chance = (float)(rand() % 100) / 100;
        return chance < this->chromosome.activation_probability(this->base);
    }
    std::vector<constitutive_origin_t> not_fired_origins;

    for (auto origin : chromosome.constitutive_origins)
        if (!chromosome.base_is_replicated(origin.base))
            not_fired_origins.push_back(origin);

    for (auto origin : not_fired_origins)
        if (this->base >= (origin.base - origins_range / 2) &&
            this->base <= (origin.base + origins_range / 2))
            return true;

    return false;
}

constitutive_origin_t *
GenomicLocation::get_constitutive_origin(int origins_range)
{
    constitutive_origin_t *ret_origin;
    std::vector<constitutive_origin_t *> not_fired_origins;
    origins_range /= 2;

    for (auto origin : chromosome.constitutive_origins)
        if (!chromosome.base_is_replicated(origin.base))
            not_fired_origins.push_back(&origin);

    for (auto origin : not_fired_origins)
    {
        if ((int)this->base >= ((int)origin->base - origins_range) &&
            (int)this->base <= ((int)origin->base + origins_range))
        {
            ret_origin = origin;
            return ret_origin;
        }
    }
    // If code reaches here, no nonfired origin was found
    throw std::range_error("There are no origins within range.");
}

bool GenomicLocation::put_fired_constitutive_origin(
    constitutive_origin_t &origin)
{
    bool found = false;
    for (auto curr_origin : chromosome.constitutive_origins)
        if (curr_origin == origin) found = true;

    if (found)
        chromosome.fired_constitutive_origins.push_back(origin);
    else
        throw std::invalid_argument("Origin not found.");

    return found;
}