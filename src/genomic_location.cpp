#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include "chromosome.hpp"
#include "genomic_location.hpp"

std::uniform_real_distribution<double> GenomicLocation::rand_distribution =
    std::uniform_real_distribution<double>(0, 1);

GenomicLocation::GenomicLocation(uint base,
                                 std::shared_ptr<Chromosome> chromosome,
                                 std::mt19937 *rand_generator)
    : chromosome(chromosome), rand_generator(rand_generator)
{
    if (base >= chromosome->size())
        throw std::invalid_argument("Base is not inside given chromosome.");

    this->base = base;
}

bool GenomicLocation::is_replicated()
{
    return this->chromosome->base_is_replicated(this->base);
}

bool GenomicLocation::will_activate(bool use_constitutive_origin,
                                    uint origins_range)
{
    if (!use_constitutive_origin)
    {
        double chance = rand_distribution(*rand_generator);
        return chance < this->chromosome->activation_probability(this->base);
    }

    bool origin_in_range = get_constitutive_origin(origins_range).base != -1;

    return origin_in_range;
}

constitutive_origin_t
GenomicLocation::get_constitutive_origin(int origins_range)
{
    constitutive_origin_t ret_origin;
    std::vector<constitutive_origin_t> not_fired_origins;
    origins_range /= 2;

    for (auto origin : *chromosome->constitutive_origins)
    {
        bool fired = false;
        for (auto fired_origin : *chromosome->fired_constitutive_origins)
        {
            if (fired_origin.base == origin.base)
                fired = true;
        }
        if (!fired)
            not_fired_origins.push_back(origin);
    }
    for (auto origin : not_fired_origins)
    {
        if ((int)this->base >= ((int)origin.base - origins_range) &&
            (int)this->base <= ((int)origin.base + origins_range))
        {
            ret_origin = origin;
            return ret_origin;
        }
    }
    // If code reaches here, no nonfired origin was found
    ret_origin.base = -1;
    return ret_origin;
}

bool GenomicLocation::put_fired_constitutive_origin(
    constitutive_origin_t origin)
{
    for (auto curr_origin : *chromosome->constitutive_origins)
    {
        if (curr_origin == origin)
        {
            //TODO: we really need to check if the origin is being fired twice?
//            for (auto fired_origin : *chromosome->fired_constitutive_origins)
//                if (fired_origin == origin) return false;
            chromosome->fired_constitutive_origins->push_back(origin);
            return true;
        }
    }
    throw std::runtime_error("Method put_fired_constitutive_origin: Could not register fired constitutive origin. Invalid base or already fired.");
}

GenomicLocation &GenomicLocation::operator+=(int bases)
{
    long tmp = (long)base + bases;
    if (tmp < 0)
        tmp = 0;
    else if (tmp >= chromosome->length)
        tmp = chromosome->length - 1;

    base = (uint)tmp;

    return *this;
}

GenomicLocation GenomicLocation::operator+(int bases)
{
    long tmp = (long)base + bases;
    if (tmp < 0)
        tmp = 0;
    else if (tmp >= chromosome->length)
        tmp = chromosome->length - 1;

    return GenomicLocation(tmp, chromosome, rand_generator);
}
