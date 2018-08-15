#include "chromosome.hpp"
#include <cmath>
#include <stdexcept>
#include <string>

Chromosome::Chromosome() {}

Chromosome::Chromosome(
    unsigned int code, unsigned int length,
    std::vector<float> &probability_landscape,
    std::vector<transcription_region_t> &transcription_regions,
    std::vector<constitutive_origin_t> &constitutive_origins)
{
    initialize(code, length, probability_landscape, transcription_regions,
               constitutive_origins);
}

void Chromosome::initialize(
    unsigned int code, unsigned int length,
    std::vector<float> &probability_landscape,
    std::vector<transcription_region_t> &transcription_regions,
    std::vector<constitutive_origin_t> &constitutive_origins)
{
    if (length <= 0)
        throw std::invalid_argument("Given length is not a positive number.");
    this->code                  = code;
    this->length                = length;
    this->probability_landscape = probability_landscape;
    this->transcription_regions = transcription_regions;
    this->constitutive_origins  = constitutive_origins;
    this->n_replicated_bases    = 0;
    this->n_fired_origins       = 0;
    this->strand.resize(length, -1);
    this->fired_constitutive_origins.resize(length, false);
}

unsigned int Chromosome::size() { return this->length; }

std::string Chromosome::to_string()
{
    std::string chromosome_string = "";
    for (int base = 0; base < length; base += 500)
    {
        chromosome_string += strand[base];
    }
    return chromosome_string;
}

bool Chromosome::base_is_replicated(unsigned int base)
{
    if (base < 0 || base >= this->length)
        throw std::out_of_range("Given base is outside Chromosome length.");
    return this->strand[base] != -1;
}

float Chromosome::activation_probability(unsigned int base)
{
    if (base < 0 || base >= this->length)
        throw std::out_of_range("Given base is outside Chromosome length.");
    return probability_landscape[base];
}

void Chromosome::set_dormant_activation_probability(unsigned int base)
{
    if (base < 0 || base >= this->length)
        throw std::out_of_range("Given base is outside Chromosome length.");
    int c          = 10000;
    int left_base  = base - 2 * c;
    int right_base = base + 2 * c;
    left_base      = left_base < 0 ? 0 : left_base;
    right_base     = right_base > this->length ? this->length : right_base;

    for (int curr_base = left_base; curr_base < right_base; curr_base++)
    {
        int offset           = curr_base - base;
        float gaussian_value = exp(-pow(offset, 2) / (2 * pow(c, 2)));
        probability_landscape[curr_base] += gaussian_value;
        if (probability_landscape[curr_base] > 1)
            probability_landscape[curr_base] = 1;
    }
}

bool Chromosome::replicate(int start, int end, int time)
{
    if (start < 0 || start > this->length)
        throw std::out_of_range("The start base is not inside the Chromosome");

    // A non normal replication refers to a replication which overlaps areas
    // already replicated or which the end base is outside the Chromosome(less
    // than zero or greater than the Chromosome itself).
    bool normal_replication = true;

    if (end < 0 || end >= strand.size())
    {
        end                = end < 0 ? 0 : strand.size() - 1;
        normal_replication = true;
    }

    for (unsigned int base = start; base != end;
         base += (end - start > 0) ? 1 : -1)
    {
        if (strand[base] == -1)
        {
            strand[base] = time;
            n_replicated_bases++;
        }
    }

    return normal_replication;
}

bool Chromosome::is_replicated()
{
    return this->n_replicated_bases == this->length;
}

unsigned int Chromosome::get_code() { return this->code; }

unsigned int Chromosome::n_constitutive_origins()
{
    return this->constitutive_origins.size();
}
