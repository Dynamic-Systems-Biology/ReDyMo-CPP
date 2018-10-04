#include "chromosome.hpp"
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

Chromosome::Chromosome() {}

Chromosome::Chromosome(
    std::string code, uint length, std::vector<double> &probability_landscape,
    std::vector<transcription_region_t> &transcription_regions,
    std::vector<constitutive_origin_t> &constitutive_origins)
{
    initialize(code, length, probability_landscape, transcription_regions,
               constitutive_origins);
}

void Chromosome::initialize(
    std::string code, uint length, std::vector<double> &probability_landscape,
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
}

uint Chromosome::size() { return this->length; }

std::string Chromosome::to_string()
{
    std::string chromosome_string = "";
    for (uint base = 0; base < length; base += 500)
    {
        chromosome_string += std::to_string(strand[base]);
        chromosome_string += "\n";
    }
    return chromosome_string;
}

bool Chromosome::base_is_replicated(uint base)
{
    if (base < 0 || base >= this->length)
        throw std::out_of_range("Given base is outside Chromosome length.");
    return this->strand[base] != -1;
}

double Chromosome::activation_probability(uint base)
{
    if (base < 0 || base >= this->length)
        throw std::out_of_range("Given base is outside Chromosome length.");
    return probability_landscape[base];
}

void Chromosome::set_dormant_activation_probability(uint base)
{
    if (base < 0 || base >= this->length)
        throw std::out_of_range("Given base is outside Chromosome length.");
    int c          = 10000;
    int left_base  = base - 2 * c;
    int right_base = base + 2 * c;
    left_base      = left_base < 0 ? 0 : left_base;
    right_base =
        right_base > (int)this->length ? (int)this->length : right_base;

    for (int curr_base = left_base; curr_base < right_base; curr_base++)
    {
        int offset            = curr_base - base;
        double gaussian_value = exp(-pow(offset, 2) / (2 * pow(c, 2)));
        probability_landscape[curr_base] += gaussian_value;
        if (probability_landscape[curr_base] > 1)
            probability_landscape[curr_base] = 1;
    }
}

bool Chromosome::replicate(int start, int end, int time)
{
    if (start < 0 || start > (int)this->length)
        throw std::out_of_range("The start base is not inside the Chromosome");

    // A non normal replication refers to a replication which overlaps areas
    // already replicated or which the end base is outside the Chromosome(less
    // than zero or greater than the Chromosome itself).
    bool normal_replication = true;

    if (end < 0 || end >= (int)strand.size())
    {
        // limit the end to the strand size
        end = end < 0 ? 0 : end;
        end = end >= (int)strand.size() ? (int)strand.size() : end;
        normal_replication = false;
    }

    bool inverted = end < start;
    for (int base = start; inverted ? base > end - 1 : base < end + 1;
         inverted ? base-- : base++)
    {
        if (strand[base] == -1)
        {
            strand[base] = time;
            n_replicated_bases++;
        }
        else if (base != start)
        {
            printf("bad replication at %d from %d to %d\n", base, start, end);
            normal_replication = false;
            break;
        }
    }

    return normal_replication;
}

bool Chromosome::is_replicated()
{
    return this->n_replicated_bases == this->length;
}

std::string Chromosome::get_code() { return this->code; }

uint Chromosome::n_constitutive_origins()
{
    return this->constitutive_origins.size();
}

uint Chromosome::get_n_replicated_bases() { return n_replicated_bases; }

uint Chromosome::get_n_fired_origins() { return n_fired_origins; }

void Chromosome::add_fired_origin() { n_fired_origins++; }

std::vector<transcription_region_t> &Chromosome::get_transcription_regions()
{
    return transcription_regions;
}

bool Chromosome::operator==(Chromosome &other)
{
    return this->code == other.get_code();
}
