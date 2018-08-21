#include "replication_fork.hpp"
#include <memory>
#include <stdexcept>

ReplicationFork::ReplicationFork(Genome &genome, uint speed)
    : genome(genome), chromosome(*(Chromosome *)nullptr)
{
    this->speed      = speed;
    this->base       = -1;
    this->direction  = 0;
}

void ReplicationFork::attach(GenomicLocation &gen_loc, int direction, uint time)
{
    if (this->is_attached()) throw "This fork is already attached.";

    this->base       = gen_loc.base;
    this->chromosome = gen_loc.chromosome;
    this->direction  = direction;
    this->chromosome.replicate(this->base, this->base, time);
}

void ReplicationFork::unattach()
{
    this->base       = -1;
    this->direction  = 0;
    this->chromosome = *(Chromosome *)nullptr;
}

bool ReplicationFork::advance(uint time)
{
    int end_base = base + speed * direction;
    if (chromosome.replicate(base, end_base, time))
    {
        unattach();
        return false;
    }

    base = end_base;
    return true;
}

bool ReplicationFork::is_attached() { return !(base == -1 || direction == 0); }

int ReplicationFork::get_direction() { return direction; }

int ReplicationFork::get_base() { return base; }

std::shared_ptr<Chromosome> ReplicationFork::get_chromosome()
{
    return (std::shared_ptr<Chromosome>)&chromosome;
}