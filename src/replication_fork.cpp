#include "replication_fork.hpp"
#include <memory>
#include <stdexcept>

ReplicationFork::ReplicationFork(Genome &genome, uint speed)
    : genome(genome), chromosome(*(Chromosome *)nullptr)
{
    this->speed         = speed;
    this->base          = -1;
    this->direction     = 0;
    this->just_detached = false;
}

void ReplicationFork::attach(GenomicLocation &gen_loc, int direction, uint time)
{
    if (this->is_attached()) throw "This fork is already attached.";

    this->base       = gen_loc.base;
    this->chromosome = gen_loc.chromosome;
    this->direction  = direction;
    this->chromosome.replicate(this->base, this->base, time);
}

void ReplicationFork::detach()
{
    this->base          = -1;
    this->direction     = 0;
    this->chromosome    = *(Chromosome *)nullptr;
    this->just_detached = true;
}

bool ReplicationFork::advance(uint time)
{
    int end_base = base + speed * direction;
    if (chromosome.replicate(base, end_base, time))
    {
        detach();
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

bool ReplicationFork::get_just_detached() { return just_detached; }

void ReplicationFork::set_just_detached(bool new_value)
{
    just_detached = new_value;
}