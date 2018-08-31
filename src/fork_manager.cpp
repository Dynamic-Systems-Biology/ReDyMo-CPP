#include "fork_manager.hpp"

ForkManager::ForkManager(uint size, Genome &genome, uint speed)
{
    this->n_forks      = size;
    this->n_free_forks = size;
    for (int i = 0; i < size; i++)
    {
        replication_forks.push_back(new ReplicationFork(genome, speed));
    }
}

uint ForkManager::check_replication_transcription_conflicts(uint time,
                                                            uint period,
                                                            bool has_dormant)
{
    uint n_collisions;

    uint RNAP_position = time % period;
    for (auto fork : replication_forks)
    {
        if (fork->is_attached())
        {
            Chromosome &chromosome = *fork->get_chromosome();
            for (auto region : chromosome.get_transcription_regions())
            {
                uint replisome_position_within_region = 0;
                uint region_size                      = 0;
                int RNAP_direction                    = 0;

                if (region.start < region.end)
                {
                    if (fork->get_base() < region.start ||
                        fork->get_base() > region.end)
                        continue;

                    replisome_position_within_region =
                        fork->get_base() - region.start;
                    RNAP_direction = 1;
                }
                else
                {
                    if (fork->get_base() > region.start ||
                        fork->get_base() < region.end)
                        continue;
                    replisome_position_within_region =
                        region.start - fork->get_base();

                    RNAP_direction = -1;
                }

                // Head to head collision!
                if (replisome_position_within_region % period ==
                        RNAP_position &&
                    fork->get_direction() != RNAP_direction)
                {
                    if (has_dormant)
                    {
                        chromosome.set_dormant_activation_probability(
                            fork->get_base());
                    }
                    fork->detach();
                    n_free_forks++;
                    n_collisions++;
                    // This fork collided, so there is no need to check other
                    // regions with it
                    break;
                }
            }
        }
    }
    return n_collisions;
}

void ForkManager::advance_attached_forks(uint time)
{
    for (auto fork : replication_forks)
    {
        if (fork->get_just_detached())
        {
            fork->set_just_detached(false);
            n_free_forks++;
        }
        else if (fork->is_attached())
            fork->advance(time);
    }
}

void ForkManager::attach_forks(GenomicLocation &location, uint time)
{
    if (n_free_forks < 2) return;

    uint n_forks_attached = 0;
    int direction         = 1;

    for (auto fork : replication_forks)
    {
        if (!fork->is_attached() && !fork->get_just_detached())
        {
            fork->attach(location, direction, time);
            n_forks_attached++;
            direction = -direction;
            if (n_forks_attached == 2) break;
        }
    }
    n_free_forks -= n_forks_attached;
}