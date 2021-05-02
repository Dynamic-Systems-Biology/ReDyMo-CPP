#include "fork_manager.hpp"
#include <iostream>

ForkManager::ForkManager(uint n_forks, std::shared_ptr<Genome> genome,
                         uint speed)
{
    this->n_forks                         = n_forks;
    this->n_free_forks                    = n_forks;
    this->metric_times_attached           = 0;
    this->metric_times_detached_normal    = 0;
    this->metric_times_detached_collision = 0;
    for (int i = 0; i < (int)n_forks; i++)
    {
        replication_forks.push_back(
            std::make_shared<ReplicationFork>(genome, this, speed));
    }
}

uint ForkManager::check_replication_transcription_conflicts(uint time,
                                                            uint period,
                                                            bool has_dormant)
{
    uint n_collisions = 0;

    uint RNAP_position = time % period;
    for (auto fork : replication_forks)
    {
        if (fork->is_attached())
        {
            std::shared_ptr<Chromosome> chromosome = fork->get_chromosome();
            for (auto region : *chromosome->get_transcription_regions())
            {
                uint replisome_position_within_region = 0;
                int RNAP_direction                    = 0;

                if (region.start < region.end)
                {
                    if (fork->get_base() < (int)region.start ||
                        fork->get_base() > (int)region.end)
                        continue;

                    replisome_position_within_region =
                        fork->get_base() - region.start;
                    RNAP_direction = 1;
                }
                else
                {
                    if (fork->get_base() > (int)region.start ||
                        fork->get_base() < (int)region.end)
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
                        chromosome->set_dormant_activation_probability(
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
    if (n_forks_attached == 2) location.chromosome->add_fired_origin();
    metric_times_attached += n_forks_attached;
    n_free_forks -= n_forks_attached;
}
