#include "s_phase.hpp"
#include <fstream>
#include <iostream>
#include <string>

SPhase::SPhase(int origins_range, int n_resources, int replication_speed,
               int timeout, int transcription_period, bool has_dormant,
               std::shared_ptr<DataManager> data, std::string organism)
    : origins_range(origins_range), n_resources(n_resources),
      replication_speed(replication_speed), timeout(timeout),
      transcription_period(transcription_period), has_dormant(has_dormant),
      data(data), organism(organism)
{
    std::vector<std::shared_ptr<Chromosome>> tmp =
        data->get_chromosome_data(organism);
    chromosomes = &tmp;
    genome      = std::make_shared<Genome>(*chromosomes);
    fork_manager =
        std::make_shared<ForkManager>(n_resources, genome, replication_speed);
}

SPhase::~SPhase() {}

void SPhase::simulate(int sim_number)
{

    int alpha                     = 1;
    int time                      = 0;
    int constitutive_origins      = (int)genome->n_constitutive_origins();
    int n_collisions              = 0;
    bool use_constitutive_origins = origins_range > 0;

    std::cout << "[INFO] Starting simulation " << sim_number << std::endl;
    while (!genome->is_replicated() && time < timeout &&
           !(constitutive_origins == 0 &&
             (int)fork_manager->n_free_forks == n_resources))
    {
        time++;

        // Advance the forks
        fork_manager->advance_attached_forks(time);

        // Check for collisions
        if (transcription_period > 0)
            n_collisions +=
                fork_manager->check_replication_transcription_conflicts(
                    time, transcription_period, has_dormant);

        // At an alpha iteration it makes one attempt for each detached fork
        if (time % alpha == 0 && !genome->is_replicated())
        {
            int n_forks = (int)fork_manager->n_free_forks;
            for (int i = 0; i < n_forks; i++)
            {
                GenomicLocation loc = *genome->random_genomic_location();

                if (!loc.is_replicated() && fork_manager->n_free_forks >= 2 &&
                    loc.will_activate(use_constitutive_origins, origins_range))
                {
                    fork_manager->attach_forks(loc, time);
                    if (use_constitutive_origins)
                    {
                        constitutive_origin_t origin =
                            loc.get_constitutive_origin(origins_range);
                        if (!loc.put_fired_constitutive_origin(origin))
                            constitutive_origins += 0;
                        // std::cout << "[WARN] Failed to add origin. "
                        //  "Simulation "
                        //   << sim_number << std::endl;
                        constitutive_origins--;
                    }
                }
            }
        }
    }

    std::cout << "[INFO] " << sim_number << " Ended simulation" << std::endl;

    if (genome->is_replicated())
        std::cout << "\t[INFO] " << sim_number << " Genome fully replicated."
                  << std::endl;
    else if (time == timeout)
        std::cout << "\t[WARN] " << sim_number
                  << " Timeout simulation: " << std::endl;

    std::cout << "\t[INFO] " << sim_number
              << " Number of Collisions: " << n_collisions << std::endl;
    if (use_constitutive_origins)
        std::cout << "\t[INFO] " << sim_number
                  << " Number of constitutive origins that did not fire: "
                  << constitutive_origins << std::endl;
    output(sim_number, time, genome->average_interorigin_distance(), genome);
}

void SPhase::output(int sim_number, int time, int iod,
                    std::shared_ptr<Genome> genome)
{
    system("mkdir -p output");
    std::string dir = "output/";
    if (has_dormant)
        dir += "true";
    else
        dir += "false";
    dir += "_";
    dir += std::to_string(n_resources);
    dir += "_";
    dir += std::to_string(transcription_period);
    dir += "/";
    system(("mkdir -p " + dir).c_str());
    std::string simulation = "simulation_" + std::to_string(sim_number);
    simulation += "/";
    system(("mkdir  -p " + dir + simulation).c_str());

    std::ofstream output_file;
    output_file.open((dir + simulation + "cell.txt").c_str());
    output_file << n_resources << "\t" << this->replication_speed << "\t"
                << time << "\t" << iod << "\t\n";
    output_file.close();
    for (auto chromosome : this->genome->chromosomes)
    {
        std::string code = chromosome->get_code() + ".txt.zst";
        FILE *out_file   = fopen((dir + simulation + code).c_str(), "w+b");
        size_t compressed_size = 0;
        void *compressed_data =
            compress_cpp_string(chromosome->to_string(), compressed_size);
        size_t written = fwrite(compressed_data, 1, compressed_size, out_file);
        if (written == 0) std::__throw_ios_failure("Failed to write to file.");
        fclose(out_file);
        free(compressed_data);
    }
}
