#include "s_phase.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

SPhase::SPhase(int origins_range, int n_resources, int replication_speed,
               int timeout, int transcription_period, bool has_dormant,
               std::shared_ptr<DataProvider> data, std::string organism,
               std::string name, std::string output_folder,
               unsigned long long seed)
    : origins_range(origins_range), n_resources(n_resources),
      replication_speed(replication_speed), timeout(timeout),
      transcription_period(transcription_period), has_dormant(has_dormant),
      data(data), organism(organism), name(name), output_folder(output_folder)
{
    checkpoint_times.start_create = std::chrono::steady_clock::now();

    std::vector<std::shared_ptr<Chromosome>> chromosomes;

    auto codes = data->get_codes();
    for (auto code = codes.begin(); code != codes.end(); code++)
    {
        chromosomes.push_back(std::make_shared<Chromosome>(*code, data));
    }

    genome = std::make_shared<Genome>(chromosomes, seed);

    fork_manager =
        std::make_shared<ForkManager>(n_resources, genome, replication_speed);

    checkpoint_times.end_create = std::chrono::steady_clock::now();
}

SPhase::SPhase(Configuration &configuration, std::shared_ptr<DataProvider> data,
               unsigned long long seed)
    : data(data)
{
    auto args = configuration.arguments();

    origins_range        = args.constitutive;
    n_resources          = args.resources;
    replication_speed    = args.speed;
    timeout              = args.timeout;
    transcription_period = args.period;
    has_dormant          = args.dormant;
    organism             = args.organism;
    name                 = args.name;
    output_folder        = args.output;

    checkpoint_times.start_create = std::chrono::steady_clock::now();

    std::vector<std::shared_ptr<Chromosome>> chromosomes;

    auto codes = data->get_codes();
    for (auto code = codes.begin(); code != codes.end(); code++)
    {
        chromosomes.push_back(std::make_shared<Chromosome>(*code, data));
    }

    genome = std::make_shared<Genome>(chromosomes, seed);

    fork_manager =
        std::make_shared<ForkManager>(n_resources, genome, replication_speed);

    checkpoint_times.end_create = std::chrono::steady_clock::now();
}

SPhase::~SPhase() {}

simulation_stats SPhase::get_stats() { return stats; }

void SPhase::simulate(int sim_number)
{

    checkpoint_times.start_sim = std::chrono::steady_clock::now();

    int alpha                     = 1;
    int time                      = 0;
    int constitutive_origins      = (int)genome->n_constitutive_origins();
    int n_collisions              = 0;
    bool use_constitutive_origins = origins_range > 0;

    std::cout << "[INFO] Starting simulation " << sim_number << std::endl
              << std::flush;
    while (!genome->is_replicated() && time < timeout &&
           !(use_constitutive_origins && constitutive_origins == 0 &&
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

    stats.time       = time;
    stats.collisions = n_collisions;

    std::cout << "[INFO] " << sim_number << " Ended simulation" << std::endl;

    if (genome->is_replicated())
        std::cout << "\t[INFO] " << sim_number
                  << " Genome fully replicated at time " << time << "."
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

    checkpoint_times.end_sim = std::chrono::steady_clock::now();

    output(sim_number, time, genome->average_interorigin_distance(), genome);
}

void SPhase::output(int sim_number, int time, int iod,
                    std::shared_ptr<Genome> genome)
{
    checkpoint_times.start_save = std::chrono::steady_clock::now();

    // Create simulation folder
    std::stringstream folder_name_stream;

    folder_name_stream << output_folder << "/" << name << "_"
                       << (has_dormant ? "true" : "false") << "_"
                       << std::to_string(n_resources) << "_"
                       << std::to_string(transcription_period) << "/";

    std::string dir = folder_name_stream.str();
    std::string simulation =
        this->name + "_" + "simulation_" + std::to_string(sim_number) + "/";

    system(("mkdir -p " + output_folder).c_str());

    system(("mkdir -p " + dir).c_str());
    system(("mkdir  -p " + dir + simulation).c_str());

    // Create Metadata File
    std::ofstream output_file;
    output_file.open((dir + simulation + "cell.txt").c_str());
    output_file << n_resources << "\t" << this->replication_speed << "\t"
                << time << "\t" << iod << "\t\n";
    output_file.close();

    // Save Chromosome data
    semantic_compression_output(sim_number, time, iod, genome,
                                dir + simulation);

    // Write chromosome data

    /*for (auto chromosome_ptr : this->genome->chromosomes)
    {
        // Get chromosome reference
        auto &chromosome = *chromosome_ptr;

        // Make filename
        std::string code = chromosome.get_code() + ".cseq";

        std::ofstream output_file;
        output_file.open((dir + simulation + code).c_str());

        output_file << chromosome.to_string();

        output_file.close();
    }*/

    checkpoint_times.end_save = std::chrono::steady_clock::now();
}

void SPhase::semantic_compression_output(int sim_number, int time, int iod,
                                         std::shared_ptr<Genome> genome,
                                         std::string path)
{
    // Makes formatted string for compression
    auto output_str = [](int start_value, int end_value, int seq_length) {
        std::string out = std::to_string(start_value);
        if (end_value != INT32_MIN && end_value != start_value)
            out += "-" + std::to_string(end_value);
        if (seq_length != 1) out += "x" + std::to_string(seq_length);

        return out;
    };

    // Write chromosome data
    for (auto chromosome_ptr : this->genome->chromosomes)
    {
        // Get chromosome reference
        auto &chromosome = *chromosome_ptr;

        // Make filename
        std::string code = chromosome.get_code() + ".cseq";

        // Open file for writing
        std::ofstream output_file;
        output_file.open((path + code).c_str());

        // Current and last two number streaks
        struct number_streak
        {
            int value   = INT32_MIN;
            int length  = INT32_MIN;
            bool in_seq = false;
        } number_streaks[3], null_streak;

        // Sequence data
        struct sequence_data
        {
            int start_value = INT32_MIN;
            int direction   = 0;
        } sequence, null_sequence;

        // Cache chromosome size
        const int chromosome_size = chromosome.size();

        for (int bp = 0; bp < chromosome_size + 2; bp++)
        {
            // Retrieve value if in range
            int value = INT32_MIN;
            if (bp < chromosome_size) value = chromosome[bp];

            // If a new number streak has started
            if (value == INT32_MIN || value != number_streaks[0].value)
            {
                // Finalize sequence if unable to continue, like when streak
                // changes size or step size
                if (sequence.start_value != INT32_MIN &&
                    (number_streaks[0].length != number_streaks[1].length ||
                     sequence.direction !=
                         number_streaks[0].value - number_streaks[1].value))
                {
                    // Write output for this sequence
                    output_file << output_str(sequence.start_value,
                                              number_streaks[1].value,
                                              number_streaks[1].length)
                                << std::endl;

                    // Zero sequence data
                    sequence = null_sequence;

                    // Set as sequence
                    number_streaks[1].in_seq = true;
                }
                // Start sequence if not in sequence, current streak is valid,
                // and start of sequence is valid
                else if (sequence.start_value == INT32_MIN &&
                         number_streaks[0].value != INT32_MIN &&
                         number_streaks[0].length != INT32_MIN &&
                         number_streaks[0].length == number_streaks[1].length &&
                         abs(number_streaks[0].value -
                             number_streaks[1].value) == 1)
                {
                    // Create sequence
                    sequence.start_value = number_streaks[1].value;
                    sequence.direction =
                        number_streaks[0].value - number_streaks[1].value;

                    // Set as sequence
                    number_streaks[1].in_seq = true;
                }
                // Set as sequence if in sequence
                else if (sequence.start_value != INT32_MIN)
                    number_streaks[1].in_seq = true;

                // If it's a unique value streak (not a sequence)
                if (!number_streaks[1].in_seq &&
                    number_streaks[1].value != INT32_MIN)
                {
                    // Write output for this value
                    output_file
                        << output_str(number_streaks[1].value, INT32_MIN,
                                      number_streaks[1].length)
                        << std::endl;
                }

                // Shift number streaks
                number_streaks[2] = number_streaks[1];
                number_streaks[1] = number_streaks[0];
                number_streaks[0] = number_streak{value, 0, false};
            }

            if (value != INT32_MIN) number_streaks[0].length++;
        }

        output_file.close();
    }
}

const s_phase_checkpoints_t SPhase::getTimes() const
{
    return checkpoint_times;
}
