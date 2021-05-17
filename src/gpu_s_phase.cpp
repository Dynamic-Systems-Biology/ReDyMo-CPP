#include "gpu_s_phase.hpp"
#include "util.hpp"
#include <fstream>
#include <iostream>
#include <string>

#ifdef GPU_ENABLED

GPUSPhase::GPUSPhase(cl::Context &clContext, int origins_range, int n_resources,
                     int replication_speed, int timeout,
                     int transcription_period, bool has_dormant,
                     std::shared_ptr<DataProvider> data, std::string organism,
                     std::string name, std::string output_folder, int seed)
    : origins_range(origins_range), n_resources(n_resources),
      replication_speed(replication_speed), timeout(timeout),
      transcription_period(transcription_period), has_dormant(has_dormant),
      data(data), organism(organism), name(name), output_folder(output_folder),
      clContext(clContext)
{
    checkpoint_times.start_create = std::chrono::steady_clock::now();

    std::vector<std::shared_ptr<Chromosome>> chromosomes;

    auto codes = data->get_codes();
    for (auto code = codes.begin(); code != codes.end(); code++)
    {
        chromosomes.push_back(std::make_shared<Chromosome>(*code, *data));
    }

    genome = std::make_shared<Genome>(chromosomes, seed);

    fork_manager =
        std::make_shared<ForkManager>(n_resources, genome, replication_speed);

    initialize_opencl();

    checkpoint_times.end_create = std::chrono::steady_clock::now();
}

GPUSPhase::~GPUSPhase() {}

void GPUSPhase::initialize_opencl()
{
    // Read OpenCL Kernel Source File
    std::ifstream kernel_source_stream("../opencl/fork.ocl");
    kernel_source_stream.seekg(0, std::ios::end);
    kernel_source.reserve(kernel_source_stream.tellg());
    kernel_source_stream.seekg(0, std::ios::beg);

    kernel_source.assign((std::istreambuf_iterator<char>(kernel_source_stream)),
                         std::istreambuf_iterator<char>());

    // Build OpenCL Kernel
    kernel_program = cl::Program(clContext, kernel_source);

    try
    {
        kernel_program.build("-cl-std=CL2.0");
    }
    catch (...)
    {
        // Print build info for all devices
        cl_int buildErr = CL_SUCCESS;
        auto buildInfo =
            kernel_program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(&buildErr);

        std::cerr << "OpenCL failed to build Program!" << std::endl;
        for (auto &pair : buildInfo)
        {
            std::cerr << pair.second << std::endl << std::endl;
        }
    }
}

void GPUSPhase::simulate(int sim_number)
{
    checkpoint_times.start_sim = std::chrono::steady_clock::now();

    int alpha                     = 1;
    int time                      = 0;
    int constitutive_origins      = (int)genome->n_constitutive_origins();
    int n_collisions              = 0;
    bool use_constitutive_origins = origins_range > 0;

    std::cout << "[INFO] Starting GPU simulation " << sim_number << std::endl;

    std::cout << "[INFO] Loading data to OpenCL" << std::endl;

    try
    {
        // Commands
        cl::CommandQueue commands(clContext);

        //////////////////////////
        // Single value buffers //
        //////////////////////////
        cl::Buffer end_time(CL_MEM_READ_WRITE, sizeof(int));
        cl::Buffer replicated(CL_MEM_READ_WRITE, sizeof(int));
        cl::Buffer free_forks(CL_MEM_READ_WRITE, sizeof(int));

        commands.enqueueFillBuffer(end_time, 0, 0, sizeof(cl_int));
        commands.enqueueFillBuffer(replicated, 0, 0, sizeof(cl_int));
        commands.enqueueFillBuffer(free_forks, n_resources, 0, sizeof(cl_int));

        ////////////////////////////////
        // Start locations/directions //
        ////////////////////////////////
        cl::Buffer start_locations(CL_MEM_READ_WRITE,
                                   sizeof(int) * n_resources);
        cl::Buffer start_directions(CL_MEM_READ_WRITE,
                                    sizeof(int) * n_resources);

        commands.enqueueFillBuffer(start_locations, -1, 0,
                                   sizeof(cl_int) * n_resources);
        commands.enqueueFillBuffer(start_locations, 0, 0,
                                   sizeof(cl_int) * n_resources);

        /////////////////////////////
        // Collision count buffers //
        /////////////////////////////
        cl::Buffer rt_collisions(clContext, CL_MEM_READ_WRITE,
                                 sizeof(int) * genome->size());

        commands.enqueueFillBuffer(rt_collisions, 0, 0,
                                   sizeof(cl_int) * genome->size());

        ///////////////////////////////////////////
        // Replication timestamps for base pairs //
        ///////////////////////////////////////////
        cl::Buffer replication_times(clContext, CL_MEM_READ_WRITE,
                                     sizeof(unsigned int) * genome->size());

        commands.enqueueFillBuffer(replication_times, (unsigned int)0, 0,
                                   sizeof(cl_uint) * genome->size());

        ///////////////////////////////////////////
        // Probability landscape for full genome //
        ///////////////////////////////////////////
        cl::Buffer probability_landscape(clContext, CL_MEM_READ_ONLY,
                                         sizeof(double) * genome->size());

        std::vector<double> probabilities(genome->size(), 2);

        // TODO: Calculate probabilities vector

        commands.enqueueWriteBuffer(probability_landscape, CL_TRUE, 0,
                                    sizeof(double) * genome->size(),
                                    probabilities.data());

        ///////////////////////////////
        // Boundaries of chromosomes //
        ///////////////////////////////
        cl::Buffer chromosome_boundaries(clContext, CL_MEM_READ_WRITE,
                                         sizeof(int) *
                                             (genome->chromosomes.size() + 1));

        std::vector<int> boundaries(genome->chromosomes.size() + 1, 0);

        for (long unsigned int c = 0, bd = -1; c < chromosomes->size(); c++)
        {
            boundaries[c] = bd;
            bd += (*chromosomes)[c]->size();
        }

        commands.enqueueWriteBuffer(
            chromosome_boundaries, CL_TRUE, 0,
            sizeof(int) * (genome->chromosomes.size() + 1), boundaries.data());

        std::cout << "[INFO] Running OpenCL Kernel" << std::endl;

        ///////////////////////////////////////////
        // Replication timestamps for base pairs //
        ///////////////////////////////////////////
        std::vector<transcription_region_t> transcription_regions_v;

        for (int c = 0; c < chromosomes->size(); c++)
        {
            auto regions = (*chromosomes)[c]->get_transcription_regions();

            for (int r = 0; r < regions.size(); r++)
            {
                transcription_regions_v.push_back(regions[r]);
            }
        }

        cl::Buffer transcription_regions(clContext, CL_MEM_READ_WRITE,
                                         sizeof(transcription_region_t) *
                                             transcription_regions_v.size());

        commands.enqueueWriteBuffer(chromosome_boundaries, CL_TRUE, 0,
                                    sizeof(transcription_region_t) *
                                        transcription_regions_v.size(),
                                    transcription_regions_v.data());

        ///////////////////
        // Create kernel //
        ///////////////////

        //
        auto forkKernel =
            cl::KernelFunctor<cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer,
                              cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer,
                              cl::Buffer, unsigned int, unsigned int,
                              cl::Buffer, unsigned int, int, int, int,
                              unsigned int, int>(kernel_program, "fork");

        forkKernel(cl::EnqueueArgs(commands, cl::NDRange(n_resources)),
                   end_time, replicated, free_forks, start_locations,
                   start_directions, rt_collisions, replication_times,
                   probability_landscape, chromosome_boundaries,
                   transcription_period, transcription_regions_v.size(),
                   transcription_regions, timeout, genome->size(),
                   chromosomes->size(), n_resources, 30, 0);

        std::cout << "[INFO] " << sim_number << " Ended simulation"
                  << std::endl;

        ////////////////////////////
        // Retrieve data from OCL //
        ////////////////////////////
        int end_time_res   = 0;
        int replicated_res = 0;
        int free_forks_res = n_resources;
        std::vector<int> replication_times_res(genome->size(), 0);
        std::vector<int> ff_colisions_res(genome->size(), 0);
        std::vector<int> fb_colisions_res(genome->size(), 0);

        commands.enqueueReadBuffer(end_time, CL_TRUE, 0, sizeof(int),
                                   &end_time_res);
        commands.enqueueReadBuffer(replicated, CL_TRUE, 0, sizeof(int),
                                   &replicated_res);
        commands.enqueueReadBuffer(free_forks, CL_TRUE, 0, sizeof(int),
                                   &free_forks_res);
        commands.enqueueReadBuffer(replication_times, CL_TRUE, 0,
                                   sizeof(int) * genome->size(),
                                   replication_times_res.data());
        commands.enqueueReadBuffer(rt_collisions, CL_TRUE, 0,
                                   sizeof(int) * genome->size(),
                                   fb_colisions_res.data());

        std::cout << "[DEBUG] " << genome->size() << std::endl;

        std::cout << "[DEBUG] End: " << end_time_res << std::endl;
        std::cout << "[DEBUG] Replicated: " << replicated_res << std::endl;
        std::cout << "[DEBUG] Free: " << free_forks_res << std::endl;
    }
    catch (cl::Error e)
    {
        std::cerr << "[ERROR] " << e.what() << ": " << getCLErrorString(e.err())
                  << std::endl;
    }

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

void GPUSPhase::output(int sim_number, int time, int iod,
                       std::shared_ptr<Genome> genome)
{
    checkpoint_times.start_save = std::chrono::steady_clock::now();
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

    checkpoint_times.end_save = std::chrono::steady_clock::now();
}

void GPUSPhase::zstd_compression_output(int sim_number, int time, int iod,
                                        std::shared_ptr<Genome> genome,
                                        std::string path)
{
    // Write chromosome data
    for (auto chromosome : this->genome->chromosomes)
    {
        std::string code       = chromosome->get_code() + ".txt.zst";
        FILE *out_file         = fopen((path + code).c_str(), "w+b");
        size_t compressed_size = 0;
        void *compressed_data =
            compress_cpp_string(chromosome->to_string(), compressed_size);
        size_t written = fwrite(compressed_data, 1, compressed_size, out_file);
        if (written == 0) std::__throw_ios_failure("Failed to write to file.");
        fclose(out_file);
        free(compressed_data);
    }
}

void GPUSPhase::semantic_compression_output(int sim_number, int time, int iod,
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

const s_phase_checkpoints_t GPUSPhase::getTimes() const
{
    return checkpoint_times;
}

#endif