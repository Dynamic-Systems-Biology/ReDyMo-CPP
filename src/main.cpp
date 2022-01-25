#include "configuration.hpp"
#include "evolution.hpp"
#include "gpu_s_phase.hpp"
#include "s_phase.hpp"
#include <algorithm>
#include <c4/yml/std/string.hpp>
#include <chrono>
#include <ctime>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <ryml.hpp>
#include <sstream>

int main(int argc, char *argv[])
{
    try
    {
        Configuration config(argc, argv);

        cl_configuration_data arg_values = config.arguments();

        srand(arg_values.seed);

        omp_set_num_threads(arg_values.threads);

        if (!arg_values.mode.compare("basic"))
        {

            std::shared_ptr<DataManager> data = std::make_shared<DataManager>(
                arg_values.organism, arg_values.data_dir + "/database.sqlite",
                arg_values.data_dir + "/MFA-Seq_" + arg_values.organism + "/",
                arg_values.probability);

            if (arg_values.gpu)
            {
#ifdef GPU_ENABLED
                std::cout << std::endl;

                // Select Platform
                std::vector<cl::Platform> platforms;
                cl::Platform::get(&platforms);
                cl::Platform platform;
                std::cout << "OpenCL Platforms:" << std::endl;
                for (uint i = 0; i < platforms.size(); i++)
                    std::cout << platforms[i].getInfo<CL_PLATFORM_NAME>()
                              << std::endl;
                std::cout << std::endl;

                for (auto &p : platforms)
                {
                    std::string platver = p.getInfo<CL_PLATFORM_VERSION>();
                    std::cout << p.getInfo<CL_PLATFORM_NAME>() << platver
                              << std::endl;
                    if ((platver.find("OpenCL 2.") != std::string::npos) ||
                        (platver.find("OpenCL 3.") != std::string::npos))
                    {
                        platform = p;
                    }
                }
                if (platform() == 0)
                {
                    std::cout << "No OpenCL 2.0 platform found." << std::endl;
                    return -1;
                }
                cl::Platform newP = cl::Platform::setDefault(platform);
                if (newP != platform)
                {
                    std::cout << "Error setting default platform." << std::endl;
                    return -1;
                }

                std::cout << "Using Platform:"
                          << platform.getInfo<CL_PLATFORM_NAME>() << std::endl;

                // Select Device
                std::vector<cl::Device> all_devices;
                platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);
                if (all_devices.size() == 0)
                {
                    std::cout
                        << " No devices found. Check OpenCL installation!\n";
                    exit(1);
                }

                cl::Device device = all_devices[0];
                std::cout << "Using Device: "
                          << device.getInfo<CL_DEVICE_NAME>() << std::endl
                          << std::endl;

                unsigned long long seed = arg_values.seed;

                // Start simulations
                // # pragma omp parallel for
                for (uint i = 0; i < arg_values.cells; i++)
                {
                    cl::Context context({device});

                    GPUSPhase s_phase(
                        context, arg_values.constitutive, arg_values.resources,
                        arg_values.speed, arg_values.timeout, arg_values.period,
                        arg_values.dormant, data, arg_values.organism,
                        arg_values.name, arg_values.output, i ^ seed);

                    s_phase.simulate(i);
                }
#endif
#ifndef GPU_ENABLED
                std::cout << "This code was not compiled with GPU support."
                          << std::endl;
#endif
            }
            else
            {
                std::vector<std::pair<int, s_phase_checkpoints_t>>
                    checkpoint_times;

                unsigned long long seed = arg_values.seed;

#pragma omp parallel for
                for (long long unsigned int i = 0; i < arg_values.cells; i++)
                {
                    // Run all simulations with the same parameters, except for
                    // seed, otherwise it would be exactly the same simulation
                    // every time.
                    SPhase *s_phase = new SPhase(
                        arg_values.constitutive, arg_values.resources,
                        arg_values.speed, arg_values.timeout, arg_values.period,
                        arg_values.dormant, data, arg_values.organism,
                        arg_values.name, arg_values.output, i ^ seed);
                    s_phase->simulate(i);

                    checkpoint_times.push_back(
                        std::pair<int, s_phase_checkpoints_t>(
                            i, s_phase->getTimes()));

                    delete s_phase;
                }

                // Calculate time statistics

                double created_sum =
                    0; // Total microseconds spent in creation phase
                double sim_sum =
                    0; // Total microseconds spent in simulation phase
                double saved_sum =
                    0; // Total microseconds spent in saving phase

                for (auto checkpoint : checkpoint_times)
                {
                    s_phase_checkpoints_t times = checkpoint.second;

                    created_sum +=
                        std::chrono::duration_cast<std::chrono::milliseconds>(
                            times.end_create - times.start_create)
                            .count();
                    sim_sum +=
                        std::chrono::duration_cast<std::chrono::milliseconds>(
                            times.end_sim - times.start_sim)
                            .count();
                    saved_sum +=
                        std::chrono::duration_cast<std::chrono::milliseconds>(
                            times.end_save - times.start_save)
                            .count();
                }

                // Calculate averages
                double created_avg = created_sum / arg_values.cells;
                double sim_avg     = sim_sum / arg_values.cells;
                double saved_avg   = saved_sum / arg_values.cells;

                std::cout << "[STAT] Average creation time   [ms] : "
                          << created_avg << std::endl;
                std::cout << "[STAT] Average simulation time [ms] : " << sim_avg
                          << std::endl;
                std::cout << "[STAT] Average saving time     [ms] : "
                          << saved_avg << std::endl;
                std::cout << "[STAT] Average s-phase time    [ms] : "
                          << created_avg + sim_avg + saved_avg << std::endl;
            }
        }
        else if (!arg_values.mode.compare("evolution"))
        {
            std::shared_ptr<EvolutionManager> evolution =
                std::make_shared<EvolutionManager>(config, arg_values.seed);

            evolution->run_all();
        }
    }
    catch (std::invalid_argument &e)
    {
        // This is how we recieve the messages from the argument parsing
        std::cout << e.what() << std::endl;
    }
}
