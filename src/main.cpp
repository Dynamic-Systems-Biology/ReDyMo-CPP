#include "s_phase.hpp"
#include <algorithm>
#include <ctime>
#include <iostream>
#include <sstream>
#include <omp.h>

int get_cmd_option(char **input, int size, const std::string &option)
{
    for (int i = 0; i < size; i++)
        if (option == input[i]) return ++i;
    return -1;
}

bool cmd_option_exists(char **begin, char **end, const std::string &option)
{
    return std::find(begin, end, option) != end;
}

int main(int argc, char *argv[])
{
    std::vector<std::string> arg_options = {"--cells",     "--organism",
                                            "--resources", "--speed",
                                            "--timeout",   "--dormant"};
    for (auto option : arg_options)
    {
        if (!cmd_option_exists(argv, argv + argc, option))
        {
            std::cout << "Incorrect arguments. For usage see README.md" << std::endl;
            exit(0);
        }
    }

    int index            = get_cmd_option(argv, argc, arg_options[0]);
    int n_cells          = atoi(argv[index]);
    index                = get_cmd_option(argv, argc, arg_options[1]);
    std::string organism = argv[index];
    index                = get_cmd_option(argv, argc, arg_options[2]);
    int n_resources      = atoi(argv[index]);
    index                = get_cmd_option(argv, argc, arg_options[3]);
    int speed            = atoi(argv[index]);
    index                = get_cmd_option(argv, argc, arg_options[4]);
    int timeout          = atoi(argv[index]);
    index                = get_cmd_option(argv, argc, arg_options[5]);

    bool dormant;
    std::stringstream ss(argv[index]);
    ss >> std::boolalpha >> dormant;

    std::shared_ptr<DataManager> data = std::make_shared<DataManager>(
        "../data/simulation.sqlite", "../data/MFA-Seq_TBrucei_TREU927/");

    int transcription_period = 0;
    if (cmd_option_exists(argv, argv + argc, "--period"))
        transcription_period =
            atoi(argv[get_cmd_option(argv, argc, "--period")]);

    int origins_range = 0;
    if (cmd_option_exists(argv, argv + argc, "--constitutive"))
        origins_range =
            atoi(argv[get_cmd_option(argv, argc, "--constitutive")]);

    srand(time(NULL));
    omp_set_num_threads(40);

#pragma omp parallel for
    for (int i = 0; i < n_cells; i++)
    {
        SPhase s_phase(origins_range, n_resources, speed, timeout,
                       transcription_period, dormant, data, organism);
        s_phase.simulate(i);
    }
}
