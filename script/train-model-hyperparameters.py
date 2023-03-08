import optuna
import subprocess
import numpy as np


# From https://github.com/seijihariki/redymo-tcruzi-analysis/blob/master/final_scripts/compressor.py
def decompress(input_it):
    'Decompresses input lines to output lines (stream compressor)'

    it = iter(input_it)

    for line in it:
        # Parse compressed file

        # Split line
        range_values, *seq_length_v = line.split('x')
        start_s, *end_v = range_values.split('-')

        # Parse values
        start = int(start_s)
        end = int(end_v[0]) if len(end_v) == 1 else start
        seq_length = int(seq_length_v[0]) if len(seq_length_v) == 1 else 1

        # Decompress
        step = 1 if end > start else - 1
        for v in range(abs(end - start) + 1):
            current_val = start + v * step
            for _ in range(seq_length):
                yield current_val


def separate_training_set(num_chromosomes):
    chromosomes = np.arange(num_chromosomes)
    return np.random.choice(chromosomes, int(len(chromosomes) * 0.7), replace=False)


def normalize(array):
    return (array - np.min(array)) / (np.max(array) - np.min(array))


def linear_interpolation(a, b, n):
    return [a + (b - a) * i / float(n) for i in range(n)]


def interpolate(array, new_size):
    old_size = len(array)
    new_array = []
    for i in range(old_size - 1):
        expanded = linear_interpolation(array[i], array[i + 1], new_size // old_size)
        print(expanded)
        new_array += expanded
    return np.array(new_array)


def mse(a, b):
    return sum((a[i] - b[i]) ** 2 for i in range(len(a))) / len(a)


def calculate_error(params):
    mse_for_each_chromosome = []
    for chrm in params['training_chromosomes_set']:
        simulations = []

        # Read mfa-seq
        with open(f"../data/MFA-Seq_TcruziCLBrenerEsmeraldo-like/TcChr{chrm}-S.txt", 'r') as f:
            mfaseq = np.array(list(f.readlines()))

        print("mfa-seq", mfaseq.shape)

        # Read all simulations
        for sim in range(params['simulations']):
            with open(
                    f"{params['name']}/round_{params['round']}_false_{params['replisomes']}_{params['replication_period']}/simulation_{sim}/TcChr{chrm}-S.cseq",
                    'r') as f:
                simulations.append(np.array(list(decompress(f.readlines()))))

        simulations = np.array(simulations)
        print("simulations", simulations.shape)

        # Normalize mfa-seq
        mfaseq = normalize(np.array(mfaseq))
        # Interpolate mfa-seq
        mfaseq = interpolate(np.array(mfaseq), len(simulations[0]))
        print("mfaseq", mfaseq.shape)
        import matplotlib.pyplot as plt
        plt.plot(mfaseq)
        plt.interactive(False)
        plt.show()

        # Calculate average of simulations
        simulations_average = np.mean(simulations, axis=0)
        print("simulations_average", simulations_average.shape)

        # Normalize simulation
        norm_simulations = normalize(simulations_average)
        print("norm_simulations", norm_simulations.shape)

        # Calculate error
        mse_for_each_chromosome.append(mse(mfaseq, norm_simulations))
    print(mse_for_each_chromosome)
    return np.mean(mse_for_each_chromosome)


# TODO: gaussian curve parameters
def objective(trial):
    params = {
        'simulations': 10,
        'timeout': 100_000_000,
        'speed': 1,
        'threads': 6,
        'organism': 'TcruziCLBrenerEsmeraldo-like',
        'num_chromosomes': 41,
        'probability': 0,
        'replisomes': trial.suggest_int('replisomes', 2, 1_000, 10),
        'replication_period': trial.suggest_int('period', 0, 1_000_000, 100),
        'round': 0,
        'name': f"trial_{trial.number}",
        'training_chromosomes_set': TRAINING_CHROMOSOMES_SET
    }
    command_str = f"nice -n 20 ../build/simulator --cells {params['simulations']} --organism {params['organism']} --resources {params['replisomes']} --data-dir ../data --speed {params['speed']} --period {params['replication_period']} --timeout {params['timeout']} --threads {params['threads']} --name round_{params['round']} --summary --output {params['name']}"

    command_arr = str.split(command_str, ' ')
    ret = subprocess.run(command_arr, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    return calculate_error(params)

TRAINING_CHROMOSOMES_SET = separate_training_set(41)
study = optuna.create_study(direction='minimize')
study.optimize(objective, n_trials=2)
print(study.trials_dataframe())

print(study.best_params)
