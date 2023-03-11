import signal

import optuna
import subprocess
import numpy as np
import pandas as pd
from multiprocessing import Pool


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
    return (array - array.min()) / (array.max() - array.min())


def linear_interpolation(a, b, n):
    return [a + (b - a) * i / float(n) for i in range(n)]


def interpolate(array, new_size):
    old_size = len(array)
    new_array = np.linspace(0, old_size - 1, new_size)
    new_array = np.interp(new_array, np.arange(old_size), array)
    return new_array


def mse(a, b):
    return sum((a[i] - b[i]) ** 2 for i in range(len(a))) / len(a)


def calculate_error_for_chrm(chrm, params):
    simulations = []

    # Read mfa-seq
    with open(f"../data/MFA-Seq_TcruziCLBrenerEsmeraldo-like/TcChr{chrm}-S.txt", 'r') as f:
        mfaseq = np.array(list(map(lambda x: float(x), list(f.readlines()))))

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
    mfaseq = normalize(mfaseq)
    # Interpolate mfa-seq
    mfaseq = interpolate(np.array(mfaseq), len(simulations[0]))
    print("mfaseq", mfaseq.shape)

    # Calculate average of simulations
    simulations_average = np.mean(simulations, axis=0)
    print("simulations_average", simulations_average.shape)

    # Normalize simulation
    norm_simulations = normalize(simulations_average)
    print("norm_simulations", norm_simulations.shape)

    # Calculate error
    return mse(mfaseq, norm_simulations)


def calculate_error(params):
    p = Pool()  # create a pool of worker processes
    mse_for_each_chromosome = p.starmap(calculate_error_for_chrm,
                                        [(chrm, params) for chrm in params['training_chromosomes_set']])
    p.close()  # close the pool
    p.join()  # wait for all processes to finish
    print(mse_for_each_chromosome)
    return np.mean(mse_for_each_chromosome)


# TODO: gaussian curve parameters
def objective(trial):
    print(f"Starting trial {trial.number}")
    params = {
        'simulations': 1000,
        'timeout': 100_000_000,
        'speed': 1,
        'threads': 70,
        'organism': 'TcruziCLBrenerEsmeraldo-like',
        'num_chromosomes': 41,
        'probability': 0,
        'replisomes': trial.suggest_int('replisomes', 2, 1_002, 10),
        'replication_period': trial.suggest_int('period', 0, 1_000_000, 100),
        'round': 0,
        'name': f"trial_{trial.number}",
        'training_chromosomes_set': TRAINING_CHROMOSOMES_SET
    }
    command_str = f"nice -n 20 ../simulator --cells {params['simulations']} --organism {params['organism']} --resources {params['replisomes']} --data-dir ../data --speed {params['speed']} --period {params['replication_period']} --timeout {params['timeout']} --threads {params['threads']} --name round_{params['round']} --summary --output {params['name']}"

    command_arr = str.split(command_str, ' ')

    sim_out = subprocess.run(command_arr, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    print(sim_out.stdout)
    print(sim_out.stderr)

    print("End of execution for trial {trial.number}")

    return calculate_error(params)


def sigint_handler(signum, frame):
    print(f"Stopping training because signal {signum} received")
    study.stop()
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print(study.trials_dataframe().to_string())
    print(study.best_params)


signal.signal(signal.SIGINT, sigint_handler)
signal.signal(signal.SIGTERM, sigint_handler)

TRAINING_CHROMOSOMES_SET = separate_training_set(41)
TRAINING_DB = 'sqlite:///opt/redymo/train.sqlite'
study = optuna.create_study(direction='minimize', study_name='redymo-no-chipseq', storage=TRAINING_DB,
                            load_if_exists=True)

study.optimize(objective, n_trials=100)

with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    print(study.trials_dataframe().to_string())
print(study.best_params)
