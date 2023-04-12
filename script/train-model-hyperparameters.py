import logging
import signal

import optuna
import subprocess
import numpy as np
import pandas as pd
from multiprocessing import Pool

# Global logger config
logging.basicConfig(level=logging.INFO, format='[%(levelname).1s %(asctime)s] %(message)s')
log = logging.getLogger(__name__)


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
    chromosomes = np.arange(num_chromosomes) + 1
    training = np.random.choice(chromosomes, int(len(chromosomes) * 0.7), replace=False)
    validation = np.setdiff1d(chromosomes, training)
    return training, validation


def normalize(array):
    return (array - array.min()) / (array.max() - array.min())


def linear_interpolation(a, b, n):
    return [a + (b - a) * i / float(n) for i in range(n)]


def interpolate(array, new_size):
    old_size = len(array)
    new_array = np.linspace(0, old_size - 1, new_size)
    new_array = np.interp(new_array, np.arange(old_size), array)
    return new_array


# Calculate SMAPE (Symmetric Mean Absolute Percentage Error)
def smape(actual, forecast):
    return (100 / len(actual)) * np.sum(np.abs(forecast - actual) / (np.abs(actual) + np.abs(forecast)))


def mse(a, b):
    return sum((a[i] - b[i]) ** 2 for i in range(len(a))) / len(a)


# Calculate MSE and SMAPE for each chromosome
def calculate_errors_for_chrm(chrm, params):
    simulations = []

    # Read mfa-seq
    with open(f"../data/MFA-Seq_{params['organism']}/Tb927_{chrm:02d}_v5.1.txt", 'r') as f:
        mfaseq = np.array(list(map(lambda x: float(x), list(f.readlines()))))

    # Read all simulations
    for sim in range(params['simulations']):
        with open(
                f"{params['name']}/round_{params['round']}_false_{params['replisomes']}_{params['replication_period']}/simulation_{sim}/Tb927_{chrm:02d}_v5.1.cseq",
                'r') as f:
            simulations.append(np.array(list(decompress(f.readlines()))))

    simulations = np.array(simulations)

    # Normalize mfa-seq
    mfaseq = normalize(mfaseq)

    # The MFA-Seq data means the inverse of what the simulated data means. High
    # MFA-Seq means sooner in the S-phase, lower simulation data means sooner. So,
    # the higher the value in the MFA-Seq data, the lower the value in the simulated
    # data should be. Since we normalized it, we cant invert it with the value 1.
    mfaseq = 1 - mfaseq

    # Interpolate mfa-seq
    mfaseq = interpolate(np.array(mfaseq), len(simulations[0]))

    # Calculate average of simulations
    simulations_average = np.mean(simulations, axis=0)

    # Normalize simulation
    norm_simulations = normalize(simulations_average)

    # Calculate error
    return mse(mfaseq, norm_simulations), smape(mfaseq, norm_simulations)


def calculate_errors(params, chrm_list):
    p = Pool()  # create a pool of worker processes
    log.info("Calculating error for each chromosome (parallel)")
    errors_for_each_chromosome = p.starmap(calculate_errors_for_chrm,
                                           [(chrm, params) for chrm in chrm_list])
    p.close()  # close the pool
    p.join()  # wait for all processes to finish

    log.info("Taking the mean of the errors of each chromosome")
    return np.mean(errors_for_each_chromosome)


# TODO: gaussian curve parameters
def objective(trial):
    log.info(f"Starting trial {trial.number}")
    params = {
        'simulations': 500,
        'timeout': 100_000_000,
        'speed': 1,
        'threads': 60,
        'organism': 'Trypanosoma brucei brucei TREU927',
        'num_chromosomes': 11,
        'probability': 0,
        'replisomes': trial.suggest_int('replisomes', 2, 1_002, 2),
        'replication_period': trial.suggest_int('period', 0, 1_000_000, 1),
        'round': 0,
        'name': f"trial_{trial.number}",
        'training_chromosomes_set': TRAINING_CHROMOSOMES_SET,
        'validation_chromosomes_set': VALIDATION_CHROMOSOMES_SET
    }
    log.info(f"Starting simulation for trial {trial.number}")
    command_str = f"nice -n 20 ../simulator --cells {params['simulations']} --organism {params['organism']} --resources {params['replisomes']} --data-dir ../data --speed {params['speed']} --period {params['replication_period']} --timeout {params['timeout']} --threads {params['threads']} --name round_{params['round']} --summary --output {params['name']}"

    command_arr = str.split(command_str, ' ')

    sim_out = subprocess.run(command_arr, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    print(sim_out.stdout)
    print(sim_out.stderr)
    if sim_out.returncode != 0:
        log.error(f"Simulation for trial {trial.number} failed")
        raise optuna.exceptions.TrialPruned(f"Simulation for trial {trial.number} failed")

    log.info(f"Ended simulation for trial {trial.number}")

    try:
        training_mse, training_smape = calculate_errors(params, params['training_chromosomes_set'])
        validation_mse, validation_smape = calculate_errors(params, params['validation_chromosomes_set'])
        trial.set_user_attr('training_mse', training_mse)
        trial.set_user_attr('validation_mse', validation_mse)
        trial.set_user_attr('training_smape', training_smape)
        trial.set_user_attr('validation_smape', validation_smape)


    except Exception as e:
        log.error(f"Error calculating error for trial {trial.number}")
        raise optuna.exceptions.TrialPruned(f"Error calculating error for trial {trial.number}")

    return training_mse


def sigint_handler(signum, frame):
    log.warning(f"Stopping training because signal {signum} received")
    study.stop()
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print(study.trials_dataframe().to_string())
    print(study.best_params)


signal.signal(signal.SIGINT, sigint_handler)
signal.signal(signal.SIGTERM, sigint_handler)

TRAINING_CHROMOSOMES_SET = [1, 5, 6, 7, 8, 9, 10, 11]
VALIDATION_CHROMOSOMES_SET = [2, 3, 4]

log.info(f"Starting trainer module")
log.info(f"Training chromosomes:   {np.sort(TRAINING_CHROMOSOMES_SET)}")
log.info(f"Validation chromosomes: {np.sort(VALIDATION_CHROMOSOMES_SET)}")

RDB_BACKEND_URL = 'sqlite:///train.sqlite'
study = optuna.create_study(direction='minimize', study_name='redymo-t-brucei', storage=RDB_BACKEND_URL,
                            load_if_exists=True)

study.optimize(objective, n_trials=100)
