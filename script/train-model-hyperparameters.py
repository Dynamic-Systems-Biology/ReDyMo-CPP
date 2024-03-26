import logging
import signal
import os

import optuna
import subprocess
import numpy as np
import pandas as pd
from multiprocessing import Pool
from sklearn.metrics import mean_squared_error as mse
from sklearn.metrics import mean_absolute_error as mae


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
        start = int(start_s) if len(start_s) > 0 else int(range_values)
        end = int(end_v[0]) if (len(end_v) == 1 and len(start_s) > 0) else start
        seq_length = int(seq_length_v[0]) if len(seq_length_v) == 1 else 1

        # Decompress
        step = 1 if end > start else - 1
        for v in range(abs(end - start) + 1):
            current_val = start + v * step
            for _ in range(seq_length):
                yield current_val


def separate_training_set(num_chromosomes):
    chromosomes = np.arange(num_chromosomes) + 1
    training = np.random.choice(chromosomes, int(num_chromosomes * 0.7), replace=False)
    chromosomes = np.setdiff1d(chromosomes, training)
    validation = np.random.choice(chromosomes, int(num_chromosomes * 0.2), replace=False)
    testing = np.setdiff1d(chromosomes, validation)
    return training, validation, testing

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
    # Add 1 to avoid division by zero
    actual = actual + 1
    forecast = forecast + 1
    return (100 / len(actual)) * np.sum(np.abs(forecast - actual) / (np.abs(actual) + np.abs(forecast)))

# Calculate MSE and SMAPE for each chromosome
def calculate_errors_for_chrm(chrm, params):
    simulations = []

    # Read mfa-seq
    with open(f"../data/MFA-Seq_TcruziCLBrenerEsmeraldo-like/TcChr{chrm}-S.txt", 'r') as f:
        mfaseq = np.array(list(map(lambda x: float(x), list(f.readlines()))))

    # Read all simulations
    for sim in range(params['simulations']):
        with open(
                f"{params['name']}/round_{params['round']}_false_{params['replisomes']}_{params['replication_period']}/simulation_{sim}/TcChr{chrm}-S.cseq",
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
    return mse(mfaseq, norm_simulations), smape(mfaseq, norm_simulations), mae(mfaseq, norm_simulations)


def calculate_errors(params, chrm_list):
    p = Pool()  # create a pool of worker processes
    errors_for_each_chromosome = p.starmap(calculate_errors_for_chrm,
                                           [(chrm, params) for chrm in chrm_list])
    p.close()  # close the pool
    p.join()  # wait for all processes to finish

    log.info("Taking the mean of the errors of each chromosome")
    return np.mean(errors_for_each_chromosome, axis=0)


# TODO: gaussian curve parameters
def objective(trial):
    log.info(f"Starting trial {trial.number}")
    params = {
        'simulations': 1000,
        'timeout': 100_000_000,
        'speed': 1,
        'threads': 60,
        'organism': 'TcruziCLBrenerEsmeraldo-like',
        'num_chromosomes': 41,
        'probability': 0,
        'replisomes': trial.suggest_int('replisomes', 500, 1_500, 2),
        'replication_period': trial.suggest_int('period', 0, 500_000, 10),
        'constitutive_range': trial.suggest_int('constitutive_range', 0, 1_000_000, 10),
        'round': 0,
        'name': f"trial_{trial.number}",
        'training_chromosomes_set': TRAINING_CHROMOSOMES_SET,
        'validation_chromosomes_set': VALIDATION_CHROMOSOMES_SET,
        'test_chromosomes_set': TEST_CHROMOSOMES_SET
    }
    log.info(f"Starting simulation for trial {trial.number}")
    command_str = f"nice -n 20 ../simulator --cells {params['simulations']} --organism organism_placeholder --resources {params['replisomes']} --data-dir ../data --speed {params['speed']} --period {params['replication_period']} --constitutive {params['constitutive_range']} --timeout {params['timeout']} --threads {params['threads']} --name round_{params['round']} --summary --output {params['name']}"

    command_arr = str.split(command_str, ' ')
    command_arr[command_arr.index('organism_placeholder')] = params['organism']
    print(command_arr)

    sim_out = subprocess.run(command_arr, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    print(sim_out.stdout)
    print(sim_out.stderr)
    if sim_out.returncode != 0:
        log.error(f"Simulation for trial {trial.number} failed")
        trial.set_user_attr('error_exit_code', sim_out.returncode)
        raise optuna.exceptions.TrialPruned(f"Simulation for trial {trial.number} failed with return code {sim_out.returncode}")

    log.info(f"Ended simulation for trial {trial.number}")

    try:
        log.info("Calculating metrics for each chromosome (parallel)")
        training_mse, training_smape, training_mae = calculate_errors(params, params['training_chromosomes_set'])
        validation_mse, validation_smape, validation_mae = calculate_errors(params, params['validation_chromosomes_set'])
        trial.set_user_attr('training_mse', training_mse)
        trial.set_user_attr('validation_mse', validation_mse)
        trial.set_user_attr('training_smape', training_smape)
        trial.set_user_attr('validation_smape', validation_smape)
        trial.set_user_attr('training_mae', training_mae)
        trial.set_user_attr('validation_mae', validation_mae)


    except Exception as e:
        log.error(f"Error calculating error for trial {trial.number}")
        raise optuna.exceptions.TrialPruned(f"Error calculating error for trial {trial.number} " + str(e))

    return training_mse, training_mae, training_smape


def sigint_handler(signum, frame):
    log.warning(f"Stopping training because signal {signum} received")
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print(study.best_trials_dataframe().to_string())
    raise optuna.exceptions.TrialPruned(f"Trial interrupted with SIGINT or SIGTERM")


signal.signal(signal.SIGINT, sigint_handler)
signal.signal(signal.SIGTERM, sigint_handler)

# 70.6% of bases in genome
TRAINING_CHROMOSOMES_SET = [1, 2, 3, 4, 5, 6, 8, 9, 11, 12, 13, 15, 19, 20, 22, 23, 24, 25, 26, 27, 28, 30, 31, 33, 34, 36, 37, 38, 39, 40]
# 20.2% of bases in genome
VALIDATION_CHROMOSOMES_SET = [7, 10, 14, 17, 29, 35, 41]
# 9.1% of bases in genome
TEST_CHROMOSOMES_SET = [16, 18, 21, 32]

log.info(f"Starting trainer module")
log.info(f"Training chromosomes:   {np.sort(TRAINING_CHROMOSOMES_SET)}")
log.info(f"Validation chromosomes: {np.sort(VALIDATION_CHROMOSOMES_SET)}")
log.info(f"Test chromosomes:       {np.sort(TEST_CHROMOSOMES_SET)}")

OPTUNA_BACKEND_URL = os.environ['OPTUNA_BACKEND_URL']
study = optuna.create_study(directions=['minimize', 'minimize', 'minimize'], study_name='redymo-tcruzi-with-chipseq-multi-objective', storage=OPTUNA_BACKEND_URL,
                            load_if_exists=True)

study.optimize(objective, n_trials=4000)