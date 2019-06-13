#!/usr/bin/python3
import os
import sys
import statistics as stat
import math
import multiprocessing as mp

# In this file a streak is a replication from a replication origin up to a
# collision with RNAP. Or if the streak has not ended yet, it is still referred
# to as a streak

testing_sample = [49484, 321, 409688, 413369, 77859, 952621, 540744, 32172, 486870, 50200, 48955, 842453, 1128739, 1165905, 159539, 216615, 26772, 5048, 502345, 1541822, 822259, 2272213, 49721, 236036, 186926, 246734, 1137875, 379739, 17575, 28574, 181984]
testing_metrics = [236036.0, 458424.4, 541215.9]

run_path = os.getcwd() + "/" + sys.argv[1]
file_model = "Tb927_{:02d}_v5.1.txt"

def find_streak(bases, collisions):
    distances = []
    streak = False
    start = -1
    for i in range(1, len(bases) - 2):
        back = i
        middle = i + 1
        front = i + 2
        if bases[back] - bases[middle] == 1 and bases[front] - bases[middle] == 1 and bases[front] - bases[back] == 0:
            streak = True
            start = middle
        elif bases[middle] - bases[back] == 1 and bases[front] - bases[middle] > 1:
            if streak:
                streak = False
                collisions += 1
                distances.append(middle - start)
                start = -1
            else:
                #print("error, not in a streak" + "=="*30)
                pass
    return distances, collisions

def process(chrm):
    collisions = 0
    chrm_distances = []
    for sim in os.listdir(run_path):
        sim_file = open(run_path + "/" + sim + "/" + file_model.format(chrm), 'r')
        file_content = sim_file.readlines()
        for i,_ in enumerate(file_content):
            file_content[i] = int(file_content[i])
        ret_distances, collisions = find_streak(file_content, collisions)
        chrm_distances += ret_distances
        file_content = file_content[::-1]
        ret_distances, collisions = find_streak(file_content, collisions)
        chrm_distances += ret_distances
        if len(chrm_distances) == 0:
            chrm_distances = []
    # print(f"Collisions: {collisions}")
    return chrm_distances


def main():
    testing = False
    if sys.argv[1] == "--test":
        global run_path
        run_path = os.getcwd() + "/../data/test_output"
        testing = True

    n_procs = 11
    pool = mp.Pool(n_procs)
    chrms = range(1, 12)

    overall_distances = pool.map(process, chrms)

    # inefficient but we only do this for 11 elements
    overall_distances = sum(overall_distances, [])

    median = stat.median(overall_distances)
    mean = stat.mean(overall_distances)
    stdev = stat.stdev(overall_distances)

    if testing:
        if overall_distances != testing_sample:
            print("Detection of streaks failed!")
            exit()
        if not math.isclose(median, testing_metrics[0], rel_tol=0.1):
            print("median failed!")
            exit()
        if not math.isclose(mean, testing_metrics[1], rel_tol=0.1):
            print("mean failed!")
            exit()
        if not math.isclose(stdev, testing_metrics[2], rel_tol=0.1):
            print("stdev failed!")
            exit()
        print("Tests passed.....OK")

    print(f"{median:.01f}\t{mean:.01f}\t{stdev:.01f}")

main()
