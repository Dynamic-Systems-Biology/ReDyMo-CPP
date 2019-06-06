#!/usr/bin/python3
import os
import sys
import statistics as stat
import multiprocessing as mp


run_path = os.getcwd() + "/" + sys.argv[1]
file_model = "Tb927_{:02d}_v5.1.txt"
collisions = 0

def find_streak(bases):
    global collisions
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
    return distances

def process(chrm):
    chrm_distances = []
    for sim in os.listdir(run_path):
        sim_file = open(run_path + "/" + sim + "/" + file_model.format(chrm), 'r')
        file_content = sim_file.readlines()
        for i,_ in enumerate(file_content):
            file_content[i] = int(file_content[i])
        chrm_distances += find_streak(file_content)
        file_content = file_content[::-1]
        chrm_distances += find_streak(file_content)
        # print(f"Collisions: {collisions}")
        collisions = 0
        if len(chrm_distances) == 0:
            chrm_distances = [0]
    return chrm_distances


def main():
    n_procs = 11
    pool = mp.Pool(n_procs)
    global collisions
    chrms = range(1, 3)
    overall_distances = pool.map(process, chrms)
    # inefficient but we only do this for 11 elements
    overall_distances = sum(overall_distances, [])
    # print(chrm_distances)
    # print(f"Median for Chromosome {chrm:02d}: {stat.median(chrm_distances):010,.01f}")

    # print(f"Overall median distance between origin and conflict: {stat.median(overall_distances):,.01f} bp ")
    print(f"{stat.median(overall_distances):.01f}\t{stat.mean(overall_distances):.01f}\t{stat.stdev(overall_distances):.01f}")


main()
