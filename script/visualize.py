#! /usr/bin/python3

import sys
import subprocess
import os
from matplotlib import pyplot as plt

if len(sys.argv) < 2:
  print("File not specified")
  print("Usage: visualize.py FILE")
  print("Example: visualize.py path/to/file.txt.zst")
  exit(0)

compressed_file_name=os.path.abspath(sys.argv[1])
uncompressed_file_name=compressed_file_name[:-4]

subprocess.call('zstd -d -f ' + compressed_file_name + ' -o ' + uncompressed_file_name, shell=True)

with open(uncompressed_file_name, "r") as file:
  # step between points to be plotted. increase it to reduce time (will mess with y axis values)
  step = 1
  times=file.readlines()
  count = 0
  reduced = []

  # plot one every step
  for t in times:
    if count >= step:
      count = 0
      reduced.append(t)
    count = count + 1

  # configure plot ticks to not bunch up
  fig = plt.figure()
  ax = fig.add_subplot(1, 1, 1)
  max_yticks = 30
  yloc = plt.MaxNLocator(max_yticks)
  ax.yaxis.set_major_locator(yloc)

  plt.plot(reduced, '.g', markersize=0.3)
  plt.show()
