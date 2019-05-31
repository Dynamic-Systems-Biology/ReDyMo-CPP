#!/usr/bin/env python3

""" This file is part of ReDyMo.

    Copyright (c) 2018  Gustavo Cayres and Marcelo Reis.

    ReDyMo is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.
    ReDyMo is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
    for more details.
    You should have received a copy of the GNU General Public License along
    with ReDyMo. If not, see <http://www.gnu.org/licenses/>.

"""

import sys
import os

aggregated_data = dict()

for result_file_name in next(os.walk(sys.argv[1]))[1]:

  result_path = '{}{}'.format(sys.argv[1], '/' + result_file_name +'/cell.txt')

  with open(result_path) as result_file:
    for line in result_file:
      if len(line) == 1:
        break

      N, speed, time, inter_dist = line.split('\t')[0:4]
      data = aggregated_data.get((int(N), int(speed)), [0, 0, 0, 0, 0])
      data[0] += float(time)
      data[1] += float(inter_dist)
      data[2] += float(time)**2
      data[3] += float(inter_dist)**2
      data[4] += 1
      aggregated_data[(int(N), int(speed))] = data

# print("F\tspeed\ttime_avg\ttime_sd\tinter_avg\tinter_sd" +\
# "\tmeasurements\t\n")

for key, value in aggregated_data.items():

  time_avg = value[0]/value[4]
  inter_avg = value[1]/value[4]

  time_sd = 0 if value[4] == 1 else ((value[2] / value[4] -\
  (value[0] / value[4]) ** 2) * (value[4] / (value[4] - 1))) ** (1 / 2)

  inter_sd = 0 if value[4] == 1 else ((value[3] / value[4] -\
  (value[1] / value[4]) ** 2) * (value[4] / (value[4] - 1))) ** (1 / 2)

  line = "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(key[0], key[1], time_avg,\
  time_sd, inter_avg, inter_sd, value[4])

  print(line, end="\t")

#-----------------------------------------------------------------------------#

