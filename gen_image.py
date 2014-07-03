#! /usr/bin/env python3

import os
import shutil
import sys
import random
import time

solution_logs_dir = 'solution_logs'
seed = int(sys.argv[1])

def make_command(seed, filename):
    solution_logs_file = solution_logs_dir + '/' + str(seed) + '/' + filename
    return "java CollageMakerVis -exec 'cat {}' -novis -seed {} -image image/{} {}.png".format(solution_logs_file, seed, seed, filename)

for filename in sorted(os.listdir(solution_logs_dir + '/' + str(seed))):
    command = make_command(seed, filename)
    print(command)
    os.popen(command).read()
