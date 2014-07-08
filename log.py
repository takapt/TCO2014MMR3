#! /usr/bin/env python

import os
import shutil
import sys
import random
import time

HOME = os.environ['HOME']
MM = os.path.join(HOME, 'mm')
exe = 'a.out' if len(sys.argv) < 2 else sys.argv[1]

result_name = None
if len(sys.argv) >= 3:
    result_name = sys.argv[2]

exe_path = os.path.join(MM, exe)
copied_exe_path = os.path.join(MM, 'copied_' + str(random.randint(0, 10**5)))
shutil.copy(exe_path, copied_exe_path)

def make_command(seed):
    command = "java CollageMakerVis -exec '{}' -novis -seed {}".format(copied_exe_path, seed)
    if result_name:
        command += ' -image results/{} {:06d}.png'.format(result_name, seed)
    return command

def get_score(seed):
    c = make_command(seed)

    start = time.time()
    try:
        output = os.popen(c).read()
    except:
        raise "exit"

    score = float(output.split()[-1])
    exe_time = time.time() - start

    return {'seed': seed, 'score': score, 'time': exe_time}

def single(seeds):
    for seed in seeds:
        result = get_score(seed)
        print('{:4d} {:3.3f} {:.3f}'.format(seed, result['score'], result['time']))
        sys.stdout.flush()

def multi(seeds):
    from multiprocessing import Pool
    pool = Pool()
    results = pool.map(get_score, seeds)
    for result in results:
        seed = result['seed']
        print('{:4d} {:3.3f} {:.3f}'.format(seed, result['score'], result['time']))
        sys.stdout.flush()

try:
    single(range(1, 20))
#     multi(range(1, 20))
finally:
    os.remove(copied_exe_path)
