#! /usr/bin/env python3

import sys

def get_scores(filename):
    scores = {}
    for line in open(filename):
        s = line.split()
        seed = int(s[0])
        score = float(s[1])
        scores[seed] = score
    return scores

a, b = sys.argv[1], sys.argv[2]
p, q = get_scores(a), get_scores(b)
seeds = list(set(p) & set(q))
sum_ratio = 0
for seed in seeds:
    x, y = p[seed], q[seed]
    ratio = x / y
    sum_ratio += ratio

    print('{:>5} {:>16} {:>16} {:>7.3f}'.format(seed, x, y, ratio))

total_ratio = sum_ratio / len(seeds)
print('total_ratio: {}'.format(total_ratio))
