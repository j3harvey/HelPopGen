#! /usr/bin/env python2.7

from modelComparison import *
import numpy as np
import matplotlib.pyplot as plt

# load MKtest results
with open("mktest_out_mel_a", 'r') as f:
    r = f.readline().rstrip().split('\t')
    s = [float(x) for x in f.readline().rstrip().split('\t')]
    p = {k:v for (k, v) in zip(r, s)}

I = [i for (i, r) in enumerate(i_rows) if r['class'] != '0']
n = [i_rows[i]['name'] for i in I]
i_mel_a = {name: a for (name, a) in zip(n, [p['alpha' + str(i)] for i in I])}

# load Obbard et al 2009 results
with open("../data/s027_all.csv", 'r') as f:
    dmel_a = {r.split(',')[0].rstrip(): float(r.split(',')[5]) for r in f.readlines()[1:]}

# load orthology groups
with open("../data/curated_orthologues_both.csv", 'r') as f:
    groups = [x.rstrip().split(',') for x in f.readlines()[1:]]

dandha = list()
for g in groups:
    if g[-1].startswith("F"):
        dandha.append(g + [dmel_a[g[-1]],])
    elif g[-1] in n:
        dandha.append(g + [i_mel_a[g[-1]],])
    else:
        dandha.append(g + [np.nan,])

groups = [(int(g[0]), g[-1]) for g in groups]

i_hmel_a = list()
i_dmel_a = list()
for i in sorted(list(set([g[0] for g in groups]))):
    # get Hmel average a for group
    hmel_ids = [g[1] for g in groups
            if g[0] == i 
            and g[1].startswith("H")
            and g[1] in n]
    if len(hmel_ids) == 0:
        i_hmel_a.append(np.nan)
    else:
        i_hmel_a.append(sum([i_mel_a[ID] for ID in hmel_ids])/len(hmel_ids))
    # get Dmel average a for group
    dmel_ids = [g[1] for g in groups if g[0] == i and g[1].startswith("F")]
    i_dmel_a.append(sum([dmel_a[ID] for ID in dmel_ids])/len(dmel_ids))
    print i, i_hmel_a[-1], i_dmel_a[-1]
i_hmel_a = np.array(i_hmel_a)
i_dmel_a = np.array(i_dmel_a)

plt.scatter(i_hmel_a, i_dmel_a)
plt.xlabel("H. melpomene a")
plt.ylabel("D. melanogaster a")
plt.show()

