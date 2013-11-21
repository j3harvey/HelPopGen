#! /usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import csv

# load list of immune loci
with open("../data/ImmunityGenes.csv", 'r') as f:
    immune_genes = [row.split(',')[0].rstrip() for row in f.readlines()]

#load results from HelPopGen run
with open("../results/results_thres0.2_150615052013.csv", 'r') as f:
    field_names = f.readline().rstrip().split(',')
    rows = [x for x in csv.DictReader(f, fieldnames=field_names)]

loci = [r['name'] for r in rows]
mel_sequences = np.array([int(r['ingroup_sequences']) for r in rows])
era_sequences = np.array([int(r['outgroup_sequences']) for r in rows])
sequence_length = np.array([float(r['usable_sequence_length']) for r in rows])
NSsites = np.array([float(r['NSsites']) for r in rows])
Ssites = np.array([float(r['Ssites']) for r in rows])
mel_P_N = np.array([float(r["in_P_N"]) for r in rows])
mel_P_S = np.array([float(r["in_P_S"]) for r in rows])
era_P_N = np.array([float(r["out_P_N"]) for r in rows])
era_P_S = np.array([float(r["out_P_S"]) for r in rows])
D_N = np.array([float(r["D_N"]) for r in rows])
D_S = np.array([float(r["D_S"]) for r in rows])

with open("mktest_out_mel", 'r') as f:
    mel_a = f.readlines()[-1].rstrip().split('\t')[4:]

mel_a = np.array([float(x) for x in mel_a])

with open("mktest_out_era", 'r') as f:
    era_a = f.readlines()[-1].rstrip().split('\t')[4:]

era_a = np.array([float(x) for x in era_a])

i_mel_a = {l: a for (l, a) in zip(loci, mel_a) if l in immune_genes}
i_era_a = {l: a for (l, a) in zip(loci, era_a) if l in immune_genes}

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
    elif g[-1] in loci:
        dandha.append(g + [i_mel_a[g[-1]],])
    else:
        dandha.append(g + [np.nan,])

groups = [(int(g[0]), g[-1]) for g in groups]

i_hera_a = list()
i_hmel_a = list()
i_dmel_a = list()
for i in sorted(list(set([g[0] for g in groups]))):
    # get Hera average a for group
    hera_ids = [g[1] for g in groups
            if g[0] == i 
            and g[1].startswith("H")
            and g[1] in loci]
    if len(hera_ids) == 0:
        i_hera_a.append(np.nan)
    else:
        i_hera_a.append(sum([i_era_a[ID] for ID in hera_ids])/len(hera_ids))
    # get Hmel average a for group
    hmel_ids = [g[1] for g in groups
            if g[0] == i 
            and g[1].startswith("H")
            and g[1] in loci]
    if len(hmel_ids) == 0:
        i_hmel_a.append(np.nan)
    else:
        i_hmel_a.append(sum([i_mel_a[ID] for ID in hmel_ids])/len(hmel_ids))
    # get Dmel average a for group
    dmel_ids = [g[1] for g in groups if g[0] == i and g[1].startswith("F")]
    i_dmel_a.append(sum([dmel_a[ID] for ID in dmel_ids])/len(dmel_ids))
    print i, i_hmel_a[-1], i_dmel_a[-1]
i_hera_a = np.array(i_hmel_a)
i_hmel_a = np.array(i_hmel_a)
i_dmel_a = np.array(i_dmel_a)

plt.scatter(i_hmel_a, i_dmel_a)
plt.show()

