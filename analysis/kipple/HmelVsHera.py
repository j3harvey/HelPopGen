#! /usr/bin/env python2.7

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("file", 
                    help="Output file from HelPopGen", 
                    nargs='?', 
                    type=str, 
                   )
args = parser.parse_args()
FILE = args.file

import numpy as np
import matplotlib.pyplot as plt
import csv

#load results from HelPopGen run
with open(FILE, 'r') as f:
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

# compare gene length to number of adaptive substitutions
mel_a = (D_N - D_S*mel_P_N/mel_P_S)/NSsites
era_a = (D_N - D_S*era_P_N/era_P_S)/NSsites

# compare alpha values for erato and melpomene
alpha_era = [1 - D_S[i]*era_P_N[i]/(D_N[i]*era_P_S[i])
           if all([D_N[i] >= 5,
                   D_S[i] >= 5,
                   mel_P_N[i] >= 5,
                   mel_P_S[i] >= 5,
                   era_P_N[i] >= 5,
                   era_P_S[i] >= 5,]) else np.nan for i in range(len(rows))]

alpha_mel = [1 - D_S[i]*mel_P_N[i]/(D_N[i]*mel_P_S[i])
           if all([D_N[i] >= 5,
                   D_S[i] >= 5,
                   mel_P_N[i] >= 5,
                   mel_P_S[i] >= 5,
                   era_P_N[i] >= 5,
                   era_P_S[i] >= 5,]) else np.nan for i in range(len(rows))]

# load list of immune loci
with open("../data/ImmunityGenes.csv", 'r') as f:
    immune_genes = [row.split(',')[0].rstrip() for row in f.readlines()]

immune_mel_P_N = np.array([a for (i, a) in enumerate(mel_P_N) if loci[i] in immune_genes])
immune_mel_P_S = np.array([a for (i, a) in enumerate(mel_P_S) if loci[i] in immune_genes])
immune_era_P_N = np.array([a for (i, a) in enumerate(era_P_N) if loci[i] in immune_genes])
immune_era_P_S = np.array([a for (i, a) in enumerate(era_P_S) if loci[i] in immune_genes])
immune_D_N = np.array([a for (i, a) in enumerate(D_N) if loci[i] in immune_genes])
immune_D_S = np.array([a for (i, a) in enumerate(D_S) if loci[i] in immune_genes])

immune_alpha_mel = 1 - immune_D_S*immune_mel_P_N/(immune_D_N*immune_mel_P_S)
immune_alpha_era = 1 - immune_D_S*immune_era_P_N/(immune_D_N*immune_era_P_S)

#i_mel_a = np.array([a for (i, a) in enumerate(mel_a) if loci[i] in immune_genes])
#i_era_a = np.array([a for (i, a) in enumerate(era_a) if loci[i] in immune_genes])

i_mel_a = [(D_N[i] - D_S[i]*mel_P_N[i]/mel_P_S[i])/NSsites[i]
           if all([D_N[i] >= 3,
                   D_S[i] >= 3,
                   mel_P_N[i] >= 3,
                   mel_P_S[i] >= 3,
                   era_P_N[i] >= 3,
                   era_P_S[i] >= 3,]) 
           and loci[i] in immune_genes else np.nan for i in range(len(rows))]

i_era_a = [(D_N[i] - D_S[i]*era_P_N[i]/era_P_S[i])/NSsites[i]
           if all([D_N[i] >= 3,
                   D_S[i] >= 3,
                   era_P_N[i] >= 3,
                   era_P_S[i] >= 3,
                   mel_P_N[i] >= 3,
                   mel_P_S[i] >= 3,]) 
           and loci[i] in immune_genes else np.nan for i in range(len(rows))]

print len(rows)

plt.scatter(i_era_a, i_mel_a)
plt.show()
