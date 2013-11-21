#! /usr/bin/env python2.7

'''
This script converts the output from HelPopGen into a suitable
input file for John Welch's software package MKtest

'''

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("file", 
                    help="Output file from HelPopGen", 
                    nargs='?', 
                    type=str, 
                   )
args = parser.parse_args()
FILE = args.file

import csv
from controlGenes import getControlGenes, mapped

# immunity genes by class
classes = {"humoral_recognition": '1',
           "cellular_recognition": '2',
           "RNAi": '3',
           "AMP": '4',
           "signalling": '5',
          }

# load results from HelPopGen run
with open("../results/results_thres0.1_150615052013.csv", 'r') as f:
    field_names = f.readline().rstrip().split(',')
    rows = [x for x in csv.DictReader(f, fieldnames=field_names)]

# load list of immunity genes
with open("../data/ImmunityGenes2.csv", 'r') as f:
    field_names = f.readline().rstrip().split(',')
    immunity_genes = [x for x in csv.DictReader(f, fieldnames=field_names)]

mapped_immunity_genes = mapped([x['Hmel_ID'] for x in immunity_genes])

control_genes = getControlGenes(mapped_immunity_genes)

def writeMKtestInputFile(FILE, lines):
    with open(FILE, 'w') as f:
        f.writelines([','.join([line['D_N'],
                                line['NSsites'],
                                line['in_P_N'],
                                line['NSsites'],
                                line['D_S'],
                                line['Ssites'],
                                line['in_P_S'],
                                line['Ssites'],
                                line['ingroup_sequences'],
                                '0', # chr
                                line[cl], # class
                               ]) + '\n' for line in lines])

# genome-wide alpha
with open("mk_genome.csv", 'w') as f:
    f.writelines([','.join([r['D_N'],
                            r['NSsites'],
                            r['in_P_N'],
                            r['NSsites'],
                            r['D_S'],
                            r['Ssites'],
                            r['in_P_S'],
                            r['Ssites'],
                            r['ingroup_sequences'],
                           ]) + '\n' for r in rows])
# Command line:
# Mktest -A mk_genome.csv

i_rows = [r for r in rows 
          if r['name'] in mapped_immunity_genes
          or r['name'] in control_genes]
for r in i_rows:
    r['class'] = '0'
for r in i_rows:
    for l in immunity_genes:
        if l['Hmel_ID'] == r['name']:
            if l['Cell_Hum'] in classes.keys():
                r['class'] = classes[l['Cell_Hum']]
            else:
                r['class'] = 999
i_rows = [x for x in i_rows if x['class'] != 999]
 
# single value of alpha for immunity and control
with open("mk_all.csv", 'w') as f:
    f.writelines([','.join([r['D_N'],
                            r['NSsites'],
                            r['in_P_N'],
                            r['NSsites'],
                            r['D_S'],
                            r['Ssites'],
                            r['in_P_S'],
                            r['Ssites'],
                            r['ingroup_sequences'],
                           ]) + '\n' for r in i_rows])
# Command line:
# Mktest -A mk_all.csv

# immunity vs control
with open("mk_immune_vs_control.csv", 'w') as f:
    f.writelines([','.join([r['D_N'],
                            r['NSsites'],
                            r['in_P_N'],
                            r['NSsites'],
                            r['D_S'],
                            r['Ssites'],
                            r['in_P_S'],
                            r['Ssites'],
                            r['ingroup_sequences'],
                            '0',
                            '1' if r['name'] in [x['Hmel_ID'] for x in immunity_genes] else '0',
                           ]) + '\n' for r in i_rows])
# Command line:
# Mktest -a4 -A mk_immune_vs_control.csv

with open("mk_classes.csv", 'w') as f:
    f.writelines([','.join([r['D_N'],
                            r['NSsites'],
                            r['in_P_N'],
                            r['NSsites'],
                            r['D_S'],
                            r['Ssites'],
                            r['in_P_S'],
                            r['Ssites'],
                            r['ingroup_sequences'],
                            '0',
                            r['class'],
                           ]) + '\n' for r in i_rows])
# Command line:
# Mktest -a4 -A mk_classes.csv

# individual alpha at each locus
with open("mk_loci.csv", 'w') as f:
    f.writelines([','.join([r['D_N'],
                            r['NSsites'],
                            r['in_P_N'],
                            r['NSsites'],
                            r['D_S'],
                            r['Ssites'],
                            r['in_P_S'],
                            r['Ssites'],
                            r['ingroup_sequences'],
                            '0',
                            r['class'],
                           ]) + '\n' for r in i_rows])
# Command line:
# Mktest -A mk_loci.csv

