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

#load results from HelPopGen run
with open(FILE, 'r') as f:
    field_names = f.readline().rstrip().split(',')
    rows = [x for x in csv.DictReader(f, fieldnames=field_names)]

with open("jw_input_melpomene.csv", 'w') as f:
    f.writelines([','.join([row['D_N'],
                            row['NSsites'],
                            row['in_P_N'],
                            row['NSsites'],
                            row['D_S'],
                            row['Ssites'],
                            row['in_P_S'],
                            row['Ssites'],
                            row['ingroup_sequences']
                           ]) + '\n' for row in rows])

with open("jw_input_erato.csv", 'w') as f:
    f.writelines([','.join([row['D_N'],
                            row['NSsites'],
                            row['out_P_N'],
                            row['NSsites'],
                            row['D_S'],
                            row['Ssites'],
                            row['out_P_S'],
                            row['Ssites'],
                            row['outgroup_sequences']
                           ]) + '\n' for row in rows])

with open("jw_input_both.csv", 'w') as f:
    f.writelines([','.join([row['D_N'],
                            row['NSsites'],
                            str(int(row['out_P_N']) + int(row['in_P_N'])),
                            row['NSsites'],
                            row['D_S'],
                            row['Ssites'],
                            str(int(row['out_P_S']) + int(row['in_P_S'])),
                            row['Ssites'],
                            row['outgroup_sequences']
                           ]) + '\n' for row in rows])

