#! /usr/bin/env python2.7

with open("all_orthomcl.out.parsed", 'r') as f:
    lines = f.readlines()

DmelNames = [l.split()[1].split('-')[0] for l in lines]
print len(set(DmelNames))
