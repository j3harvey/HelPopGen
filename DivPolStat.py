#! /usr/bin/env python2.7

import Bio.SeqIO
from Bio.Data import CodonTable
import csv
from collections import Counter
import itertools


###
# 
# HELPER FUNCTIONS AND CONSTANTS
# 

default_codon_table = CodonTable.standard_dna_table

'''
D is a dictionary of distances between codons.

For codons c and c, D[(c,d)] is a tuple:
the first element is the number of non-synonymous mutations needed to change
c to d, and the second element if the number of synonymous mutations needed.

'''
D = dict()
_FT=default_codon_table.forward_table
for c in _FT.keys():
    for i in range(3):
        for b in "ACGT":
            d = ''.join([c[j] if j != i else b for j in range(3)])
            if d == c:
                D[(c,d)] = (0,0)
            elif d not in _FT.keys():
                continue
            else:
                D[(c,d)] = (0,1) if _FT[c] == _FT[d] else (1,0)
for c1 in _FT.keys():
   for c2 in _FT.keys():
        nequal = sum([c1[i] != c2[i] for i in range(3)])
        if nequal < 2:
            continue
        elif nequal == 2:
            D[(c1,c2)] = min(
              [(D[(c1,d)][0] + D[(d,c2)][0], D[(c1,d)][1] + D[(d,c2)][1]) 
              for d in _FT.keys() if (c1,d) in D.keys() and (d,c2) in D.keys()]
              )
        else:
            continue
for c1 in _FT.keys():
    for c2 in _FT.keys():
        nequal = sum([c1[i] != c2[i] for i in range(3)])
        if nequal == 3:
            D[(c1,c2)] = min(
              [(D[(c1,d)][0] + D[(d,c2)][0], D[(c1,d)][1] + D[(d,c2)][1]) 
              for d in _FT.keys() if (c1,d) in D.keys() and (d,c2) in D.keys()]
              )

def countMutations(site, ct=default_codon_table):
    '''
    Uses a heuristic method (Kruskal's algorithm) to find a set of mutations 
    that can explain the polymorphism observed at the site.

    The return value is a tuple: the first element is the number of
    non-synonymous mutations, and the second element is the number of
    sysnonymous mutations.

    This is an example of a Steiner Tree Problem, which is known to be
    NP-hard.

    '''
    E = set([tuple(sorted((c,d))) for c in site for d in site if d != c])
    E = sorted(E, key=lambda e: D[e])
    T = list()
    sets = {c : set([c,]) for c in site}
    if len(E) == 0:
        return (0,0)
    for (u,v) in E:
        if sets[u] != sets[v]:
            for (key,val) in sets.items():
                if u in val or v in val:
                    sets[key] = sets[u] | sets[v]
            T.append((u,v))
    return (
        sum([D[e][0] for e in set([tuple(sorted(t)) for t in T])]), 
        sum([D[e][1] for e in set([tuple(sorted(t)) for t in T])])
        )

###
# 
# MUTATIONAL OPPORTTUNITY STATISTICS
# 

def meanNonSynPos(record, codonTable=default_codon_table):
    '''Calculates the number of non-synonymous positions in an alignment.

    Calculates the mean number of non-synonymous sites in a set of aligned 
    sequences. The return value is the sum of the siteMeanNonSynPos for
    each site in the alignment.

    '''
    ls = len(record[0][1]) - len(record[0][1]) % 3
    return sum([siteMeanNonSynPos(s, codonTable) 
               for s in [[x[1][i:i+3] for x in record] for i in range(0,ls,3)]])

def meanSynPos(record, codonTable=default_codon_table):
    '''Calculates the number of synonymous positions in an alignment.

    Calculates the mean number of synonymous sites in a set of aligned 
    sequences. The return value is the sum of the siteMeanSynPos for
    each site in the alignment.

    '''
    ls = len(record[0][1]) - len(record[0][1]) % 3
    return sum([siteMeanSynPos(s, codonTable) 
               for s in [[x[1][i:i+3] for x in record] for i in range(0,ls,3)]])

def siteMeanNonSynPos(site, codonTable=default_codon_table):
    '''Calculates the mean number of non-synonymous positions at a site.
    
    The return value is the average number of non-synonymous substitutions
    possible at that site.

    '''
    if any([(x in site) for x in codonTable.stop_codons]):
        return 0
    if any([("N" in x) for x in site]):
        return 0
    freqs = Counter(site)
    return sum([freqs[c]*nonSynPos(c,codonTable) for c in freqs.keys()])/float(len(site))

def siteMeanSynPos(site, codonTable=default_codon_table):
    '''Calculates the number of synonymous positions at a site.

    The return value is the average number of synonymous substitutions
    possible at that site.

    '''
    if any([(x in site) for x in codonTable.stop_codons]):
        return 0
    if any([("N" in x) for x in site]):
        return 0
    freqs = Counter(site)
    return sum([freqs[c]*synPos(c,codonTable) for c in freqs.keys()])/float(len(site))

def nonSynPos(c, codonTable=default_codon_table):
    '''Calculates the number of non-synonymous positions in a codon.
    
    The return value is the number of substitutions that will alter the amino 
    acid coded for by this codon.
    
    '''
    count = 0
    for i in [0,1,2]:
        for base in "ACTG":
            d = ''.join([c[j] if j != i else base for j in [0,1,2]])
            if codonTable.forward_table.get(d, c) != codonTable.forward_table.get(c, c):
                count += 1
    return count

def synPos(c, codonTable=default_codon_table):
    '''Calculates the number of synonymous positions in a codon.
    
    The return value is the number of substitutions that do not alter the 
    amino acid coded for by this codon.
    
    '''
    return 9 - nonSynPos(c,codonTable)

###
# 
# DIVERGENCE STATISTICS
# 

def div(record, codonTable=default_codon_table):
    '''Counts synonymous polymorphisms in a set of aligned sequences.'''
    ls = len(record[0][1]) - len(record[0][1]) % 3
    groups = [x[2] for x in record]
    if len(set(groups)) < 2:
        return 0    # Record has only one group or is empty
    sites = [[x[1][i:i+3] for x in record] for i in range(0,ls,3)]
    divS = 0    # Synonymous divergence
    divNS = 0   # Non-synonymous divergence
    for s in sites:
        if any([c in codonTable.stop_codons for c in s]) or any(['N' in c for c in s]):
            continue
        divS += countMutations(s, codonTable)[1]
        divNS += countMutations(s, codonTable)[0]
    return {'S':divS, 'NS': divNS}

def nonSynDiv(record, codonTable=default_codon_table):
    '''Counts synonymous polymorphisms in a set of aligned sequences.'''
    return div(record,codonTable)['NS']

def synDiv(record, codonTable=default_codon_table):
    '''Counts synonymous polymorphisms in a set of aligned sequences.'''
    return div(record,codonTable)['S']

###
# 
# POLYMORPHISM STATISTICS
# 

def synPol(record, codonTable=default_codon_table):
    '''Counts synonymous polymorphisms in a set of aligned sequences.'''
    ls = len(record[0][1]) - len(record[0][1]) % 3
    sites = [[x[1][i:i+3] for x in record] for i in range(0,ls,3)]
    usableSites = [s for s in sites if not 
            (any(['N' in c for c in s]) or any([c in s for c in codonTable.stop_codons]))]
    return sum([siteSynPol(s, [x[2] for x in record], codonTable) for s in usableSites])

def nonSynPol(record, codonTable=default_codon_table):
    '''Counts non-synonymous polymorphisms in a set of aligned sequences.'''
    ls = len(record[0][1]) - len(record[0][1]) % 3
    sites = [[x[1][i:i+3] for x in record] for i in range(0,ls,3)]
    usableSites = [s for s in sites if not 
            (any(['N' in c for c in s]) or any([c in s for c in codonTable.stop_codons]))]
    return sum([siteNonSynPol(s, [x[2] for x in record], codonTable) for s in usableSites])

def isPolymorphic(site, groups):
    '''Returns False if a site is polymorphic, otherwise True.
    
    A site is polymorphic if any two sequences in the same group differ at
    that site.

    '''
    return True if len(Counter(zip(site,groups))) > len(Counter(groups)) else False

def siteSynPol(site, groups, codonTable):
    '''Counts synonymous polymorphism at a site.'''
    count = 0
    for h in Counter(groups).keys():
        s = [c for (c,g) in zip(site, groups) if g == h]
        count += countMutations(s, codonTable)[1]
    return count

def siteNonSynPol(site, groups, codonTable):
    '''Counts non-synonymous polymorphism at a site.'''
    count = 0
    for h in Counter(groups).keys():
        s = [c for (c,g) in zip(site, groups) if g == h]
        count += countMutations(s, codonTable)[0]
    return count

