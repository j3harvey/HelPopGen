#! /usr/bin/env python2.7

'''
This script finds appropriate genes to use as a control group for
carrying out McDonald-Kreitman tests on immunity genes.

A RAD map gives a chromosome and approximate location for each
Heliconius scaffold (about 80% of scaffolds are mapped). This script
maps coding sequences to chromosomes based on the RAD map, and
selects the 20 loci either side of a given immune locus to use as
a control.

'''

from re import search, findall

fRADmap = "../data/Hmel1-1_chromosome.agp"
fImmune = "../data/ImmunityGeneNames.csv"
fAllLoci = "../data/ag5am5era4_allGenes_ambig.fasta"

# read locus id's for H melpomene and extract for each one
# the scaffold name and position within the scaffold

with open(fAllLoci, 'r') as f:
    idlines = [l for l in f.readlines() if l.startswith(">")][::16]

idlines = [(search(r"HMEL[0-9]+", x).group(),
            search(r"HE[0-9]+", x).group(),
            int(x.split(';')[-1].split('-')[0]),
            int(x.split(';')[-1].split('-')[1]),) if "HMEL" in x else
           (search(r"HmGr[0-9]+", x).group(),
            search(r"HE[0-9]+", x).group(),
            int(x.split(';')[-1].split('-')[0]),
            int(x.split(';')[-1].split('-')[1]),) for x in idlines]

# load the RAD map data
with open(fRADmap, 'r') as f:
    lines = [l for l in f.readlines() if not l.startswith("#")]

mapped_scaffolds = [r.split('\t')[5] for r in lines if 
        "unmapped" not in r.split('\t')[0] and r.split('\t')[5].startswith("HE")]

def getChr(n):
    return [(l.split('\t')[5], int(l.split('\t')[1]), int(l.split('\t')[2]), l.split('\t')[8])
            for l in lines 
            if l.split('\t')[0] == "chr" + str(n)
            and l.split('\t')[5].startswith("HE")]

chromosomes_scaffolds = {n:tuple(getChr(n)) for n in range(1,21) + ["Z",]}

chromosomes = dict()
for Chr in chromosomes_scaffolds.keys():
    mapped_loci = list()
    for v in chromosomes_scaffolds[Chr]:
        loci = [l for l in idlines if l[1] == v[0]]
        if v[3] == "+":
            loci = sorted(loci, key=lambda x: x[3])
            mapped_loci.extend([(x[0], v[1] + x[2], v[1] + x[3]) for x in loci])
        else: # v[3] == "-"
            loci = sorted(loci, key=lambda x: -x[3])
            mapped_loci.extend([(x[0], v[2] - x[3], v[2] - x[2]) for x in loci])
    chromosomes[Chr] = mapped_loci

# get list of immune loci and for each one use the RAD map to fetch
# appropriate control genes
with open(fImmune, 'r') as f:
    immune_loci = [line.split(',')[0].rstrip() for line in f.readlines()[1:]]

def getPosition(locus):
    for k in chromosomes.keys():
        for l in chromosomes[k]:
            if locus == l[0]:
                return (k, l[1], l[2])

def mapped(loci):
    if isinstance(loci, str): loci = [loci,]
    mapped_loci = list()
    for locus in loci:
        if any([locus in [x[0] for x in ch] for ch in chromosomes.values()]):
            mapped_loci.append(locus)
    return mapped_loci

def getControlGenes(loci, n):
    '''Get n loci either side of each locus'''
    if isinstance(loci, str): loci = [loci,]
    control_genes = list()
    for locus in loci:
        for Chr in chromosomes.values():
            if locus in [x[0] for x in Chr]:
                index = [i for (i, l) in enumerate([x[0] for x in Chr]) if l == locus][0]
                control_genes.extend([x[0] for x in Chr][slice(index - n, index + n)])
    return list(set(control_genes) - set(loci))

