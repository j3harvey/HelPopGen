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

fRADmap = "Hmel1-1_chromosome.agp"
fImmune = "ImmunityGeneNames.csv"
fAllLoci = "ag5am5era4_allGenes_ambig.fasta"

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
    return tuple([tuple(l.split('\t')) for l in lines if
            l.split('\t')[0] == "chr" + str(n) and
            l.split('\t')[5].startswith("HE")])

chromosomes = {n:getChr(n) for n in range(1,21) + ["Z",]}

for Chr in chromosomes.keys():
    mapped_loci = list()
    for v in chromosomes[Chr]:
        loci = [l for l in idlines if l[1] == v[5]]
        if v[8] == "+":
            loci = sorted(loci, key=lambda x: x[3])
        else: # v[8] == "-"
            loci = sorted(loci, key=lambda x: -x[3])
        mapped_loci.extend([x[0] for x in loci])
    chromosomes[Chr] = mapped_loci

# get list of immune loci and for each one use the RAD map to fetch
# appropriate control genes
with open(fImmune, 'r') as f:
    immune_loci = [line.split(',')[0].rstrip() for line in f.readlines()[1:]]

def getControlgenes(locus):
    '''I the locus is mapped, get 20 loci either side as a control'''
    for Chr in chromosomes.values():
        if locus in Chr:
            index = [i for (i, l) in enumerate(Chr) if l == locus][0]
            lower = Chr[slice(index-20, index)]
            upper = Chr[slice(index+1, index+21)]
            return list(lower) + list(upper)
    print "Not mapped"
    return list()

