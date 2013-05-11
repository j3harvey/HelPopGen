#! /usr/bin/env python2.7

'''
Parse arguments before doing anything else.

This avoids lengthy imports when all you want is
the 'usage' information './HelPolGen.py --help'

'''
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("threshold", 
                    help="Minor allele frequency cutoff", 
                    nargs='?', 
                    type=float, 
                    default=0.0)
args = parser.parse_args()
THRESHOLD = args.threshold
INPUT_FILE = "data/ag5am5era4_allGenes_ambig.fasta"

import Bio.SeqIO
from Bio.Data import CodonTable
import csv
import sys
from itertools import product, chain
from collections import Counter
import DivPolStat

def dataRecords():
    '''A generator of aligned sequences for each locus in the data set.
    
    This generator yields lists of (name,sequence) tuples
    e.g. [(name1,seq1,group1),(name2,seq2,group2),...]
    
    Each such list contains a set of aligned sequences for 
    caulculation of divergence and/or polymorphism statistics.
    
    This method should be modified by the user to correctly read their own data!
    
    '''
    # We just load the entire fasta file and split records into groups of 16.
    # Your data may be arranged differently!
    # If so, you will have to implement your own function below or rearrange
    # your data.
    seqs = Bio.SeqIO.parse(INPUT_FILE, 'fasta')
    record = []
    for seq in seqs:
        if len(record) == 16:
            yield [x for x in record if x[2] != 0]
            record = [(seq.id, seq.seq.upper().tostring(), assignGroup(seq))]
        else:
            record.append((seq.id, seq.seq.upper().tostring(), assignGroup(seq)))
    if len(record) == 16:
        yield [x for x in record if x[2] != 0]

def assignGroup( seq ):
    '''Assigns a group number to a sequence object.
    
    Assigns a group number to seq, an object of class Bio.SeqRecord.SeqRecord.
    The group number is one of:
      0:                  ignored
      999:                outgroup
      1,2,3...:           ingroups
    
    Sequences with group number 0 are not used for any population genetic
    statistics. Sequences with group number 999 are excluded from the
    calculation of certain statistics. Other groups are used in all 
    calculations.
    
    You will have to write a custom method to suit your own data.
    
    '''
    # Ignore pardalinus (an alternative outgroup)
    # Also ignore the three inbred individuals.
    ignored = ["pardalinus", "Hmel", "aglaope.1_", "amaryllis.1_"]
    # All other amaryllis and aglaope are "in"
    ingroups = ["aglaope", "amaryllis"]
    # Erato is the outgroup
    outgroups = ["erato"]
    
    if any([seq.id.startswith(x) for x in ignored]):
        return 0
    elif any([seq.id.startswith(x) for x in ingroups]):
        return 1
    elif any([seq.id.startswith(x) for x in outgroups]):
        return 999
    else:
        #print("Unable to resolve group: ", seq.id)
        return 0

def phaseRecord(record):
    '''Takes a record as input and yields a record of phased data.'''
    # Repeat each id and group twice
    newIds    = chain(*[[x[0],x[0]] for x in record])
    newGroups = chain(*[[x[2],x[2]] for x in record])
    # Perform phasing at each site in the alignment and form a new alignment of the phased sequences
    ls        = len(record[0][1]) - len(record[0][1]) % 3
    sites     = [[x[1][i:i+3] for x in record] for i in range(0,ls,3)]
    newSeqs   = [''.join(newSeq) for newSeq in zip( *[phased(site) for site in sites])]
    # Return a new record of phased data, i.e. a list of tuples:
    # [(id1,seq1a,group1), (id1,seq1b,group1), (id1,seq2a,group2),...]
    return zip( newIds, newSeqs, newGroups )

_ambiguousBases = {"K": ["G","T"],
                   "M": ["A","C"],
                   "R": ["A","G"],
                   "S": ["C","G"],
                   "W": ["A","T"],
                   "Y": ["C","T"],}

def phased( site ):
    '''Takes a site (list of N codons) and outputs a phased version of that site (2N codons)'''
    newSite = list()
    for c in site:
        # For each codon in the site, count how many ambiguous bases are present
        numAmbig =  sum([b in _ambiguousBases.keys() for b in c])
        if numAmbig == 0:
            newSite.extend([c]*2)
        elif numAmbig == 1:
            # Simple phasing  e.g. "ATK" --> "ATG" and "ATT"
            newSite.extend(map(lambda z: ''.join(z),
                product(*[_ambiguousBases.get(x,x) for x in c])))
        elif numAmbig == 2:
            # Harder phasing  e.g. "KKT" --> ("GGT" and "TTT") or ("GTT" and "TGT")
            # This is yet to be implemented
            newSite.extend(['NNN','NNN'])
        else:
            # No phasing  e.g. "KYK" --> 4 possible genotypes
            # Give up, unless one possible phasing of the data matches other sequences 
            # at that site. This might be implemented later, but for now we throw out 
            # the small number of sites where this applies.
            newSite.extend(['NNN','NNN'])
    return newSite

def missingData(record):
    '''Returns a subset of a record with maximum usable data.
    
    missingData chooses a theshold t and removes all sequences with a proportion t or greater
    of data missing. Site with missing data are those that contain N's or stop codons.
    
    A range of values of t are tested, and a best value t* is chosen based on the total number
    of usable codons i.e. (number of sites with no missing data)*(number of sequences).
    
    The return value is a record with all sequences with proportion >t* missing data removed.
    
    '''
    thresholds = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
    counts     = [usableSites(record, t) for t in thresholds]
    thres      = max([t for (t, c) in zip(thresholds, counts) if c == max(counts)])
    return [x for x in record if x[2] != 0 and x[1].count("N")/float(len(x[1])) < thres]

def usableSites(record, thres=1.0):
    '''Calculates the amount of usable data in a record
    
    The amount of usable data in a record is defined as the number of sequences multiplied
    by the number of sites (aligned codons) with no N's and no stop codons, calculated
    after all sequences with a proportion >thres missing data are reomved.
    
    '''
    record = [x for x in record if x[2] != 0 and x[1].count("N")/float(len(x[1])) < thres]
    ns = len(record)
    count = 0
    if ns > 0:
        for i in range(0, len(record[0][1]), 3):
            if all([(x[1][i:i+3] not in ["TGA","TAG","TAA"] and 
                "N" not in x[1][i:i+3]) for x in record]):
                count += ns
    return count

def qualityCheck(record):
    '''Checks that a record is suitable for calculation of polymorphism and divergence statistics.
    
    Takes a record (a list of aligned sequences, with their identifiers and group
    numbers) as input. Returns 'True' if the data is suitable for calculation of
    polymorphism and divergence statistics, and 'False' otherwise.
    
    At present, this method only checks that there is more than one sequence in
    each group, so that we can detect polymorphism and divergence.
    
    '''
    groups = Counter([x[2] for x in record])
    # If there are too few 'in' or 'out' sequences to detect polymorphism, return False
    if (groups[1] < 2) or (groups[999] < 2):
        return False
    else:
        return True

def removeMinorAlleles(record, threshold):
    '''Replaces rare alleles with the major allele at each site in a record.
    
    For each site in the record, alleles with a frequency less than or equal to 
    threshold are repalced with the major allele at that site.
    
    '''
    ls = len(record[0][1]) - len(record[0][1]) % 3
    inSeqs   = [x for x in record if x[2] == 1]
    outSeqs  = [x for x in record if x[2] == 999]
    inSites  = [[x[1][i:i+3] for x in inSeqs] for i in range(0,ls,3)]
    outSites = [[x[1][i:i+3] for x in outSeqs] for i in range(0,ls,3)]
    newIn    = [''.join(newSeq) for newSeq in 
                  zip(*[removeMinorCodons(site,threshold) for site in inSites])]
    newOut   = [''.join(newSeq) for newSeq in 
                  zip(*[removeMinorCodons(site,threshold) for site in outSites])]
    ids      = [x[0] for x in inSeqs] + [x[0] for x in outSeqs]
    groups   = [x[2] for x in inSeqs] + [x[2] for x in outSeqs]
    return zip(ids, newIn + newOut, groups)

def removeMinorCodons(site,threshold):
    '''Replaces rare alleles with the major allele at a given site'''
    counts = Counter(site)
    codons = [c for c in site if c not in ['TAA','TAG','TGA'] and 'N' not in c]
    if len(codons) == 0:
        # Site has only stops and/or codons with missing data.
        return site
    majorAllele = max(set(codons), key=site.count)
    newSite = list()
    for c in site:
        if c in codons and site.count(c)/float(len(codons)) < threshold:
            newSite.append(majorAllele)
        else:
            newSite.append(c)
    return newSite

def polDivStats( record ):
    '''Calculates a dictionary of population statistics for a set of aligned sequences'''
    return { "name"                   : record[0][0],
             "ingroup_sequences"      : sum([x[2]==1 for x in record]),
             "outgroup_sequences"     : sum([x[2]==999 for x in record]),
             "sequence_length"        : len(record[0][1]),
             "usable_sequence_length" : 3*usableSites(record)/float(len(record)),
             "Ssites"                 : DivPolStat.meanNumSynPos(record),
             "NSsites"                : DivPolStat.meanNumNonSynPos(record),
             "P_N"                    : DivPolStat.numNonSynPol(record),
             "P_S"                    : DivPolStat.numSynPol(record),
             "D_N"                    : DivPolStat.numNonSynDiv(record),
             "D_S"                    : DivPolStat.numSynDiv(record),
           }

def trtv():
    '''Genome-wide estimate of Transition:Transversion (tr:tv) ratio'''
    tr = 0
    tv = 0
    for r in dataRecords():
        for i in range(len(r[0][1])):
              bases = Counter([x[1][i] for x in r]).keys()
              if 'A' in bases and 'G' in bases:
                  tr += 1
              if 'C' in bases and 'T' in bases:
                  tr += 1
              if ('A' in bases or 'G' in bases) and ('C' in bases or 'T' in bases):
                  tv += 1
    # Result: tr = 575644,  tv = 266020
    return (tr, tv)

###
# 
# MAIN ROUTINE
# 

def main():
    '''Calculate statistics for each set of aligned sequences in the input file.
    
    The main routine runs as follows:
    
    1) Load unphased data in an object 'records'
    2) Remove sequences with group number 0 (i.e. all 'ignored' sequences)
    3) Attempt to phase the data using the 'phased' function. Insert 'N's 
       whenever phasing fails.
    4) Attempt to maximise the amount of usable data ( = number of sequences 
       x number of usable sites) for each record. The strategy used is:
          i)   Pick a threshold t between 0 and 1.
          ii)  Discard sequences with a proportion of missing data greater 
               than t.
          iii) The amount of usable data = number of remaining sequences x 
               number of sites with no 'N's
          iv)  Repeat for other values of t and choose the best value.
       This search for an optimum threshold for each site is carried out by 
       the 'missingData' function. The sequences remaining after step (ii) 
       are stored in an object called 'usableData'.
    5) For each alignment in 'usableData' that passes certain quality control 
       tests ( the function 'qualityCheck' ) we calculate some polymorphism 
       and divergence statistics:
          i)   Number of synonymous and non-synonymous sites
          ii)  Number of polymorphic sites
          iii) Number of divergent sites
    
    '''
    
    records        = (record               for record in dataRecords())
    phased_records = (phaseRecord(record)  for record in records)
    #no_stops       = (stopsRemoved(record) for record in phased_records)
    major_alleles  = (removeMinorAlleles(record, THRESHOLD) for record in phased_records)
    usable_data    = (missingData(record)  for record in major_alleles)
    results        = (polDivStats(record)  for record in usable_data if qualityCheck(record))
    
    fieldNames = ["name",
                  "ingroup_sequences",
                  "outgroup_sequences",
                  "sequence_length",
                  "usable_sequence_length",
                  "Ssites",
                  "NSsites",
                  "P_N",
                  "P_S",
                  "D_N",
                  "D_S",]
    resultsWriter = csv.DictWriter(sys.stdout, fieldNames)
    resultsWriter.writerow({k:k for k in fieldNames})
    for line in results:
        resultsWriter.writerow(line)

if __name__ == "__main__":
    main()

