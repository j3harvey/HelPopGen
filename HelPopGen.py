#! /usr/bin/env python2.7

import Bio.SeqIO
from Bio.Data import CodonTable
import csv
from collections import Counter

inputFile = "ag5am5par1era4_allGenes_StopCodonsRemoved.fasta"

# Set to True if data is already phased, otherwise set to False
dataIsPhased = True

# Set to True if the data has already been quality controlled: this will
# cause the script to attempt to use all data. Otherwise set to False.
skipQualityControl = False

###
#
# Define generators for the data sets.
# 
# These are a memory-efficient and safe way to iterate over the 
# alignments in our data files. Several generators are defined here:
#
# 1) `dataRecords` yields successive lists of aligned sequences. Each 
#     element in each list is a (name, sequence, group) tuple, where
#     group is one of "in", "out" or "ignore", and is determined by the
#     method whichGroup
# 2) `phasedRecords` performs phasing of the data (if necessary),
#     assuming parsimony with respect to non-synonymous sites, and that
#     there have been no stop codons
# 3) `qualityRecords` performs quality control check on the data. It 
#     yields only high quality data. It may remove some elements from each
#     record or throw out the record entirely.
#
# Having processed and filtered the data, we then move codon wise
# through each alignment, calculating diversion and polymorphism statistics.
#

def dataRecords():
  """
  fastaFile: a file containing input sequence data
  
  This generator yields lists of (name,sequence) tuples
  e.g. [(name1,seq1),(name2,seq2),...]
  
  Each such list contains a set of sequences to be aligned for 
  caulculation of divergence and/or polymorphism statistics.
  
  This method should be modified by the user to correctly read their own data!
  """
  # We just load the entire fasta file and split records into groups of 16.
  # Your data may be arranged differently!
  # If so, you will have to implement your own function below or rearrange
  # your data.
  seqs = Bio.SeqIO.parse( inputFile, 'fasta' )
  record = []
  for seq in seqs:
    if len(record) == 16:
      yield record
      record = [(seq.id, seq.seq.upper().tostring(), assignGroup(seq))]
    else:
      record.append((seq.id, seq.seq.upper().tostring(), assignGroup(seq)))
  if len(record) == 16:
    yield record

def assignGroup( seq ):
  '''
  Assigns a group to the given sequence. This is either:
    0:                  ignored
    999:                outgroup
    1,2,3...:           ingroups
  '''
  # Ignore pardalinus (an alternative outgroup)
  # Also ignore the three inbred individuals.
  ignored = ["pardalinus", "Hmel", "aglaope.1_", "amaryllis.1_"]
  # All other amaryllis and aglaope are "in"
  ingroups = ["amaryllis", "aglaope"]
  # Erato is the outgroup
  outgroups = ["erato"]
  
  if any( [seq.id.startswith(x) for x in ignored] ):
    return 0
  elif any( [seq.id.startswith(x) for x in ingroups] ):
    return 1
  elif any( [seq.id.startswith(x) for x in outgroups] ):
    return 999
  else:
    print "Unable to resolve group: ", seq.id
    return 0

def phasedData():
  '''
  Takes a generator of locus records as input and yields records of 
  phased data.
  '''
  for record in dataRecords:
    if dataIsPhased:
      yield record
    else:
      yield phased( record )

'''
If some individuals are missing for a site, then I would either 

  (i)   exclude the site
  (ii)  subsample the number of individuals for that whole allele 
        (such that, say, 6/8 alleles were available for a given gene), or 
  (iii) use an "average number of alleles" for that gene (e.g., 7.8 alleles).

(iii) is clearly a bodge, but it can be used.
I wouldn't use the frequency spectrum for this purpose as it is certain 
to vary between sites and loci, and I think (iii) is a preferable hack.
'''

def qualityControlledData():
  '''
  Takes sequence records as input and yields only high quality data.
  '''
  for record in dataRecords:
    if skipQualityControl:
      yield record
    else:
      yield qced( record )

###
# 
# POLYMORPHISM STATISTICS COPIED FROM BPP
# 

def numberChangesBPP( site ):
  '''
  A copy of the Bio++ method CodonSiteTools::numberOfSubsitutions 
  from the bpp-seq library. This is used in egglib's MKtable method.
  
  I don't get it either.
  '''
  Scodon = len(set(site)) - 1
  Sbases = (len(set([x[0] for x in site])) +
            len(set([x[1] for x in site])) +
            len(set([x[2] for x in site])) - 3)
  return max( Scodon, Sbases )

def numberNonSynonymousChangesBPP( site, codonTable=CodonTable.standard_dna_table ):
  '''
  A copy of the Bio++ method CodonSiteTools::numberOfNonSynonymousSubstitutions
  for counting non-synonymous changes, used by egglib's MKtable method.
  '''
  return len( set( [codonTable.forward_table[x] for x in site 
                                    if x in codonTable.forward_table.keys()] ) ) - 1

def numberSynonymousChangesBPP( site, codonTable=CodonTable.standard_dna_table ):
  '''
  A copy of the Bio++ method for counting synonymous changes
  '''
  return numberChangesBPP(site) - numberNonSynonymousChangesBPP(site, codonTable)

def meanNumSynPos( site, codonTable=CodonTable.standard_dna_table ):
  if any([(x in site) for x in codonTable.stop_codons]):
    return 0
  if any([("N" in x) for x in site]):
    return 0
  freqs = Counter(site)
  return float(sum( [freqs[c]*numSynPos(c,codonTable) for c in freqs.keys()] )/len(site))

def numNonSynPos( c, codonTable=CodonTable.standard_dna_table ):
  return sum([sum([codonTable.forward_table[''.join([c[j] if j!=i else base for j in [0,1,2]])] != codonTable.forward_table[c] for base in "ACTG"]) for i in [0,1,2]])

def numSynPos( c, codonTable=CodonTable.standard_dna_table ):
  return 9 - numNonSynPos(c,codonTable)

def numDiv( record ):
  ls = len(record[0][1])
  groups = set( x[2] for x in record )
  if len(groups) < 2:
    return 0
  return sum([len(set(zip([x[1][i:i+3] for x in record], x[2] for x in record]))) == len(groups) for in range(0,ls,3)])

def numNonSynDiv( record, codonTable=CodonTable.standard_dna_table ):
  ls = len(record[0][1])
  groups = set( x[2] for x in record )
  if len(groups) < 2:
    return 0
  codonGroupPairs = [ set( zip( [x[2] for x in record], x[1][i:i+3] for x in record] ) ) for i in range(0,ls,3) ]
  return sum( [ (len(x) == len(groups) and len(set([y[1] for y in x])) != 1) for x in codonGroupPairs ] )

def numSynDiv( record, codonTable=CodonTable.standard_dna_table ):
  return numDiv(record) - numNonSynDiv(record, codonTable)

###
# 
# MAIN ROUTINE
# 

def main():
  with open( inputFile, 'r' ) as fastaFile:
    records = [record for record in dataRecords()]

if __name__ == "__main__":
  main()

