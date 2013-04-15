#! /usr/bin/env python2.7

import Bio.SeqIO
from Bio.Data import CodonTable
import csv
from itertools import product, chain
from collections import Counter

INPUT_FILE = "ag5am5par1era4_allGenes_StopCodonsRemoved.fasta"
OUTPUT_FILE = INPUT_FILE.rstrip(".fasta") + ".csv"

# Set to True if data is already phased, otherwise set to False
dataIsPhased = True

trtv = {transitions:   [("A","G"),("C","T"),("G","A"),("T","C")],
        transversions: [("A","C"),("A","T"),("C","A"),("C","G"),("G","C"),("G","T"),("T","A"),("T","G")],
       }

class Align:
  '''
  The main data class.

  Acts like a list of aligned sequences, or a 2-dimensional array of codons.
  '''
  def __init__(self,ids,seqs,groups,description=None):
    self.ids = ids
    self.seqs = seqs
    self.groups = groups
    self.desc = description
    self.ns = len(ids)
    self.ls = len(seqs[0])
  def __iter__(self):
    return iter(zip(self.ids, self.seqs, self.groups))
  def sites(self):
    return [[x[i:i+3] for x in self.seqs] for i in range(0,self.ls,3)]
  def __getsites__(self,*args):
    return Align(self.ids, [''.join(x) for x in zip(*self.sites()[args[0]])], self.groups, self.desc)
  def __getitem__(self, *args):
    if len(args[0]) > 2:
      raise IndexError("Align indices must be of length 1 or 2")
    if type(args[0]) in (int, slice):
      return self.seqs[args[0]]
    elif args[0][0] == slice(None,None,None):
      return self.__getsites__(args[0][1])
    else:
      s, t = args[0]
      return Align(self.ids[s], self.seqs[s], self.groups[s],self.desc).__getsites__(t)

def dataRecords():
  '''
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
  seqs = Bio.SeqIO.parse( INPUT_FILE, 'fasta' )
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
  Assigns a group number to seq, an object of class Bio.SeqRecord.SeqRecord.
  The group number is one of:
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
    #print("Unable to resolve group: ", seq.id)
    return 0

def phaseRecord(record):
  '''
  Takes a record as input and yields a record of phased data.
  '''
  # Repeat each id and group twice
  newIds    = chain( *[[x[0],x[0]] for x in record] )
  newGroups = chain( *[[x[2],x[2]] for x in record] )
  # Perform phasing at each site in the alignment and form a new alignment of the phased sequences
  ls        = len(record[0][1])
  sites     = [[x[1][i:i+3] for x in record] for i in range(0,ls,3)]
  newSeqs   = [''.join(newSeq) for newSeq in zip( *[phased(site) for site in sites] ) ]
  # Return a new record of phased data, i.e. a list of tuples:
  # [(id1,seq1a,group1), (id1,seq1b,group1), (id1,seq2a,group2),...]
  return zip( newIds, newSeqs, newGroups )

ambiguousBases = {"K": ["G","T"],
                  "M": ["A","C"],
                  "R": ["A","G"],
                  "S": ["C","G"],
                  "W": ["A","T"],
                  "Y": ["C","T"],}

def phased( site ):
  newSite = list()
  for c in site:
    # For each codon in the site, count how many ambiguous bases are present
    numAmbig =  sum( [b in ambiguousBases.keys() for b in c] )
    if numAmbig == 0:
      newSite.extend([c]*2)
    if numAmbig == 1:
      # Simple phasing  e.g. "ATK" --> "ATG" and "ATT"
      newSite.extend( map( lambda z: ''.join(z),
                              product(*[ambiguousBases.get(x,x) for x in c]) ) )
    if numAmbig == 2:
      # Harder phasing  e.g. "KKT" --> ("GGT" and "TTT") or ("GTT" and "TGT")
      # This is yet to be implemented
      continue
    else:
      # No phasing  e.g. "KYK" --> 4 possible genotypes
      # Give up, unless one possible phasing of the data matches other sequences 
      # at that site. This might be implemented later, but for now we throw out 
      # the small number of sites where this applies.
      continue
  return newSite

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

def missingData(record):
  thresholds = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
  counts     = [usableSites(record, t) for t in thresholds]
  thres      = max( [t for (t,c) in zip(thresholds,counts) if c == max(counts)] )
  return [x for x in record if x[2] != 0 and x[1].count("N")/float(len(x[1])) < thres]

def usableSites(record, thres=1.0):
  record = [x for x in record if x[2] != 0 and x[1].count("N")/float(len(x[1])) < thres]
  ns = len(record)
  count = 0
  if ns > 0:
    for i in range(0,len(record[0][1]),3):
      if all( [(x[1][i:i+3] not in ["TGA","TAG","TAA"] and "N" not in x[1][i:i+3]) for x in record]):
        count += ns
  return count

def qualityCheck(record):
  '''
  Takes a record (a list of aligned sequences, with their identifiers and group
  numbers) as input. Returns 'True' if the data is suitable for calculation of
  polymorphism and divergence statistics, and 'False' otherwise.

  At present, this method only checks that there is more than one sequence in
  each group, so that we can detect polymorphism and divergence.
  '''
  groups = collections.Counter([x[2] for x in record])
  # If there are too few 'in' or 'out' sequences to detect polymorphism, return False
  if (len(groups[1]) < 2) or (len(groups[999]) < 2):
    return False
  else:
    return True

def polDivStats( record ):
  return { "name"                   : record[0][0],
           "ingroup_sequences"      : len( [x[2]==1 for x in record] ),
           "outgroup_sequences"     : len( [x[2]==999 for x in record] ),
           "sequence_length"        : len( record[0][1] ),
           "usable_sequence_length" : 3*usableSites(record),
           "Ssites"                 : DivPolStat.meanNumSynPos(record),
           "NSsites"                : DivPolStat.meanNumNonSynPos(record),
           "P_N"                    : DivPolStat.numNonSynPol(record),
           "P_S"                    : DivPolStat.numSynPol(record),
           "D_N"                    : DivPolStat.numNonSynDiv(record),
           "D_S"                    : DivPolStat.numSynDiv(record),
         }

###
# 
# MAIN ROUTINE
# 

def main():
  '''
  The main routine runs as follows:

  1) Load unphased data in an object 'records'
  2) Remove sequences with group number 0 (i.e. all 'ignored' sequences)
  3) Attempt to phase the data using the 'phased' function. Insert 'N's whenever phasing fails.
  4) Attempt to maximise the amount of usable data ( = number of sequences x number of usable
     sites) for each record. The strategy used is:
        i)   Pick a threshold t between 0 and 1.
        ii)  Discard sequences with a proportion of missing data greater than t.
        iii) The amount of usable data = number of remaining sequences x number of sites with no 'N's
        iv)  Repeat for other values of t and choose the best value.
     This search for an optimum threshold for each site is carried out by the 'missingData' function.
     The sequences remaining after step (ii) are stored in an object called 'usableData'.
  5) For each alignment in 'usableData' that passes certain quality control tests ( the function 
     'qualityCheck' ) we calculate some polymorphism and divergence statistics:
        i)   Number of synonymous and non-synonymous sites
        ii)  Number of polymorphic sites
        iii) Number of divergent sites
  '''
  records        = (record               for record in dataRecords())
  phased_records = (phased(record)       for record in records)
  no_stops       = (stopsRemoved(record) for record in phased_records)
  usable_data    = (missingData(record)  for record in phased_records)
  results        = (polDivStats(record)  for record in usable data if qualityCheck(record))

  with open( OUTPUT_FILE, 'w' ) as f:
    resultsWriter = csv.DictWriter( open(OUTPUT_FILE,'w') )
    resultsWriter.writerows(results)

if __name__ == "__main__":
  main()

