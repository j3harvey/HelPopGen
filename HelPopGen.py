#! /usr/bin/env python2.7

import Bio.SeqIO
import csv

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

def assignGroup( sequence ):
  '''
  Assigns a group to the given sequence. This is either:
    0:                  ignored
    999:                outgroup
    1:                  ingroup
  '''
  if seq.id.upper().startswith("ERATO"):
    return 999
  elif seq.id.upper().startswith("AMARYLLIS"):
    return 1
  elif seq.id.upper().startswith("AGLAOPE"):
    return 1
  else:
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

def qualityControlledData():
  '''
  Takes sequence records as input and yields only high quality data.
  '''
  for record in dataRecords:
    if skipQualityControl:
      yield record
    else:
      yield qced( record )

with open( inputFile, 'r' ) as fastaFile:
  for record in dataRecords( fastaFile ):
    qualityControl( record )
    myAlign = Align( string=str(record) )
    writeResults( popGenStats( myAlign ) )

