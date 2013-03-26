#! /usr/bin/env python2.7

import Bio.SeqIO

#HOPEFULLY NOT
#import egglib

def fastaRecords( fastaFile ):
  """
  fastaFile: a file containing input sequence data
  
  This generator yields lists of (name,sequence) tuples
  e.g. [(name1,seq1),(name2,seq2),...]
  
  Each such list contains a set of sequences to be aligned for 
  caulculation of divergence and/or polymorphism statistics.
  
  This method should be modified by the user to correctly read their own data!
  """
  # We just load the entire fasta file and split records into groups of 32.
  # Your data may be arranged differently!
  # If so, you will have to implement your own function below or rearrange
  # your data.
  seqs = Bio.SeqIO.parse( fastaFile )
  record = []
  for seq in seqs:
    if len(record) == 32:
      yield record
      record = []
    else:
      record.append((seq.id, seq.sequence.upper.tostring()))
  if len(record) == 32:
    yield record


with open( inputFile ) as fastaFile:
  for record in fastaRecords( fastaFile ):
    qualityControl( record )
    myAlign = egglib.Align( string=str(record) )
    writeResults( popGenStats( myAlign ) )
  
