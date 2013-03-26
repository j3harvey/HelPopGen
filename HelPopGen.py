#! /usr/bin/env python2.7

import Bio.SeqIO
import csv

#HOPEFULLY NOT
#import egglib

inputFile = "ag5am5par1era4_allGenes_StopCodonsRemoved.fasta"

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
  seqs = Bio.SeqIO.parse( fastaFile, 'fasta' )
  record = []
  for seq in seqs:
    if len(record) == 32:
      yield record
      record = [(seq.id, seq.seq.upper().tostring())]
    else:
      record.append((seq.id, seq.seq.upper().tostring()))
  if len(record) == 32:
    yield record


with open( inputFile, 'r' ) as fastaFile:
  for record in fastaRecords( fastaFile ):
    qualityControl( record )
    myAlign = egglib.Align( string=str(record) )
    writeResults( popGenStats( myAlign ) )

