#! /usr/bin/env python2.7

import Bio.SeqIO
import time

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
    if len(record) == 2:
      yield record
      record = [(seq.id, seq.seq.upper().tostring())]
    else:
      record.append((seq.id, seq.seq.upper().tostring()))
  if len(record) == 2:
    yield record

def anyMultiChange( seq1, seq2 ):
  diffs = [ (x!=y) & (x!="N") & (y!="N")  for (x,y) in zip( list(seq1), list(seq2) ) ]
  return sum( [ int(sum(diffs[i:i+3]) > 1) for i in xrange(0,len(seq1),3) ] )

def numHet( seq1, seq2 ):
  diffs = [ (x!=y) & (x!="N") & (y!="N")  for (x,y) in zip( list(seq1), list(seq2) ) ]
  return sum( [ int(sum(diffs[i:i+3]) > 0) for i in xrange(0,len(seq1),3) ] )

t0 = time.time()
with open( inputFile, 'r' ) as fastaFile:
  nMultiChange = sum( [anyMultiChange( record[0][1], record[1][1] ) for record in fastaRecords( fastaFile )] )

t1 = time.time()
with open( inputFile, 'r' ) as fastaFile:
  nSites = sum( [len(record[0][1]) for record in fastaRecords( fastaFile )] )

t2 = time.time()
with open( inputFile, 'r' ) as fastaFile:
  nHet = sum( [numHet( record[0][1], record[1][1] ) for record in fastaRecords( fastaFile )] )

t3 = time.time()
