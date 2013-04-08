#! /usr/bin/env python2.7

import Bio.SeqIO
import time

inputFile = "ag5am5par1era4_allGenes_StopCodonsRemoved.fasta"

def fastaRecords():
  """
  This generator yields lists of (name,sequence) tuples from the input file
  e.g. [(name1,seq1),(name2,seq2),...]
  
  Each such list contains a set of sequences to be aligned for 
  caulculation of divergence and/or polymorphism statistics.
  
  This method should be modified by the user to correctly read their own data!
  """
  # We just load the entire fasta file and split records into groups of 32.
  # Your data may be arranged differently!
  # If so, you will have to implement your own function below or 
  # preprocess your data.
  seqs = Bio.SeqIO.parse( inputFile, 'fasta' )
  record = []
  for seq in seqs:
    if len(record) == 2:
      yield record
      record = [(seq.id, seq.seq.upper().tostring())]
    else:
      record.append((seq.id, seq.seq.upper().tostring()))
  if len(record) == 2:
    yield record

codes = { "AA": "A",
          "AC": "M",
          "AG": "R",
          "AT": "W",
          "CA": "M",
          "CC": "C",
          "CG": "S",
          "CT": "Y",
          "GA": "R",
          "GC": "S",
          "GG": "G",
          "GT": "K",
          "TA": "W",
          "TC": "Y",
          "TG": "K",
          "TT": "T",
          "NN": "N"}

def main():
  count = 0
  for record in fastaRecords():
    #if count > 100:
    #  break
    print count, record[0][0]
    try:
      ''.join([ codes[(x+y)] for (x,y) in zip(record[0][1], record[1][1]) ])
    except KeyError:
      print record
      (r1, r2) = record
      for i in range(len(r1[1])):
        if not codes.has_key( r1[1][i] + r2[1][i] ):
          print i
          yield r1,r2
          break
      break
    count += 1

if __name__ == "__main__":
  main()


