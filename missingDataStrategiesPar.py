#! /usr/bin/env python2.7

import Bio.SeqIO
from Bio.Data import CodonTable
import csv
import time
from multiprocessing import Pool

inputFile = "ag5am5era4_allGenes_ambig.fasta"

def dataRecords():
  # We just load the entire fasta file and split records into groups of 16.
  # Your data may be arranged differently!
  # If so, you will have to implement your own function below or rearrange
  # your data.
  seqs = Bio.SeqIO.parse( inputFile, 'fasta' )
  record = []
  records = []
  for seq in seqs:
    if len(record) == 16:
      records.append(record)
      record = [(seq.id, seq.seq.upper().tostring())]
    else:
      record.append((seq.id, seq.seq.upper().tostring()))
  if len(record) == 16:
    records.append(record)
  return records

def usableSites( record ):
  count = 0
  for x in record:
    for i in range( 0, len(x[1]), 3 ):
      if x[1][i:i+3] not in ["TGA", "TAA", "TAG"] and "N" not in x[1][i:i+3]:
        count += 1
  return count

def main():
  t0 = time.time()
  meanSeqNums = list()
  for record in dataRecords():
    meanSeqNums.append( usableSites(record) )
  print time.time() - t0
  return meanSeqNums

def mainParallel():
  t0 = time.time()
  pool = Pool()
  results = pool.map(usableSites, dataRecords())
  pool.close()
  pool.join()
  print time.time() - t0
  return results

if __name__ == "__main__":
  x = main()
  print sum(x)
  with open('notP.txt','w') as f:
    f.write('\n'.join([str(z) for z in x]))
  y = mainParallel()
  print sum(y)
  with open('notnotP.txt','w') as f:
    f.write('\n'.join([str(z) for z in y]))
 
