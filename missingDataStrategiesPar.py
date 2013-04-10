#! /usr/bin/env python2.7


'''

This script is intended for comparing strategies for dealing with 
missing data.

JJWelch writes:

If some individuals are missing for a site, then I would either 

  (i)   exclude the site
  (ii)  subsample the number of individuals for that whole allele 
        (such that, say, 6/8 alleles were available for a given gene), or 
  (iii) use an "average number of alleles" for that gene (e.g., 7.8 alleles).

(iii) is clearly a bodge, but it can be used.
I wouldn't use the frequency spectrum for this purpose as it is certain 
to vary between sites and loci, and I think (iii) is a preferable hack.

'''


import Bio.SeqIO
from Bio.Data import CodonTable
import csv
import time
from multiprocessing import Pool

inputFile = "ag5am5era4_allGenes_ambig.fasta"

def dataRecords():
  """
  This generator yields lists of (name,sequence,group) tuples
  e.g. [(name1,seq1,0),(name2,seq2,0),...] from the inputFile.
  
  Each such list contains a set of sequences to be aligned for 
  caulculation of divergence and/or polymorphism statistics.
  
  WARNING:
  You should check that this method will correctly read your data.
  If your data are arranged differently to our, you will have to modify
  this function.
  For our data, each record is obtained by reading in the next 16 
  fasta records from the file.
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
    1:                  ingroup
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

def usableSitesStar( args ):
  return usableSites( *args )

'''
def usableSites( record, excludeSites=False, thres=1.0 ):
  record = [x for x in record if x[2] != 0 and x[1].count("N")/float(len(x[1]))<thres]
  excluded = list()
  if excludeSites == True:
    for x in record:
      for i in range( 0, len(x[1]), 3):
        if x[1][i:i+3] in ["TGA", "TAA", "TAG"] or "N" in x[1][i:i+3]:
          excluded.append(i)
  count = 0
  for x in record:
    for i in range( 0, len(x[1]), 3 ):
      if i not in excluded:
        if x[1][i:i+3] not in ["TGA", "TAA", "TAG"] and "N" not in x[1][i:i+3]:
          count += 1
  return count
'''

def usableSites(record, excludeSites=False, thres=1.0):
  record = [x for x in record if x[2] != 0 and x[1].count("N")/float(len(x[1])) < thres]
  ns = len(record)
  count = 0
  if ns > 0:
    for i in range(0,len(record[0][1]),3):
      if excludeSites:
        if all( [(x[1][i:i+3] not in ["TGA","TAG","TAA"] and "N" not in x[1][i:i+3]) for x in record]):
          count += ns
      else:
        for x in record:
          if x[1][i:i+3] not in ["TGA", "TAG", "TAA"] and "N" not in x[1][i:i+3]:
            count += 1
  return count

def main(number=float('inf')):
  pool = Pool()
  
  t0 = time.time()
  
  def gen( excludeSites=False, thres=1.0 ):
    #for x in dataRecords():
    #  yield (x, excludeSites, thres )
    return [(x,excludeSites,thres) for x in dataRecords()]
  
  def countSites( thres ):
    x = pool.map( usableSitesStar, gen(True,thres))
    print time.time() - t0
    return x

  try:
    excludeSeqs = [countSites( thres ) for thres in [1.0, 0.9, 0.8, 0.7, 0.6, 0.4, 0.3, 0.2, 0.1] ]
    meanSeqNums = pool.map( usableSites, dataRecords() )
    print time.time() - t0
  except:
    pool.close()
    pool.join()
    raise
  
  return excludeSeqs + [meanSeqNums]

if __name__ == "__main__":
  results = main()
  print [sum(y) for y in results]
  with open("results2.csv", 'abw') as f:
    fwriter = csv.writer(f)
    for i in range(len(results[0])):
      fwriter.writerow( [results[j][i] for j in range(len(results))] )

