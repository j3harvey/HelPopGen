#! /usr/bin/env python2.7

import Bio.SeqIO
import time

inputFile = "ag5am5par1era4_allGenes_StopCodonsRemoved.fasta"

def fastaRecords():
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

def printHardSites():
  count = 0
  for record in fastaRecords():
    for i in range(0,len(record[0][1]),3):
      site = [x[1][i:i+3] for x in record]
      if any( [ sum( [b in "YKMSRW" for b in c] ) > 1 for c in site ] ):
        count += 1
        print count, ":  ", site

def anyMultiChange( seq1, seq2 ):
  diffs = [ (x!=y) & (x!="N") & (y!="N")  for (x,y) in zip( list(seq1), list(seq2) ) ]
  return sum( [ int(sum(diffs[i:i+3]) > 1) for i in xrange(0,len(seq1),3) ] )

def anyThreeChange( seq1, seq2 ):
  diffs = [ (x!=y) & (x!="N") & (y!="N")  for (x,y) in zip( list(seq1), list(seq2) ) ]
  return sum( [ int(sum(diffs[i:i+3]) > 2) for i in xrange(0,len(seq1),3) ] )

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
with open( inputFile, 'r' ) as fastaFile:
  nThree = sum( [anyThreeChange( record[0][1], record[1][1] ) for record in fastaRecords( fastaFile )] )

t4 = time.time()


def threeChange( seq1, seq2 ):
  diffs = [ (x!=y) & (x!="N") & (y!="N")  for (x,y) in zip( list(seq1), list(seq2) ) ]
  for i in xrange(0,len(seq1),3):
    if int(sum(diffs[i:i+3]) > 2):
      return i


with open( inputFile, 'r' ) as fastaFile:
  for record in fastaRecords():
    if threeChange( record[0][1], record[1][1] ):
      i = threeChange( record[0][1], record[1][1] )
      print record[0][0], i, len(record[0][1]), record[0][1][i:i+3], record[1][1][i:i+3]


### Count number of usable sites


def noNs(c):
    if any([y in c for y in "NKMRSWY"]):
        return False
    elif c in ["TAA","TAG","TGA"]:
        return False
    else:
        return True

count = 0
for r in dataRecords():
    ls = len(r[0][1]) - len(r[0][1]) % 3
    grps = [x[2] for x in r]
    for i in range(0,ls,3):
        s = phased([x[1][i:i+3] for x in r])
        if all([noNs(c) for c in s]):
             count += 1

#Result: 4269019

### Count number of hard to phase sites

def Sbases(site):
    return (len(set([x[0] for x in site])) + 
            len(set([x[1] for x in site])) +
            len(set([x[2] for x in site])) - 3)

def Scodon(site):
    return len(set(site)) - 1

def noNs(c):
    if any([y in c for y in "NKMRSWY"]):
        return False
    elif c in ["TAA","TAG","TGA"]:
        return False
    else:
        return True

count = 0
for r in dataRecords():
    count += 1
    ls = len(r[0][1]) - len(r[0][1]) % 3
    grps = [x[2] for x in r]
    for i in range(0,ls,3):
        s = phased([x[1][i:i+3] for x in r])
        if all([noNs(c) for c in s]):
            if len(set(zip(s,grps))) > len(set(grps)):
                if Sbases(s) != Scodon(s):
                    count += 1

#Result: 14294 out of 5785206 total sites

### Count number of polymorphic sites

def noNs(c):
    if any([y in c for y in "NKMRSWY"]):
        return False
    elif c in ["TAA","TAG","TGA"]:
        return False
    else:
        return True

count = 0
for r in dataRecords():
    ls = len(r[0][1]) - len(r[0][1]) % 3
    grps = [x[2] for x in r]
    for i in range(0,ls,3):
        s = phased([x[1][i:i+3] for x in r])
        if all([noNs(c) for c in s]):
          if len(set(s)) > 1:
            count += 1

#Result: 1114119 out of 5785206 total sites


### Count spectrum of polymorphic sites

def noNs(c):
    if any([y in c for y in "NKMRSWY"]):
        return False
    elif c in ["TAA","TAG","TGA"]:
        return False
    else:
        return True

count = dict()
for r in dataRecords():
    ls = len(r[0][1]) - len(r[0][1]) % 3
    grps = [x[2] for x in r]
    for i in range(0,ls,3):
        s = phased([x[1][i:i+3] for x in r])
        if all([noNs(c) for c in s]):
          if len(set(s)) > 3:
            print len(set(s)), list(set([(g[0],g[1]) for g in zip(itertools.chain(*[[x,x] for x in grps]), s)]))
          if len(set(s)) in count.keys():
            count[len(set(s))] += 1
          else:
            count[len(set(s))] = 1

#Result: {1: 3154900, 2: 1005742, 3: 101029, 4: 7051, 5: 287, 6: 9, 7: 1}


def mutPath(c1, c2, CT=CodonTable.standard_dna_table):
  '''
  Calculates a parsimonious sequence of of mutations from
  codon c1 to codon c2.
  '''
  FT = CT.forward_table
  nequal = [c1[i] != c2[i] for i in range(3)]
  if sum(nequal) == 0:
    return (list(),0)
  elif sum(nequal) == 1:
    return ([c1,c2],1)
  elif sum(nequal) == 2:
    d1 = ''.join([c2[0],c1[1],c1[2]]) if nequal[0] == 1 else ''.join([c1[0],c2[1],c1[2]])
    d2 = ''.join([c1[0],c1[1],c2[2]]) if nequal[2] == 1 else ''.join([c1[0],c2[1],c1[2]])
    numNonSyn1 = len(list(itertools.groupby( [FT.get(x,'*') for x in [c1,d1,c2]] )))
    numNonSyn2 = len(list(itertools.groupby( [FT.get(x,'*') for x in [c1,d2,c2]] )))
    if numNonSyn1 >= numNonSyn2 and d2 in FT.keys():
      return ([c1,d2,c2], numNonSyn2-1)
    elif d1 in FT.keys():
      return ([c1,d1,c2], numNonSyn1-1)
    else:
      print c1,c2, ':   ', "Codons not valid"
  else:
    interCodons = [(''.join([c1[0],c1[1],c2[2]]),''.join([c1[0],c2[1],c2[2]])),
                   (''.join([c1[0],c1[1],c2[2]]),''.join([c2[0],c1[1],c2[2]])),
                   (''.join([c1[0],c2[1],c1[2]]),''.join([c1[0],c2[1],c2[2]])),
                   (''.join([c1[0],c2[1],c1[2]]),''.join([c2[0],c2[1],c1[2]])),
                   (''.join([c2[0],c1[1],c1[2]]),''.join([c2[0],c1[1],c2[2]])),
                   (''.join([c2[0],c1[1],c1[2]]),''.join([c2[0],c2[1],c1[2]])),]
    numNonSyn = [len(list(itertools.groupby( [FT.get(x,'*') for x in [c1,]+list(interCodons[i])+[c2,]] ))) for i in range(6)]
    numNonSyn = sorted(zip(interCodons,numNonSyn),key=lambda x: x[1])
    for x in numNonSyn:
      if '*' not in x[0]:
        return ([c1,] + list(x[0]) + [c2,], x[1]-1)
    return None
        
Muts = dict()
for c1 in CodonTable.standard_dna_table.forward_table.keys():
  for c2 in CodonTable.standard_dna_table.forward_table.keys():
    Muts[ (c1,c2) ] = mutPath(c1,c2)

D = dict()
FT=CodonTable.standard_dna_table.forward_table
for c in FT.keys():
  for i in range(3):
    for b in "ACGT":
      d = ''.join([c[j] if j != i else b for j in range(3)])
      if d == c:
        D[(c,d)] = (0,0)
      elif d not in FT.keys():
        continue #D[(c,d)] = float('Inf')
      else:
        D[(c,d)] = (0,1) if FT[c] == FT[d] else (1,0)
for c1 in FT.keys():
  for c2 in FT.keys():
    nequal = sum([c1[i] != c2[i] for i in range(3)])
    if nequal < 2:
      continue
    elif nequal == 2:
      D[(c1,c2)] = min( [(D[(c1,d)][0] + D[(d,c2)][0], D[(c1,d)][1] + D[(d,c2)][1]) for d in FT.keys() if (c1,d) in D.keys() and (d,c2) in D.keys()] )
    else:
      continue
for c1 in FT.keys():
  for c2 in FT.keys():
    nequal = sum([c1[i] != c2[i] for i in range(3)])
    if nequal == 3:
      D[(c1,c2)] = min( [(D[(c1,d)][0] + D[(d,c2)][0], D[(c1,d)][1] + D[(d,c2)][1]) for d in FT.keys() if (c1,d) in D.keys() and (d,c2) in D.keys()] )

def mutPath(c1, c2, CT=CodonTable.standard_dna_table):
  '''
  Calculates a parsimonious sequence of of mutations from
  codon c1 to codon c2.
  '''
  FT = CT.forward_table
  nequal = [c1[i] != c2[i] for i in range(3)]
  if sum(nequal) == 0:
    return (list(),0)
  elif sum(nequal) == 1:
    return ([c1,c2],1)
  elif sum(nequal) == 2:
    d1 = ''.join([c2[0],c1[1],c1[2]]) if nequal[0] == 1 else ''.join([c1[0],c2[1],c1[2]])
    d2 = ''.join([c1[0],c1[1],c2[2]]) if nequal[2] == 1 else ''.join([c1[0],c2[1],c1[2]])
    numNonSyn1 = len(list(itertools.groupby( [FT.get(x,'*') for x in [c1,d1,c2]] )))
    numNonSyn2 = len(list(itertools.groupby( [FT.get(x,'*') for x in [c1,d2,c2]] )))
    if numNonSyn1 >= numNonSyn2 and d2 in FT.keys():
      return ([c1,d2,c2], numNonSyn2-1)
    elif d1 in FT.keys():
      return ([c1,d1,c2], numNonSyn1-1)
    else:
      print c1,c2, ':   ', "Codons not valid"
  else:
    interCodons = [(''.join([c1[0],c1[1],c2[2]]),''.join([c1[0],c2[1],c2[2]])),
                   (''.join([c1[0],c1[1],c2[2]]),''.join([c2[0],c1[1],c2[2]])),
                   (''.join([c1[0],c2[1],c1[2]]),''.join([c1[0],c2[1],c2[2]])),
                   (''.join([c1[0],c2[1],c1[2]]),''.join([c2[0],c2[1],c1[2]])),
                   (''.join([c2[0],c1[1],c1[2]]),''.join([c2[0],c1[1],c2[2]])),
                   (''.join([c2[0],c1[1],c1[2]]),''.join([c2[0],c2[1],c1[2]])),]
    numNonSyn = [len(list(itertools.groupby( [FT.get(x,'*') for x in [c1,]+list(interCodons[i])+[c2,]] ))) for i in range(6)]
    numNonSyn = sorted(zip(interCodons,numNonSyn),key=lambda x: x[1])
    for x in numNonSyn:
      if '*' not in x[0]:
        return ([c1,] + list(x[0]) + [c2,], x[1]-1)
    return None

def countMutations(site, ct=CodonTable.standard_dna_table):
  '''
  Uses a heuristic method to find a set of mutations that can explain
  the polymorphism at the site.

  The return value is a tuple: the first element is the number of
  non-synonymous mutations, and the second element is the number of
  sysnonymous mutations.

  This is an example of a Steiner Tree Problem, which is know to be
  NP-hard.
  '''
  E = set( [tuple(sorted((c,d))) for c in site for d in site if d != c] )
  E = sorted(E, key=lambda e: D[e])
  T = list()
  sets = { c : set([c,]) for c in site }
  if len(E) == 0:
    return (0,0)
  for (u,v) in E:
    print u,v,sets[u],sets[v]
    if sets[u] == sets[v]:
      print "Discarding edge:  ",u,v
    else:
      print "Adding edge:  ", u,v
      for (key,val) in sets.items():
        if u in val or v in val:
          sets[key] = sets[u] | sets[v]
      T.append((u,v))
  return (sum( [D[e][0] for e in set([tuple(sorted(t)) for t in T])] ), sum( [D[e][1] for e in set([tuple(sorted(t)) for t in T])] ))

Muts = dict()
for c1 in CodonTable.standard_dna_table.forward_table.keys():
  for c2 in CodonTable.standard_dna_table.forward_table.keys():
    Muts[ (c1,c2) ] = mutPath(c1,c2)

D = dict()
FT=CodonTable.standard_dna_table.forward_table
for c in FT.keys():
  for i in range(3):
    for b in "ACGT":
      d = ''.join([c[j] if j != i else b for j in range(3)])
      if d == c:
        D[(c,d)] = (0,0)
      elif d not in FT.keys():
        continue #D[(c,d)] = float('Inf')
      else:
        D[(c,d)] = (0,1) if FT[c] == FT[d] else (1,0)
for c1 in FT.keys():
  for c2 in FT.keys():
    nequal = sum([c1[i] != c2[i] for i in range(3)])
    if nequal < 2:
      continue
    elif nequal == 2:
      D[(c1,c2)] = min( [(D[(c1,d)][0] + D[(d,c2)][0], D[(c1,d)][1] + D[(d,c2)][1]) for d in FT.keys() if (c1,d) in D.keys() and (d,c2) in D.keys()] )
    else:
      continue
for c1 in FT.keys():
  for c2 in FT.keys():
    nequal = sum([c1[i] != c2[i] for i in range(3)])
    if nequal == 3:
      D[(c1,c2)] = min( [(D[(c1,d)][0] + D[(d,c2)][0], D[(c1,d)][1] + D[(d,c2)][1]) for d in FT.keys() if (c1,d) in D.keys() and (d,c2) in D.keys()] )

D = dict()
P = dict()
FT=CodonTable.standard_dna_table.forward_table
for c in FT.keys():
  for i in range(3):
    for b in "ACGT":
      d = ''.join([c[j] if j != i else b for j in range(3)])
      if d == c:
        D[(c,c)] = (0,0)
        P[(c,c)] = list()
      elif d not in FT.keys():
        continue #D[(c,d)] = float('Inf')
      else:
        D[(c,d)] = (0,1) if FT[c] == FT[d] else (1,0)
        P[(c,d)] = [c,]
for c1 in FT.keys():
  for c2 in FT.keys():
    nequal = sum([c1[i] != c2[i] for i in range(3)])
    if nequal < 2:
      continue
    elif nequal == 2:
      D[(c1,c2)] = min( [(D[(c1,d)][0] + D[(d,c2)][0], D[(c1,d)][1] + D[(d,c2)][1]) for d in FT.keys() if (c1,d) in D.keys() and (d,c2) in D.keys()] )
      for d in FT.keys():
        if (c1,d) in D.keys() and (d,c2) in D.keys() and d not in (c1,c2):
          if (D[(c1,d)][0] + D[(d,c2)][0], D[(c1,d)][1] + D[(d,c2)][1]) == D[(c1,c2)]:
            try:
              P[(c1,c2)] = P[(c1,d)] + P[(d,c2)]
            except:
              print c1,d,c2
              raise
            break
    else:
      continue
for c1 in FT.keys():
  for c2 in FT.keys():
    nequal = sum([c1[i] != c2[i] for i in range(3)])
    if nequal == 3:
      D[(c1,c2)] = min( [(D[(c1,d)][0] + D[(d,c2)][0], D[(c1,d)][1] + D[(d,c2)][1]) for d in FT.keys() if (c1,d) in D.keys() and (d,c2) in D.keys()] )
      for d in FT.keys():
        if (c1,d) in D.keys() and (d,c2) in D.keys() and d not in (c1,c2):
          if (D[(c1,d)][0] + D[(d,c2)][0], D[(c1,d)][1] + D[(d,c2)][1]) == D[(c1,c2)]:
            try:
              P[(c1,c2)] = P[(c1,d)] + P[(d,c2)]
            except:
              print c1,d,c2
            break
 

