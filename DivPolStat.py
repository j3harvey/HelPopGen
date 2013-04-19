#! /usr/bin/env python2.7

import Bio.SeqIO
from Bio.Data import CodonTable
import csv
from collections import Counter
import itertools

###
# 
# POLYMORPHISM STATISTICS COPIED FROM BPP
# 

def numberChangesAtSiteBPP( site ):
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

def numberNonSynonymousChangesAtSiteBPP( site, codonTable=CodonTable.standard_dna_table ):
  '''
  A copy of the Bio++ method CodonSiteTools::numberOfNonSynonymousSubstitutions
  for counting non-synonymous changes, used by egglib's MKtable method.
  '''
  return len( set( [codonTable.forward_table.get(x) for x in site ] ) ) - 1

def numberSynonymousChangesAtSiteBPP( site, codonTable=CodonTable.standard_dna_table ):
  '''
  A copy of the Bio++ method for counting synonymous changes
  '''
  return numberChangesAtSiteBPP(site) - numberNonSynonymousChangesAtSiteBPP(site, codonTable)

def numberChangesBPP( record ):
  ls = len(record[0][1]) - len(record[0][1]) % 3
  return sum( [numberChangesAtSiteBPP(site) for site in [[x[1][i:i+3] for x in record] for i in range(0,ls,3)]] )

def numberSynonymousChangesBPP( record, codonTable=CodonTable.standard_dna_table ):
  ls = len(record[0][1]) - len(record[0][1]) % 3
  return sum( [numberSynonymousChangesAtSiteBPP(site) for site in [[x[1][i:i+3] for x in record] for i in range(0,ls,3)]] )

def numberNonSynonymousChangesBPP( record, codonTable=CodonTable.standard_dna_table ):
  ls = len(record[0][1]) - len(record[0][1]) % 3
  return sum( [numberNonSynonymousChangesAtSiteBPP(site) for site in [[x[1][i:i+3] for x in record] for i in range(0,ls,3)]] )

###
# 
# MUTATIONAL OPPORTTUNITY STATISTICS
# 

def meanNumNonSynPos( record, codonTable=CodonTable.standard_dna_table ):
  ls = len(record[0][1]) - len(record[0][1]) % 3
  return sum( [siteMeanNumNonSynPos(s, codonTable) for s in [[x[1][i:i+3] for x in record] for i in range(0,ls,3)]] )

def meanNumSynPos( record, codonTable=CodonTable.standard_dna_table ):
  ls = len(record[0][1]) - len(record[0][1]) % 3
  return sum( [siteMeanNumSynPos(s, codonTable) for s in [[x[1][i:i+3] for x in record] for i in range(0,ls,3)]] )

def siteMeanNumNonSynPos( site, codonTable=CodonTable.standard_dna_table ):
  if any([(x in site) for x in codonTable.stop_codons]):
    return 0
  if any([("N" in x) for x in site]):
    return 0
  freqs = Counter(site)
  return sum( [freqs[c]*numNonSynPos(c,codonTable) for c in freqs.keys()] )/float(len(site))

def siteMeanNumSynPos( site, codonTable=CodonTable.standard_dna_table ):
  if any([(x in site) for x in codonTable.stop_codons]):
    return 0
  if any([("N" in x) for x in site]):
    return 0
  freqs = Counter(site)
  return sum( [freqs[c]*numSynPos(c,codonTable) for c in freqs.keys()] )/float(len(site))

def numNonSynPos( c, codonTable=CodonTable.standard_dna_table ):
  try:
    return sum([sum([codonTable.forward_table.get(''.join([c[j] if j!=i else base for j in [0,1,2]]),c) != codonTable.forward_table[c] for base in "ACTG"]) for i in [0,1,2]])
  except:
    print c
    raise

def numSynPos( c, codonTable=CodonTable.standard_dna_table ):
  try:
    return 9 - numNonSynPos(c,codonTable)
  except:
    print c
    raise

###
# 
# DIVERGENCE STATISTICS
# 

def numDiv( record ):
  ls = len(record[0][1]) - len(record[0][1]) % 3
  groups = set( [x[2] for x in record] )
  if len(groups) < 2:
    return 0
  return sum([len(set(zip([x[1][i:i+3] for x in record], [x[2] for x in record]))) == len(groups) for i in range(0,ls,3)])

def numNonSynDiv( record, codonTable=CodonTable.standard_dna_table ):
  ls = len(record[0][1]) - len(record[0][1]) % 3
  groups = set( x[2] for x in record )
  if len(groups) < 2:
    return 0
  codonGroupPairs = [ set( zip( [x[2] for x in record], [x[1][i:i+3] for x in record] ) ) for i in range(0,ls,3) ]
  return sum( [ (len(x) == len(groups) and len(set([y[1] for y in x])) != 1) for x in codonGroupPairs ] )

def numSynDiv( record, codonTable=CodonTable.standard_dna_table ):
  return numDiv(record) - numNonSynDiv(record, codonTable)

###
# 
# POLYMORPHISM STATISTICS
# 

def numSynPol( record, codonTable=CodonTable.standard_dna_table ):
  ls = len(record[0][1]) - len(record[0][1]) % 3
  sites = [[x[1][i:i+3] for x in record] for i in range(0,ls,3)]
  usableSites = [s for s in sites if not ( any( ['N' in c for c in s] ) or any( [c in s for c in codonTable.stop_codons] ) ) ]
  return sum( [siteNumSynPol(s, [x[2] for x in record], codonTable) for s in usableSites] )

def numNonSynPol( record, codonTable=CodonTable.standard_dna_table ):
  ls = len(record[0][1]) - len(record[0][1]) % 3
  sites = [[x[1][i:i+3] for x in record] for i in range(0,ls,3)]
  usableSites = [s for s in sites if not ( any( ['N' in c for c in s] ) or any( [c in s for c in codonTable.stop_codons] ) ) ]
  return sum( [siteNumNonSynPol(s, [x[2] for x in record], codonTable) for s in usableSites] )

def isPolymorphic( site, groups ):
  return True if len(Counter(zip(site,groups))) > len(Counter(groups)) else False

def siteNumSynPol(site, groups, codonTable):
  return sum( [countMutations( [c for (c,g) in zip(site,groups) if g == h], codonTable )[1] for h in Counter(groups).keys()] )

def siteNumNonSynPol(site, groups, codonTable):
  return sum( [countMutations( [c for (c,g) in zip(site,groups) if g == h], codonTable )[0] for h in Counter(groups).keys()] )

def countMutations(site, ct=CodonTable.standard_dna_table):
  '''
  Uses a heuristic method (Kruskal's algorithm to find a set of mutations 
  that can explain the polymorphism observed at the site.

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
    if sets[u] != sets[v]:
      for (key,val) in sets.items():
        if u in val or v in val:
          sets[key] = sets[u] | sets[v]
      T.append((u,v))
  return (sum( [D[e][0] for e in set([tuple(sorted(t)) for t in T])] ), sum( [D[e][1] for e in set([tuple(sorted(t)) for t in T])] ))

'''
D is a dictionary of distances between codons.

For codons c and c, D[(c,d)] is a tuple:
the first element is the number of non-synonymous mutations needed to change
c to d, and the second element if the number of synonymous mutations needed.
'''
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

