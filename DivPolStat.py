#! /usr/bin/env python2.7

import Bio.SeqIO
from Bio.Data import CodonTable
import csv
from collections import Counter

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
  return len( set( [codonTable.forward_table[x] for x in site 
                                    if x in codonTable.forward_table.keys()] ) ) - 1

def numberSynonymousChangesAtSiteBPP( site, codonTable=CodonTable.standard_dna_table ):
  '''
  A copy of the Bio++ method for counting synonymous changes
  '''
  return numberChangesAtSiteBPP(site) - numberNonSynonymousChangesAtSiteBPP(site, codonTable)

def numberChangesBPP( record ):
  ls = len(record[0][1])
  return sum( [numberChangesAtSiteBPP(site) for site in [[x[1][i:i+3] for x in record] for i in range(0,ls,3)]] )

def numberSynonymousChangesBPP( record, codonTable=CodonTable.standard_dna_table ):
  ls = len(record[0][1])
  return sum( [numberSynonymousChangesAtSiteBPP(site) for site in [[x[1][i:i+3] for x in record] for i in range(0,ls,3)]] )

def numberNonSynonymousChangesBPP( record, codonTable=CodonTable.standard_dna_table ):
  ls = len(record[0][1])
  return sum( [numberNonSynonymousChangesAtSiteBPP(site) for site in [[x[1][i:i+3] for x in record] for i in range(0,ls,3)]] )

###
# 
# MUTATIONAL OPPORTTUNITY STATISTICS
# 

def meanNumNonSynPos( record, codonTable=CodonTable.standard_dna_table ):
  if any([(x in site) for x in codonTable.stop_codons]):
    return 0
  if any([("N" in x) for x in site]):
    return 0
  freqs = Counter(site)
  return float(sum( [freqs[c]*numNonSynPos(c,codonTable) for c in freqs.keys()] )/len(site))

def meanNumSynPos( record, codonTable=CodonTable.standard_dna_table ):
  if any([(x in site) for x in codonTable.stop_codons]):
    return 0
  if any([("N" in x) for x in site]):
    return 0
  freqs = Counter(site)
  return float(sum( [freqs[c]*numSynPos(c,codonTable) for c in freqs.keys()] )/len(site))

def numNonSynPos( c, codonTable=CodonTable.standard_dna_table ):
  return sum([sum([codonTable.forward_table[''.join([c[j] if j!=i else base for j in [0,1,2]])] != codonTable.forward_table[c] for base in "ACTG"]) for i in [0,1,2]])

def numSynPos( c, codonTable=CodonTable.standard_dna_table ):
  return 9 - numNonSynPos(c,codonTable)

###
# 
# DIVERGENCE STATISTICS
# 

def numDiv( record ):
  ls = len(record[0][1])
  groups = set( [x[2] for x in record] )
  if len(groups) < 2:
    return 0
  return sum([len(set(zip([x[1][i:i+3] for x in record], x[2] for x in record]))) == len(groups) for in range(0,ls,3)])

def numNonSynDiv( record, codonTable=CodonTable.standard_dna_table ):
  ls = len(record[0][1])
  groups = set( x[2] for x in record )
  if len(groups) < 2:
    return 0
  codonGroupPairs = [ set( zip( [x[2] for x in record], x[1][i:i+3] for x in record] ) ) for i in range(0,ls,3) ]
  return sum( [ (len(x) == len(groups) and len(set([y[1] for y in x])) != 1) for x in codonGroupPairs ] )

def numSynDiv( record, codonTable=CodonTable.standard_dna_table ):
  return numDiv(record) - numNonSynDiv(record, codonTable)


