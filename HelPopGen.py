#! /usr/bin/env python2.7


'''
This script takes a set of genes in fasta format, and calculates
population genetics statistics using egglib.
'''


'''
CHANGES TODO:

 1 add if __name__ == "__main__" idiom                    DONE 15/3/13
 2 add !#                                                 DONE 15/3/13
 3 make the logging pythony: use csv module
 4 fix locus name extraction
 5 triallelic site handling
 6 remove references fo species names
   only distinguish ingrup/outgroup
 7 fix stop codon removal
 8 rewrite everything

'''

## Import dependencies

import os
import egglib

## Defining variables

# Specify outgroup + Define variables to assign population numbers.
# N.B. This does not apply to the pairwise analysis, for which each individual
# is given a separate population number
outgroup_sp = 'erato'
number_ingroup_seqs = 16 # We have 4 amaryllis and 4 aglaope individuals each with 2 alleles
number_outgroup_seqs = 8

# A case-insensitive list of stop codons
StopCodons = set( ['TAA', 'TGA', 'TAG'] )

# Define low frequency polymorphism removal threshold
# All polymorphisms with frequency less than or equal to this value are
# replaced by major variant at that site.
# thresholds 0, 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5
# Set at 0 if don't want any removed
low_freq_poly_threshold = 0

# Define proportion of missing data a sequence must have before it is chucked out.
missing_data_threshold = 0.5

# Input file
Input_file = 'test.fasta'

# Results file (#3)
Results_file = open ('ag5am5par1era4_outgroup-' + 
                     outgroup_sp + '_allGenes_LowFreqThreshold-' + 
                     str(low_freq_poly_threshold) + 
                     '_StopCodonsRemoved_RESULTS.txt', 'w')

ColHeaders = ['Coordinates',            # locus
              'Alleles_ingroup',        # number of ingroup alleles
              'Ln_ingroup',             # number of non-synonymous sites
              'Ls_ingroup',             # number of sunonymous sites
              'Pn_ingroup',             # number of non-synonymous polymorphisms
              'Ps_ingroup',             # number of synonymous polymorphisms
              'FDn',                    # number of non-synonymous divergences
              'FDs',                    # number of synonymous divergences
              'pairwiseDn',             # average pairwise non-synonymous divergence -- BORING
              'pairwiseDs',             # average pairwise synonymous divergence -- BORING
              'PiN_ingroup',            # ingroup non-synonymous heterozygocity
              'PiS_ingroup',            # ingroup synonymous heterozygocity
              'Alleles_outgroup',       # number of outgroup alleles
              'Ln_outgroup',            # number of non-synonymous sites
              'Ls_outgroup',            # number of synonymous sites
              'Pn_outgroup',            # number of non-synonymous polymorphisms
              'Ps_outgroup',            # number of synonymous polymorphisms
              'PiN_outgroup',           # non-synonymous heterozygocity
              'PiS_outgroup']           # non-synonymous heterozygocity
Results_file.write( '\t'.join(ColHeaders) + '\n' )



##################################################
## Defining functions

def AlignByGroupNumbers(align,groupNumbers):
  '''Creates new alignment containing only desired populations'''
  newAlign = align.slice(0,0)
  for seqNumber in range(len(align)):
    if align[seqNumber][2] in groupNumbers:
      newAlign.addSequences([align[seqNumber]])
  return newAlign

#def inbase_colFreqs(align, columnNumber):
#  '''Only used in low_freq_removal so can probably be removed'''
#  bases = [ x for (i,x) in enumerate( align.column(columnNumber) ) 
#              if align[i][2] != 999          # ignore outgroup sequences
#              and x not in "Nn- " ]          # only count valid bases
#  return { "A" : bases.count("A") / len( inBases ),
#           "C" : bases.count("C") / len( inBases ),
#           "G" : bases.count("G") / len( inBases ),
#           "T" : bases.count("T") / len( inBases ) }


def low_freq_removal (align, low_freq_poly_threshold, ingroup):
  '''
  This function takes an alignment with ingroup (1) and outgroup (999) 
  populations defined, and removes polymorphisms below a certain threshold 
  frequency from the ingroup population, replacing them with the major variant 
  at that locus.
  Returns an alignment consisting of only ingroup and outgroup populations 
  with low freq polymorphisms replaced 
  '''
  Ingroup_alignment = AlignByGroupNumbers (align, [1])
  Polymorphic_sites = Ingroup_alignment.polymorphism () ['siteIndices']
  #Create alignment containing only ingroup and outgroup individuals
  Ingroup_outgroup_alignment = AlignByGroupNumbers (align, [1,999])
  
  # Loop across each polymorphic site, removing any base below threshold
  # frequency and replacing it with the major variant at that locus
  for position in Polymorphic_sites:
    freq_dict = inbase_colFreqs (Ingroup_alignment, position) 
    bases = [ x for x in Ingroup_alignment.column(position) if x not in "Nn- " ]
    freq_dict =  { "A" : bases.count("A") / len( inBases ),
                   "C" : bases.count("C") / len( inBases ),
                   "G" : bases.count("G") / len( inBases ),
                   "T" : bases.count("T") / len( inBases ) }
    base_freqs = [ x for x in freq_dict.values() if x > 0 ]
    if len (base_freqs) > 2:      #
      print ('triallelic site')   # -- TRIALLELIC SITE HANDLING (#5)
      continue                    #
    low_freq_base_ID = [ x for x in "ACTG" if freq_dict[x] == min(base_freqs) ]
    high_freq_base_ID = [ x for x in "ACTG" if freq_dict[x] == max(base_freqs) ]
    newCol = [ high_freq_base_ID[0] if x in low_freq_base_ID else x for x in align.column(position) ]
    for i in range (Ingroup_alignment.ns()):
        # i specifies which sequence
        # position specifies which column
        # seq[i] specifies which element of list col align should replace with
        Ingroup_outgroup_alignment.set (i, position, newCol[i])
    return Ingroup_outgroup_alignment


def remove_stop_codons (align):
  '''This does not do what it claims'''
  #for x in range (align.ls(), 3):
  #  position1 = align.column (x)
  #  position2 = align.column (x+1)
  #  position3 = align.column (x+2)
  #  codon = position1 + position2 + position3
  #  if codon in ['TAA', 'TGA', 'TAG']:
  #    position1 = 'N'
  #    position2 = 'N'
  #    position3 = 'N'
  #    for i in range (align.ns()):
  #      align.set (i, x, position1)
  #      align.set (i, x+1 , position2)
  #      align.set (i, x+2, position3)
  return align

##################################################
### Main body of script

def main():
  '''
  The main body of the script is run when the script is called from a shell.
  If the script is imported from within python this function is not run
  but appears as main()
  '''
  
  ### Making a list in which each element is a gene locus
  with open( Input_file, 'rU') as f:
    separate_seqs = "".join( f.readlines() ).split("\n\n")
  
  ### The loop for each gene locus in fasta file
  
  for locus in separate_seqs:
    # quality control:
    #   choose maximinally informative subset of sequences
    #   remove low frequency polymorphisms
    #   phase data taking into account multi-change codons
    #   remove stops
    #   calculate stats
    #   write results
  
  for n in separate_seqs:
    if len (n) > 0:
      ## Extracting name of locus from file -- THIS DOESN'T SEEM TO WORK (#4)
      
      # Each seq_line is a single line consisting of genetic info for given locus for all individuals.
      # Here we make a new list seq_lines for which each element is a line like in the original multifasta file.
      # This allows us to get at the name of the file in subsequent steps
      seq_lines = n.split ('\n')
      
      # We can now extract the name of the file. Split on underscore.  |
      # Take all but first and second element. Then join elements with |-- THIS IS NOT AT ALL WHAT IS HAPPENING (#4)
      # an underscore to restore original name format.                 |
      seq_name = "".join (seq_lines [0].split ('_') [0])
      
      # Remove the '>' from the start of the line -- THIS SHOULD NOT BE NECESSARY (#4)
      seq_name_stripped = seq_name.strip ('>') 
      
      # Contains only the positional info for the locus (not which individual)
      # This control statement takes into account the different name formatting 
      # of GR genes compared to the others (GR gene names have an extra '_' in them)
      if 'HmGr' in seq_lines [0]:
        seq_position = "_".join (seq_lines [0].split ('_') [2:])
      else:
        seq_position = "".join (seq_lines [0].split ('_') [2:])
      
      # Writing new fasta file containing all info from single locus
      with open (seq_name_stripped + '.fasta', 'w') as f:
        f.write(n)
  
      # Identifying individuals with large amounts of missing data and 
      # appending their names to a list 
      with open (seq_name_stripped + '.fasta', 'rU') as gene_sequences:
        gene_sequences_lines = gene_sequences.readlines()
        Missing_data_individuals = []
        Missing_data_individuals_names =[]
        # picks out all the lines with sequence information
        for i in  range (1, len (gene_sequences_lines), 2):
          if float (gene_sequences_lines [i].count ('N'))/len (gene_sequences_lines [i]) > missing_data_threshold: 
            individual_line = gene_sequences_lines [i-1]
            individual_name = individual_line.strip ('>\n')
            individual_sp = individual_name.split ('.')[0]
            # Append any individuals who have more than a threshold fraction of Ns
            # in the sequence to the Missing_data_individuals list.
            Missing_data_individuals.append (individual_name)
            Missing_data_individuals_names.append (individual_sp)
  
      # Count number of ingroup and outgroup sequences with too much missing data
      # NO use group numbers...
      ingroup_missing_count = 0
      outgroup_missing_count = 0
      for i in Missing_data_individuals_names:
        if 'amaryllis' in i:
          ingroup_missing_count += 1
        elif 'aglaope' in i:
          ingroup_missing_count += 1
        elif outgroup_sp in i:
          outgroup_missing_count += 1
      
      # Loop is stopped and the next iteration started if:
      #  1) Only 0 or 1 ingroup sequences have sufficient coverage
      #  2) Only 0 or 1 outgroup seqs have sufficient coverage
      #  3) Sequences are not of same length
      
      # Stop if only 1 or 0 ingroup sequences have sufficient coverage. 
      # Cannot calculate polymorphism data with just 1 sequence
      if ingroup_missing_count >= number_ingroup_seqs - 1: 
        print('Ingroup individuals have insufficient coverage on ' + seq_position )
        Results_file.write(seq_position + 
                            ' - Ingroup individuals have insufficient coverage\n') # (#3)
        continue
      # Stop if only 1 or 0 outgroup sequences have sufficient coverage
      elif outgroup_missing_count >= number_outgroup_seqs - 1:
        print('Outgroup individuals have insufficient coverage on ' + seq_position)
        Results_file.write(seq_position + 
                            ' - Outgroup individuals have insufficient coverage\n') # (#3)
        continue
         
      ## Calculate pop gen statistics

      # Convert fasta file into egglib object
      try:
        my_alignment1 = egglib.Align (seq_name_stripped + '.fasta')
      except ValueError:
        print( "Alignment failed at " + seq_position )
        Results_file.write(seq_position + " - alignment failed.")
        continue
      
      ## Assigning population numbers
      a = 0           # a is the population number of melpomene
      b = 1           # b                             aglaope
      c = 1           # c                             amaryllis
      d = 0           # d                             pardalinus
      e = 999         # e                             erato  
          
      #Create dictionary of population names and their associated population number
      popDict = {"Hmel1-1_1_" + seq_position : a, 
                 "Hmel1-1_2_" + seq_position : a, 
                 "aglaope.108_1_" + seq_position : b, 
                 "aglaope.108_2_" + seq_position : b, 
                 "aglaope.572_1_" + seq_position : b, 
                 "aglaope.572_2_" + seq_position : b,
                 "aglaope.112_1_"+ seq_position : b, 
                 "aglaope.112_2_"+ seq_position : b,
                 "aglaope.569_1_"+ seq_position : b, 
                 "aglaope.569_2_"+ seq_position : b, 
                 "aglaope.1_1_" + seq_position : 0, # not b?
                 "aglaope.1_2_" + seq_position : 0, # not b?
                 "amaryllis.216_1_"+ seq_position : c, 
                 "amaryllis.216_2_"+ seq_position : c, 
                 "amaryllis.160_1_"+ seq_position : c, 
                 "amaryllis.160_2_"+ seq_position : c, 
                 "amaryllis.48_1_"+ seq_position : 1, # not c?
                 "amaryllis.48_2_"+ seq_position : c, 
                 "amaryllis.293_1_"+ seq_position : c, 
                 "amaryllis.293_2_"+ seq_position : c, 
                 "amaryllis.1_1_"+ seq_position : 0, # not c?
                 "amaryllis.1_2_"+ seq_position : 0, # not c?
                 "pardalinus.371_1_" + seq_position : d, 
                 "pardalinus.371_2_" + seq_position : d, 
                 "erato.2979_1_" + seq_position: e, 
                 "erato.2979_2_" + seq_position: e, 
                 "erato.618_1_" + seq_position: e, 
                 "erato.618_2_" + seq_position: e, 
                 "erato.2980_1_" + seq_position: e, 
                 "erato.2980_2_"+ seq_position: e, 
                 "erato.2981_1_" + seq_position: e, 
                 "erato.2981_2_" + seq_position: e
                }
  
      #Assign population number to each individual.
      for i in my_alignment1:
        i.group = popDict[i.name]
      
      ## Removing low frequency polymorphisms from ingroup (aglaope/amaryllis)
      Ingroup_outgroup_alignment = low_freq_removal (my_alignment1, low_freq_poly_threshold, 'aglaope_amaryllis')
      remove_stop_codons (Ingroup_outgroup_alignment)
      
      ## Removing individuals from the alignment who have too much missing data
      
      for i in Missing_data_individuals:
        # Removal only carried out if sequence with missing data is in the alignment
        if i in Ingroup_outgroup_alignment.names():
          Ingroup_outgroup_alignment.remove (i)
        
      # Create new alignments with desired populations
      # N.B. these new alignments are formed from Ingroup_outgroup_alignment
      # which has already had any low freq polys removed
      Align_ingroup_pop = AlignByGroupNumbers(Ingroup_outgroup_alignment, [1])
      Align_outgroup_pop = AlignByGroupNumbers (Ingroup_outgroup_alignment, [999])
      # Re-assign group from outgroup only alignment to 1 
      # This allows us to use the methods .polymorphism and .polymorphismBPP on this alignment
      for n in range (Align_outgroup_pop.ns()):
        Align_outgroup_pop.group(n, group = 1)
      
      
      ## Calculate statistics for each gene locus
      
      try:
        # Set minimumExploitableData to 1 so that are including same amount
        # of info in polymorphism and polymorphismBPP
        data1 = Align_ingroup_pop.polymorphism (minimumExploitableData = 1) 
        # DataType = 4 means that the data is standard codons
        data1BPP = Align_ingroup_pop.polymorphismBPP(dataType = 4) 
        data2 = Ingroup_outgroup_alignment.polymorphism (minimumExploitableData = 1)
        data2BPP = Ingroup_outgroup_alignment.polymorphismBPP (dataType = 4)
        data_outgrp = Align_outgroup_pop.polymorphism(minimumExploitableData = 1)
        data_outgrpBPP = Align_outgroup_pop.polymorphismBPP (dataType = 4)
      except RuntimeError as error: # (#3)
        Results_file.write (seq_position + 
                            '- Problem in using polymorphism methods, ingroup with insufficient data -' + 
                            str(error) + '\n')
        continue
      
      
      ## With erato as outgroup
      
      # Reassign population numbers
      a = 0           # a is the population number of melpomene
      b = 999         # b                             aglaope
      c = 999         # c                             amaryllis
      d = 0           # d                             pardalinus
      e = 1           # e                             erato
      
      # This alignment will be for aglaope/amaryllis as outgroup
      New_outgroup_alignment = my_alignment1 
      
      #Create dictionary of population names and their associated population number
      popDict = {"Hmel1-1_1_" + seq_position : a,
                 "Hmel1-1_2_" + seq_position : a,
                 "aglaope.108_1_" + seq_position : b,
                 "aglaope.108_2_" + seq_position : b,
                 "aglaope.572_1_" + seq_position : b,
                 "aglaope.572_2_" + seq_position : b,
                 "aglaope.112_1_"+ seq_position : b,
                 "aglaope.112_2_"+ seq_position : b,
                 "aglaope.569_1_"+ seq_position : b,
                 "aglaope.569_2_"+ seq_position : b,
                 "aglaope.1_1_" + seq_position : b,
                 "aglaope.1_2_" + seq_position : b,
                 "amaryllis.216_1_"+ seq_position : c,
                 "amaryllis.216_2_"+ seq_position : c,
                 "amaryllis.160_1_"+ seq_position : c,
                 "amaryllis.160_2_"+ seq_position : c,
                 "amaryllis.48_1_"+ seq_position : 1, # Why is this 1, not c?
                 "amaryllis.48_2_"+ seq_position : c,
                 "amaryllis.293_1_"+ seq_position : c,
                 "amaryllis.293_2_"+ seq_position : c,
                 "amaryllis.1_1_"+ seq_position : c,
                 "amaryllis.1_2_"+ seq_position : c,
                 "pardalinus.371_1_" + seq_position : d,
                 "pardalinus.371_2_" + seq_position : d,
                 "erato.2979_1_" + seq_position: e,
                 "erato.2979_2_" + seq_position: e,
                 "erato.618_1_" + seq_position: e,
                 "erato.618_2_" + seq_position: e,
                 "erato.2980_1_" + seq_position: e,
                 "erato.2980_2_" + seq_position: e,
                 "erato.2981_1_" + seq_position: e,
                 "erato.2981_2_" + seq_position: e
                }
      
      #Assign population number to each individual.
      for i in New_outgroup_alignment:
        i.group = popDict[i.name]
             
      # Removing low frequency polymorphisms from ingroup (now erato is ingroup)
      New_ingroup_outgroup_alignment = low_freq_removal (New_outgroup_alignment, low_freq_poly_threshold, 'erato')
      remove_stop_codons (New_ingroup_outgroup_alignment) # -- ERK!! This is not working
      # Calculate statistics with aglaope/amaryllis outgroup
      MK_test_table = New_ingroup_outgroup_alignment.polymorphismBPP (dataType = 4) ['MK']
      outgroup_PN = MK_test_table [0]
      outgroup_PS = MK_test_table [1]
      outgroup_DN = MK_test_table [2]
      outgroup_DS = MK_test_table [3]
      
      # Write statistics to results file- Use BPP for certain stat
      # use 1 when only ingroup required, 2 when outgroup also required (#3)
      Results_file.write( '\t'.join( [seq_position,
                                      str(Align_ingroup_pop.ns()),
                                      str(data1BPP['NSsites']),
                                      str(data1BPP['Ssites']),
                                      str(data2BPP['MK'][0]),
                                      str(data2BPP['MK'][1]),
                                      str(data2BPP ['MK'][2]),
                                      str(data2BPP ['MK'][3]),
                                      str(data1BPP ['PiNS']),
                                      str(data1BPP['PiS']),
                                      str(Align_outgroup_pop.ns()),
                                      str(data_outgrpBPP['NSsites']),
                                      str(data_outgrpBPP['Ssites']),
                                      str(outgroup_PN),
                                      str(outgroup_PS),
                                      str(data_outgrpBPP ['PiNS']),
                                      str(data_outgrpBPP['PiS'])) + '\n')
      
      # Delete the fasta file that has been created for this gene locus alone
      # Prevents accumulation of large number of intermediary fasta files
      os.system ('rm ' + seq_name_stripped + '.fasta')
  
  Results_file.close ()

#When run from a shell as "python script.py" the function main() is executed
#When imported from within the python interpreter, all functions are loaded but main() is not run
if __name__ == "__main__":
  main()

