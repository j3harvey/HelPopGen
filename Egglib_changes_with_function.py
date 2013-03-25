#################################################
#Script takes a set of genes in fasta format, pipes them into Egglib, and writes population genetic statistics to results file

## Defining variables
outgroup_sp = 'erato' # Specify outgroup + Define variables to assign population numbers. N.B. This does not apply to the pairwise analysis, for which each individual is given a separate population number
number_ingroup_seqs = 16 # We have 4 amaryllis and 4 aglaope individuals each with 2 alleles
number_outgroup_seqs = 8
low_freq_poly_threshold = 0 #  thresholds 0, 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5 # Define low frequency polymorphism removal threshold (All polymorphisms with frequency less than or equal to this value are replaced by major variant at that site. Set at 0 if don't want any removed
missing_data_threshold = 0.5 # Define proportion of missing data a sequence must have before it is chucked out.

## Input file
i = 'ag5am5par1era4_allGenes_StopCodonsRemoved.fasta' # input file

## Results file
Results_file = open ('ag5am5par1era4_outgroup-' + outgroup_sp + '_allGenes_LowFreqThreshold-' + str(low_freq_poly_threshold) + '_StopCodonsRemoved_RESULTS.txt', 'w')
Results_file.write ('Coordinates\tAlleles_ingroup\tLn_ingroup\tLs_ingroup\tPn_ingroup\tPs_ingroup\tFDn aglaope/amaryllis ingroup\tFDs aglaope/amaryllis ingroup\tpairwiseDn\tpairwiseDs\tPiN_ingroup\tPiS_ingroup\tTajima\'s D_ingroup\tFu and Li\'s F\tAlleles_outgroup\tLn_outgroup\tLs_outgroup\tPn_outgroup\tPs_outgroup\tPiN_outgroup\tPiS_outgroup\tTajima\'s D_outgroup\tAlpha- aglaope/amaryllis ingroup\n')

##################################################
### Main body of script
## Defining functions
def AlignByGroupNumbers(align,groupNumbers): # Creates new alignment containing only desired populations
  newAlign = align.slice(0,0)
  for seqNumber in range(len(align)):
    if align[seqNumber][2] in groupNumbers:
      newAlign.addSequences([align[seqNumber]])
  return newAlign

def unique(things): # Returns a list of the unique elements of the argument
  output = []
  for x in things:
    if x not in output:
      output.append(x)
  output.sort()
  return output
  
def inbase_colFreqs(align, columnNumber): # As Simon Martin's colFreqs function but only counting frequency for ingroup
  output = {}
  bases = align.column(columnNumber)
  inBases = []
  for seqNumber in range(len(align)):
    if align[seqNumber][2] != 999:
      inBases.append(bases[seqNumber])
  Acount = float(inBases.count("A"))
  Ccount = float(inBases.count("C"))
  Gcount = float(inBases.count("G"))
  Tcount = float(inBases.count("T"))
  total = Acount + Ccount + Gcount + Tcount
  output["A"] = Acount/total
  output["C"] = Ccount/total
  output["G"] = Gcount/total
  output["T"] = Tcount/total
  return output

## This function takes an alignment with ingroup (1) and outgroup (999) populations defined, and removes polymorphisms below a certain threshold frequency from the ingroup population, replacing them with the major variant at that locus.
# The function also creates a text file that notes down all the changes and the positions that have been made
# Returns an alignment consisting of only ingroup and outgroup populations with low freq polys replaced 
def low_freq_removal (align, low_freq_poly_threshold, ingroup):
  Ingroup_alignment = AlignByGroupNumbers (align, [1])
  Ingroup_outgroup_alignment = AlignByGroupNumbers (align, [1,999]) #Create alignment containing only ingroup and outgroup individuals - not any of the sequences with group no. 0
  Polymorphic_sites = Ingroup_alignment.polymorphism () ['siteIndices'] 
  # Loop across each polymorphic site, removing any base below threshold frequency and replacing it with the major variant at that locus
  #Low_freq_removal_results = open (ingroup + '_Genes_with_low_freq_polys_replaced.txt', 'w')
  #Low_freq_removal_results.write ('Gene\tPosition in sequence\tMinor allele\tMinor allele frequency\tMajor allele\n')
  for position in Polymorphic_sites: # Loop carried out for each polymorphic site along the length of the sequence
    base_freqs = []
    low_freq_base_ID = [] 
    high_freq_base_ID = []
    freq_dict = inbase_colFreqs (Ingroup_alignment, position) # Creates dictionary of base freqs at site for ingroup
    for base in ["A", "C", "G", "T"]:
      if freq_dict[base] != 0: 
        base_freqs.append (freq_dict[base]) # Appends frequencies of each base at given site to list. Only appends bases with non-zero frequency
    if len (base_freqs) > 2:
      print ('triallelic site')
      #Low_freq_removal_results.write (seq_position + '\t' + str(position) + '\ttriallelic site\n')
      continue  #If encounter a triallelic site at which more than 2 bases have a non-zero frequency then end this loop proceed to next iteration.
    if min(base_freqs) <= low_freq_poly_threshold: # carry out following commands if frequency of lowest frequency base at site is less than or equal to the threshold frequency
      # ID low freq base
      for base in ["A", "C", "G", "T"]:
        if freq_dict[base] == min(base_freqs):
          low_freq_base_ID.append (base)
      # ID high freq base
      for base in ["A", "C", "G", "T"]:
        if freq_dict[base] == max(base_freqs):
          high_freq_base_ID.append (base)
      #Produce list of bases for each sequence at given site
      col_align = Ingroup_alignment.column (position) 
      # Replace low freq base with high freq base in list
      for x in range(len(col_align)):
        if low_freq_base_ID[0] in col_align[x]:
          col_align[x] = high_freq_base_ID [0]
      for i in range (Ingroup_alignment.ns()): # Do this for each sequence in the ingroup
        Ingroup_outgroup_alignment.set (i, position, col_align[i]) # i specifies which sequence, position specifies which column, seq[i] specifies which element of list col align should replace with
      print ('Replaced ' + low_freq_base_ID[0] + ' with ' + high_freq_base_ID [0])
      #Low_freq_removal_results.write (seq_position + '\t' + str(position) + '\t' + low_freq_base_ID[0] + '\t' + str(freq_dict[low_freq_base_ID[0]]) + '\t' + high_freq_base_ID [0] + '\n') #Write details of replacements to file- gene, position, low freq base, low freq base frequency, high freq base replacement
  #Low_freq_removal_results.close()
  return Ingroup_outgroup_alignment

def remove_stop_codons (align):
  for x in range (align.ls(), 3):  
    position1 = align.column (x)
    position2 = align.column (x+1)
    position3 = align.column (x+2)
    codon = position1 + position2 + position3
    if codon in ['TAA', 'TGA', 'TAG']:
      position1 = 'N'
      position2 = 'N'
      position3 = 'N'
      for i in range (align.ns()):
        align.set (i, x, position1)
        align.set (i, x+1 , position2)
        align.set (i, x+2, position3)
  return align  
  
import os
import egglib 


### Making a list in which each element is a gene locus
input = open (i, 'rU')
genelist = input.readlines() #Reads all the lines of the file into the list genelist, with each line being a new element
single_element_genelist = "".join (genelist) #Joins the elements of the list genelist into a single element. Call this new list consisting of 1 element single_element_genelist
separate_seqs = single_element_genelist. split ("\n\n") #Splits the single element in the list single_element_genelist into multiple elements by presence of \n\n. This should produce each gene locus as a separate element in list

### The loop for each gene locus in fasta file

## Results file
Results_file = open ('ag5am5par1era4_outgroup-' + outgroup_sp + '_allGenes_LowFreqThreshold-' + str(low_freq_poly_threshold) + '_StopCodonsRemoved_RESULTS.txt', 'w')
Results_file.write ('Coordinates\tAlleles_ingroup\tLn_ingroup\tLs_ingroup\tPn_ingroup\tPs_ingroup\tFDn aglaope/amaryllis ingroup\tFDs aglaope/amaryllis ingroup\tpairwiseDn\tpairwiseDs\tPiN_ingroup\tPiS_ingroup\tTajima\'s D_ingroup\tFu and Li\'s F\tAlleles_outgroup\tLn_outgroup\tLs_outgroup\tPn_outgroup\tPs_outgroup\tPiN_outgroup\tPiS_outgroup\tTajima\'s D_outgroup\tAlpha- aglaope/amaryllis ingroup\n')

for n in separate_seqs:
  if len (n) > 0:
    # Extracting name of locus from file
    seq_lines = n.split ('\n') # n is an element of list separate_seqs. Each element is a single line consisting of genetic info for given locus for all 6 individuals. Here we make a new list seq_lines for which each element is a line like in the original multifasta file. This allows us to get at the name of the file in subsequent steps
    seq_name = "".join (seq_lines [0].split ('_') [0]) #Having produced a list which has as elements each line (as orignally presented) of a single locus. We can now extract the name of the file. Split on underscore. Take all but first and second element. Then join elements with an underscore to restore original name format.
    seq_name_stripped = seq_name.strip ('>') #Remove the '>' from the start of the line
    # This control statement takes into account the different name formatting of GR genes compared to the others (GR gene names have an extra '_' in them)
    if 'HmGr' in seq_lines [0]:
      seq_position = "_".join (seq_lines [0].split ('_') [2:]) #Contains only the positional info for the locus (not which individual)
    else:
      seq_position = "".join (seq_lines [0].split ('_') [2:])
    # Writing new fasta file containing single locus
    output = open (seq_name_stripped + '.fasta', 'w') # Create a write file with the name "seq_position + .fasta". Here we have made this file into a python object.
    output.write (n) #Write n (which represents all info from single locus for all individs) to output file. This loop will repeat itself, each time writing to a different output file depending on its name.    
    output.close ()

    # Identifying individuals with large amounts of missing data and appending there names to a list 
    gene_sequences = open (seq_name_stripped + '.fasta', 'rU')
    gene_sequences_lines = gene_sequences.readlines() # Creates list in which each line is an element
    Missing_data_individuals = []
    Missing_data_individuals_names =[]
    for i in  range (1, len (gene_sequences_lines), 2): # picks out all the lines with sequence information
      if float (gene_sequences_lines [i].count ('N'))/len (gene_sequences_lines [i]) > missing_data_threshold: 
        individual_line = gene_sequences_lines [i-1]
        individual_name = individual_line.strip ('>\n')
        individual_sp = individual_name.split ('.')[0]
        Missing_data_individuals.append (individual_name) #Append any individuals who have more than a threshold fraction of Ns in the sequence to the Missing_data_individuals list. The full name and position is appended to this list so it can later be removed
        Missing_data_individuals_names.append (individual_sp) # The specific name of taxon only is appended to this list

    # Count number of ingroup and outgroup sequences with too much missing data
    ingroup_missing_count = 0
    outgroup_missing_count = 0
    for i in Missing_data_individuals_names:
      if 'amaryllis' in i:
        ingroup_missing_count += 1
      elif 'aglaope' in i:
        ingroup_missing_count += 1
      elif outgroup_sp in i:
        outgroup_missing_count += 1
    
    # Loop is stopped and the next iteration started if: 1) All ingroup sequences have insufficient coverage 2) All outgroup seqs with insufficient coverage 3) Sequences are not of same length
    if ingroup_missing_count >= number_ingroup_seqs - 1: #i.e. if only 1 or 0 ingroup sequences remain. Cannot calculate polymorphism data with just 1 ingroup sequence
      print ('Ingroup individuals have insufficient coverage on ' + seq_position )
      Results_file.write (seq_position + '- Ingroup individuals have insufficient coverage\n')
      continue #Continues to next iteration of loop
    elif outgroup_missing_count >= number_outgroup_seqs - 1: # if only 1 or 0 outgroup sequences have sufficient coverage
      print ('Outgroup individuals have insufficient coverage on ' + seq_position)
      Results_file.write (seq_position + '- Outgroup individuals have insufficient coverage\n')
      continue
       
    ## If none of the above conditions about missing sequences are met the loop continues to calculate pop gen statistics
    # Convert fasta file into egglib object
    my_alignment1 = egglib.Align (seq_name_stripped + '.fasta')
    
    ## Assigning population numbers
    if outgroup_sp == 'pardalinus':
      a = 0 #a is the population number of melpomene
      b = 1 #b is pop no. of aglaope
      c = 1 #c               amaryllis
      d = 999 #d           pardalinus
      e = 0  # e is pop no. of erato
    elif outgroup_sp == 'erato':
      a = 0 #a is the population number of melpomene
      b = 1 #b is pop no. of aglaope
      c = 1 #c               amaryllis
      d = 0 #d           pardalinus
      e = 999 # e is pop no. of erato  
        
    #Create dictionary of population names and their associated population number
    popDict = {"Hmel1-1_1_" + seq_position : a, "Hmel1-1_2_" + seq_position : a, "aglaope.108_1_" + seq_position : b, "aglaope.108_2_" + seq_position : b, "aglaope.572_1_" + seq_position : b, "aglaope.572_2_" + seq_position : b,"aglaope.112_1_"+ seq_position : b, "aglaope.112_2_"+ seq_position : b,"aglaope.569_1_"+ seq_position : b, "aglaope.569_2_"+ seq_position : b, "aglaope.1_1_" + seq_position : 0, "aglaope.1_2_" + seq_position : 0, "amaryllis.216_1_"+ seq_position : c, "amaryllis.216_2_"+ seq_position : c, "amaryllis.160_1_"+ seq_position : c, "amaryllis.160_2_"+ seq_position : c, "amaryllis.48_1_"+ seq_position : 1, "amaryllis.48_2_"+ seq_position : c, "amaryllis.293_1_"+ seq_position : c, "amaryllis.293_2_"+ seq_position : c, "amaryllis.1_1_"+ seq_position : 0, "amaryllis.1_2_"+ seq_position : 0, "pardalinus.371_1_" + seq_position : d, "pardalinus.371_2_" + seq_position : d, "erato.2979_1_" + seq_position: e, "erato.2979_2_" + seq_position: e, "erato.618_1_" + seq_position: e, "erato.618_2_" + seq_position: e, "erato.2980_1_" + seq_position: e, "erato.2980_2_"+ seq_position: e, "erato.2981_1_" + seq_position: e, "erato.2981_2_" + seq_position: e}

    #Assign population number to each individual. Have included this exception handling because 1 gene: separate_seqs[329] produces KeyError when try to assign numbers to individuals. I cannot see why this is so for the moment have allowed script to deal with error
    for i in my_alignment1: # For each individual in alignment:
      i.group = popDict[i.name] # i.group assigns a group number to each i. - i.name specifies the name of i
    
    ## Removing low frequency polymorphisms from ingroup (aglaope/amaryllis)
    Ingroup_outgroup_alignment = low_freq_removal (my_alignment1, low_freq_poly_threshold, 'aglaope_amaryllis') # This function is defined at start of script
    remove_stop_codons (Ingroup_outgroup_alignment)
    
    ## Removing individuals from the alignment who have too much missing data
    if len (Missing_data_individuals) > 0:
      for i in Missing_data_individuals:
        if i in Ingroup_outgroup_alignment.names(): # Removal only carried out if sequence with missing data is in the alignment - this prevents sequences with pop. no. = 0, which are absent from this alignment from having an attempt at being removed
          Ingroup_outgroup_alignment.remove (i)   # Remove sequences with too much missing data from my_alignment1
      
    # Create new alignments with desired populations - N.B. these new alignments are formed from Ingroup_outgroup_alignment which has already had any low freq polys removed
    Align_ingroup_pop = AlignByGroupNumbers(Ingroup_outgroup_alignment, [1]) # Contains only ingroup pop
    Align_outgroup_pop = AlignByGroupNumbers (Ingroup_outgroup_alignment, [999]) # Contains only outgroup pop
    # Re-assign group from outgroup only alignment to 1 - This allows us to use the methods .polymorphism and .polymorphismBPP on this alignment
    for n in range (Align_outgroup_pop.ns()):
      Align_outgroup_pop.group(n, group = 1)
    
    # Calculate statistics for each gene locus
    try:
      data1 = Align_ingroup_pop.polymorphism (minimumExploitableData = 1) #Set minimumExploitableData to 1 so that are including same amount of info in polymorphism and polymorphismBPP
      data1BPP = Align_ingroup_pop.polymorphismBPP(dataType = 4) #DataType = 4 means that the data is standard codons
      data2 = Ingroup_outgroup_alignment.polymorphism (minimumExploitableData = 1)
      data2BPP = Ingroup_outgroup_alignment.polymorphismBPP (dataType = 4)
      data_outgrp = Align_outgroup_pop.polymorphism(minimumExploitableData = 1)
      data_outgrpBPP = Align_outgroup_pop.polymorphismBPP (dataType = 4)
    except RuntimeError as error:
      Results_file.write (seq_position + '- Problem in using polymorphism methods, ingroup with insufficient data -' + str(error) + '\n') # 
      continue
     
    # Calculate alpha statistic from MK test results table
    if ((float(data2BPP ['MK'][2]))*(float(data2BPP ['MK'][1]))) != 0: # Only calculate if the denominator != 0
      Alpha = float (1) - ((float((data2BPP ['MK'][3]))*(float(data2BPP ['MK'][0])))/((float(data2BPP ['MK'][2]))*(float(data2BPP ['MK'][1]))))
    else:
      Alpha = 'NA'    
    
    ## With erato as outgroup
    # Reassign population numbers
    a = 0
    b = 999 #aglaope pop no.
    c = 999 #amaryllis pop no.
    d = 0 #pardalinus pop no.
    e = 1 #erato pop no.
        
    New_outgroup_alignment = my_alignment1 # This alignment will be for aglaope/amaryllis as outgroup
    
    #Create dictionary of population names and their associated population number
    popDict = {"Hmel1-1_1_" + seq_position : a, "Hmel1-1_2_" + seq_position : a, "aglaope.108_1_" + seq_position : b, "aglaope.108_2_" + seq_position : b, "aglaope.572_1_" + seq_position : b, "aglaope.572_2_" + seq_position : b,"aglaope.112_1_"+ seq_position : b, "aglaope.112_2_"+ seq_position : b,"aglaope.569_1_"+ seq_position : b, "aglaope.569_2_"+ seq_position : b, "aglaope.1_1_" + seq_position : b, "aglaope.1_2_" + seq_position : b, "amaryllis.216_1_"+ seq_position : c, "amaryllis.216_2_"+ seq_position : c, "amaryllis.160_1_"+ seq_position : c, "amaryllis.160_2_"+ seq_position : c, "amaryllis.48_1_"+ seq_position : 1, "amaryllis.48_2_"+ seq_position : c, "amaryllis.293_1_"+ seq_position : c, "amaryllis.293_2_"+ seq_position : c, "amaryllis.1_1_"+ seq_position : c, "amaryllis.1_2_"+ seq_position : c, "pardalinus.371_1_" + seq_position : d, "pardalinus.371_2_" + seq_position : d, "erato.2979_1_" + seq_position: e, "erato.2979_2_" + seq_position: e, "erato.618_1_" + seq_position: e, "erato.618_2_" + seq_position: e, "erato.2980_1_" + seq_position: e, "erato.2980_2_"+ seq_position: e, "erato.2981_1_" + seq_position: e, "erato.2981_2_" + seq_position: e}

    #Assign population number to each individual. Have included this exception handling because 1 gene: separate_seqs[329] produces KeyError when try to assign numbers to individuals. I cannot see why this is so for the moment have allowed script to deal with error
    for i in New_outgroup_alignment: # For each individual in alignment:
      i.group = popDict[i.name] # i.group assigns a group number to each i. - i.name specifies the name of i
           
    # Removing low frequency polymorphisms from ingroup (now erato is ingroup)
    New_ingroup_outgroup_alignment = low_freq_removal (New_outgroup_alignment, low_freq_poly_threshold, 'erato')
    remove_stop_codons (New_ingroup_outgroup_alignment)
    # Calculate statistics with aglaope/amaryllis outgroup
    MK_test_table = New_ingroup_outgroup_alignment.polymorphismBPP (dataType = 4) ['MK']
    outgroup_PN = MK_test_table [0]
    outgroup_PS = MK_test_table [1]
    outgroup_DN = MK_test_table [2]
    outgroup_DS = MK_test_table [3]
    
    ## Obtaining average pairwise distance between sequences 
    my_alignment2 = egglib.Align (seq_name_stripped + '.fasta')
    popDict_pairwise = {"Hmel1-1_1_" + seq_position : 0, "Hmel1-1_2_" + seq_position : 0, "aglaope.108_1_" + seq_position :1, "aglaope.108_2_" + seq_position :2, "aglaope.572_1_" + seq_position :3, "aglaope.572_2_" + seq_position :4,"aglaope.112_1_"+ seq_position :5, "aglaope.112_2_"+ seq_position :6,"aglaope.569_1_"+ seq_position :7, "aglaope.569_2_"+ seq_position :8, "aglaope.1_1_" + seq_position :0, "aglaope.1_2_" + seq_position :0, "amaryllis.216_1_"+ seq_position :9, "amaryllis.216_2_"+ seq_position :10, "amaryllis.160_1_"+ seq_position :11, "amaryllis.160_2_"+ seq_position :12, "amaryllis.48_1_"+ seq_position :13, "amaryllis.48_2_"+ seq_position :14, "amaryllis.293_1_"+ seq_position :15, "amaryllis.293_2_"+ seq_position :16, "amaryllis.1_1_"+ seq_position :0, "amaryllis.1_2_"+ seq_position :0, "pardalinus.371_1_" + seq_position :0, "pardalinus.371_2_" + seq_position :0, "erato.2979_1_" + seq_position: 17, "erato.2979_2_" + seq_position: 18, "erato.618_1_" + seq_position: 19, "erato.618_2_" + seq_position: 20, "erato.2980_1_" + seq_position: 21, "erato.2980_2_"+ seq_position: 22, "erato.2981_1_" + seq_position: 23, "erato.2981_2_" + seq_position: 24}
    for i in my_alignment2: # For each individual in alignment:
      i.group = popDict_pairwise[i.name] # i.group assigns a group number to each i. - i.name specifies the name of i
    
    if outgroup_sp == 'pardalinus':
      outgroups = [] #Need to enter in pop numbers of pardalinus if include this
    elif outgroup_sp == 'erato':
      outgroups = range (17,25) #creates list of numbers from 17 to 24 inclusive
    
    ingroups = range (1,17) #creates list of numbers from 1 to 16 inclusive
    Non_syn_polys = []
    Syn_polys = []
    # Calculate Dn and Ds between 2 sequences by placing the 2 sequences together in an alignment and calculating no. of syn and non syn polymorphic sites
    for outgroup in outgroups:
      for ingroup in ingroups:
        pairwise_alignment = AlignByGroupNumbers(my_alignment2, [outgroup, ingroup])
        data_pairwiseBPP = pairwise_alignment.polymorphismBPP(dataType = 4)
        Non_syn_polys.append (data_pairwiseBPP ['SNS'])
        Syn_polys.append (data_pairwiseBPP ['SS'])
    Non_syn_polys_avg = float (sum(Non_syn_polys))/ float (len(Non_syn_polys))
    Syn_polys_avg = float (sum(Syn_polys))/float(len(Syn_polys))

    # Write statistics to results file- Use BPP for certain stats, use 1 when only ingroup required, 2 when outgroup also required
    Results_file.write (seq_position + '\t' + str(Align_ingroup_pop.ns()) + '\t' + str(data1BPP['NSsites']) + '\t' + str(data1BPP['Ssites']) + '\t' + str(data2BPP['MK'][0]) + '\t' + str(data2BPP['MK'][1]) + '\t' + str(data2BPP ['MK'][2]) + '\t' + str(data2BPP ['MK'][3]) + '\t' + str(Non_syn_polys_avg) + '\t' + str(Syn_polys_avg) + '\t' + str(data1BPP ['PiNS']) + '\t' + str(data1BPP['PiS']) + '\t' + str(data1 ['D']) + '\t' +  str(data2BPP['F'])  + '\t' + str(Align_outgroup_pop.ns()) + '\t' + str(data_outgrpBPP['NSsites']) + '\t' + str(data_outgrpBPP['Ssites']) + '\t' + str(outgroup_PN) + '\t' + str(outgroup_PS) + '\t' + str(data_outgrpBPP ['PiNS']) + '\t' + str(data_outgrpBPP['PiS']) + '\t' + str(data_outgrp ['D']) + '\t' + str(Alpha) + '\n')
    # Delete the fasta file that has been created for this gene locus alone - Prevents accumulation of large number of intermediary fasta files
    os.system ('rm ' + seq_name_stripped + '.fasta')
  
Results_file.close ()

  
  
# 'Coordinates\tAlleles_ingroup\tLn_ingroup\tLs_ingroup\tPn_ingroup\tPs_ingroup\tFDn aglaope/amaryllis ingroup\tFDs aglaope/amaryllis ingroup\tpairwiseDn\tpairwiseDs\tPiN_ingroup\tPiS_ingroup\tTajima\'s D_ingroup\tFu and Li\'s F\tAlleles_outgroup\tLn_outgroup\tLs_outgroup\tPn_outgroup\tPs_outgroup\tFDN erato ingroup\tFDS erato ingroup\tPiN_outgroup\tPiS_outgroup\tTajima\'s D_outgroup\tAlpha\n'




##NB The order in which the MK test result is given is (Pn, Ps, Dn, Ds)

## expand names to take into account new individuals  + set inbred individuals population numbers to 0