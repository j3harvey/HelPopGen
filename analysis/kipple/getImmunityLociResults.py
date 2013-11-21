#! /usr/bin/env python2.7

import csv
from controlGenes import getControlGenes
import matplotlib.pyplot as plt
import numpy.numarray as na
from pylab import *

'''
with open("../data/ImmunityGenes.csv", 'r') as f:
    lines = f.readlines()

ids = [l.split(',')[0] for l in lines]
names = [l.split(',')[1] for l in lines]
dictNames = {i:n for (i, n) in zip(ids, names)}

print ids
print names

count = 0
with open("../results/third_run_00.csv", 'r') as f:
    for l in f.readlines():
        for i in ids:
            if (i + '-') in l:
                (P_N, P_S, D_N, D_S) = l.split(',')[-4:]
                if int(D_N) != 0 and int(P_S) != 0:
                    print dictNames[i], ':  ', 1 - float(P_N)*float(D_S)/(float(P_S)*float(D_N))
                else:
                    print dictNames[i], ':  ', P_N, P_S, D_N, D_S
                count += 1
# 91 immune genes have results
# print count
'''

def aBar(loci):
    with open("../results/results_004314052013.csv", 'r') as f:
        results = [l.split(',') for l in f.readlines() if any([i in l for i in loci])]
        P_N = [float(x[-6]) for x in results]
        P_S = [float(x[-5]) for x in results]
        D_N = [float(x[-2]) for x in results]
        D_S = [float(x[-1]) for x in results]
        L_NS = [float(x[4]) for x in results]
        return sum([(D_N[i] - P_N[i]*D_S[i]/(P_S[i] +1))/L_NS[i] for i in range(len(results))])/len(results)

def alphaBar(loci):
    with open("../results/results_004314052013.csv", 'r') as f:
        results = [l.split(',') for l in f.readlines() if any([i in l for i in loci])]
        P_N = [float(x[-6]) for x in results]
        P_S = [float(x[-5]) for x in results]
        D_N = [float(x[-2]) for x in results]
        D_S = [float(x[-1]) for x in results]
        return 1 - sum([P_N[i] / (P_S[i] + 1) for i in range(len(P_N))])*sum(D_S)/(sum(D_N)*len(P_S))

# load list of immunity genes
with open("../data/hmel_dmel_orthologues.csv", 'r') as f:
    field_names = f.readline().rstrip().split(',')
    immunity_genes = [x for x in csv.DictReader(f, fieldnames=field_names)]

mapped_immunity_genes = mapped([x['Hmel_ID'] for x in immunity_genes])
control_genes = getControlGenes(mapped_immunity_genes)

# load Obbard 2009 gene table
with open('../data/s027_all.csv', 'r') as f:
    field_names = f.readline().rstrip().split(',')
    Obb_genes = [x for x in csv.DictReader(f, fieldnames=field_names)]

# get classes
for i in immunity_genes:
    d = i['Dmel_ID']
    for g in Obb_genes:
        if g['FBgn'] == d:
            i['Class'] = g['Class']
            i['Cell_Hum'] = g['Cell_Hum']

#with open("../data/ImmunityGenes.csv", 'r') as f:
#    rows = [x for x in csv.DictReader(f)]

# all immune genes
immune = [row["Hmel_ID"] for row in immunity_genes
        if row["Class"] != "RNAi"
        and row["Cell_Hum"] != "AMP"
        and len(getControlGenes(row["Hmel_ID"])) > 0]

immune_control = getControlGenes(immune)

#print "immune:                        ", alphaBar(immune), '\t', aBar(immune)
#print "immune_control:                ", alphaBar(immune_control), '\t', aBar(immune_control)
#print ''

# humoral genes
humoral = [row["Hmel_ID"] for row in immunity_genes 
                  if row["Class"] == "humoral"
                  and len(getControlGenes(row["Hmel_ID"])) > 0]
# cellular genes
cellular = [row["Hmel_ID"] for row in immunity_genes 
                  if row["Class"] == "cellular"
                  and len(getControlGenes(row["Hmel_ID"])) > 0]
# RNAi genes
RNAi_class = [row["Hmel_ID"] for row in immunity_genes 
                  if row["Class"] == "RNAi"
                  and len(getControlGenes(row["Hmel_ID"])) > 0]

humoral_control = getControlGenes(humoral)
cellular_control = getControlGenes(cellular)
RNAi_class_control = getControlGenes(RNAi_class)

#
#print "Gene class                      AlphaBar"
#print ''
#print "humoral:                       ", alphaBar(humoral), '\t', aBar(humoral)
#print "humoral_control:               ", alphaBar(humoral_control), '\t', aBar(humoral_control)
#print ''
#print "cellular:                       ", alphaBar(cellular), '\t', aBar(humoral)
#print "cellular_control:              ", alphaBar(cellular_control), '\t', aBar(cellular_control)
#print ''
#print "RNAi_class:                    ", alphaBar(RNAi_class), '\t', aBar(RNAi_class)
#print "RNAi_class_control:            ", alphaBar(RNAi_class_control), '\t', aBar(RNAi_class_control)
#print ''
#print ''

# signalling genes
signalling = [row["Hmel_ID"] for row in immunity_genes 
                  if row["Cell_Hum"] == "signalling"
                  and len(getControlGenes(row["Hmel_ID"])) > 0]
# anti-microbial peptided
amp = [row["Hmel_ID"] for row in immunity_genes 
                  if row["Cell_Hum"] == "AMP"
                  and len(getControlGenes(row["Hmel_ID"])) > 0]
# cellular recognition genes
cellular_recognition = [row["Hmel_ID"] for row in immunity_genes 
                  if row["Cell_Hum"] == "cellular_recognition"
                  and len(getControlGenes(row["Hmel_ID"])) > 0]
# humoral recognition genes
humoral_recognition = [row["Hmel_ID"] for row in immunity_genes 
                  if row["Cell_Hum"] == "humoral_recognition"
                  and len(getControlGenes(row["Hmel_ID"])) > 0]
# RNAi genes
RNAi_cell_hum = [row["Hmel_ID"] for row in immunity_genes 
                  if row["Cell_Hum"] == "RNAi"
                  and len(getControlGenes(row["Hmel_ID"])) > 0]

signalling_control = getControlGenes(signalling)
amp_control = getControlGenes(amp)
cellular_recognition_control = getControlGenes(cellular_recognition)
humoral_recognition_control = getControlGenes(humoral_recognition)
RNAi_cell_hum_control = getControlGenes(RNAi_cell_hum)

#print "Cell hum                        AlphaBar"
#print ''
#print "signalling:                    ", alphaBar(signalling), '\t', aBar(signalling)
#print "signalling_control:            ", alphaBar(signalling_control), '\t', aBar(signalling_control)
#print ''
#print "amp:                           ", alphaBar(amp), '\t', aBar(amp)
#print "amp_control:                   ", alphaBar(amp_control), '\t', aBar(amp_control)
#print ''
#print "cellular_recognition:           ", alphaBar(cellular_recognition), '\t', aBar(cellular_recognition)
#print "cellular_recognition_control:  ", alphaBar(cellular_recognition_control), '\t', aBar(cellular_recognition_control)
#print ''
#print "humoral_recognition:           ", alphaBar(humoral_recognition), '\t', aBar(humoral_recognition)
#print "humoral_recognition_control:   ", alphaBar(humoral_recognition_control), '\t', aBar(humoral_recognition_control)
#print ''
#print "RNAi_cell_hum:                 ", alphaBar(RNAi_cell_hum), '\t', aBar(RNAi_cell_hum_control)
#print "RNAi_cell_hum_control:         ", alphaBar(RNAi_cell_hum_control), '\t', aBar(RNAi_cell_hum_control)
#print ''

print "amps:  ", amp # amps are bulshit
print "RNAi:  ", RNAi_class

# Plot the results

labels = ["Immune", 
          "Control", 
          '', 
          'Cellular', 
          #'Control',
          'Humoral', 
          #'Control',
          #'RNAi', 
          #'Control',
          '', 
          'Signalling', 
          #'Control',
          #'AMP', 
          #'Control',
          'Cellular\nrecognition', 
          #'Control',
          'Humoral\nrecognition', 
          #'Control',
          #'RNAi',
          #'Control',
         ]
data =   [alphaBar(immune), 
          alphaBar(immune_control),
          0,
          alphaBar(cellular),
          #alphaBar(cellular_control),
          alphaBar(humoral),
          #alphaBar(humoral_control),
          #alphaBar(RNAi_class),
          #alphaBar(RNAi_class_control),
          0,
          alphaBar(signalling),
          #alphaBar(signalling_control),
          #alphaBar(amp),
          #alphaBar(amp_control),
          alphaBar(cellular_recognition),
          #alphaBar(cellular_recognition_control),
          alphaBar(humoral_recognition),
          #alphaBar(humoral_recognition_control),
          #alphaBar(RNAi_cell_hum),
          #alphaBar(RNAi_cell_hum_control),
         ]

xlocations = na.array(range(len(data)))+0.5
width = 0.7
plot1 = bar(xlocations, data, width=width, color=['b','r','b','b','b','b','b','b','b'])
subplots_adjust(top=0.9, bottom=0.25)
yticks(na.array([-0.4, -0.2, 0, 0.2, 0.4]))
xticks(xlocations+ width/2, labels, rotation='vertical')
xlim(0, xlocations[-1]+width*2)
title("AlphaBar estimates")
axhline(data[0], ls='--', zorder=0)
axhline(data[1], ls='--', color='r', zorder=0)
gca().get_xaxis().tick_bottom()
gca().get_yaxis().tick_left()

plt.savefig("alphaBar_erato.png")
plt.close()

data =   [aBar(immune), 
          aBar(immune_control),
          0,
          aBar(cellular),
          #aBar(cellular_control),
          aBar(humoral),
          #aBar(humoral_control),
          #aBar(RNAi_class),
          #aBar(RNAi_class_control),
          0,
          aBar(signalling),
          #aBar(signalling_control),
          #aBar(amp),
          #aBar(amp_control),
          aBar(cellular_recognition),
          #aBar(cellular_recognition_control),
          aBar(humoral_recognition),
          #aBar(humoral_recognition_control),
          #aBar(RNAi_cell_hum),
          #aBar(RNAi_cell_hum_control),
         ]

xlocations = na.array(range(len(data)))+0.5
width = 0.7
plot1 = bar(xlocations, data, width=width, color=['b','r','b','b','b','b','b','b','b'])
subplots_adjust(top=0.9, bottom=0.25)
yticks(na.array([-0.002, -0.001, 0, 0.001, 0.002, 0.003]))
xticks(xlocations+ width/2, labels, rotation='vertical')
xlim(0, xlocations[-1]+width*2)
title("aBar estimates")
axhline(data[0], ls='--', zorder=0)
axhline(data[1], ls='--', color='r', zorder=0)
gca().get_xaxis().tick_bottom()
gca().get_yaxis().tick_left()

plt.savefig("aBar_erato.png")
plt.close()
