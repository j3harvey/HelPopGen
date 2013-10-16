import csv
from controlGenes import getControlGenes, mapped
import subprocess

###
#
# Load orthology relationships from BioMart
#

with open("../data/dmelHmelOrthologues_allSources.csv", 'r') as f:
    ortho_names = f.readline().rstrip().split(',')
    orthologues = [x for x in csv.DictReader(f, ortho_names)]

###
#
# Load Obbard data
#

with open("../data/s027_all.csv", 'r') as f:
    obb_names = f.readline().rstrip().split(',')
    obb = [x for x in csv.DictReader(f, obb_names)]

###
#
# Load HelPopGen output
#

with open("../results/results_thres0.2.csv", 'r') as f:
    hpg_names = f.readline().rstrip().split(',')
    hpg = [x for x in csv.DictReader(f, hpg_names)]

###
#
# Load completer list of immunity genes
#

with open("../data/immunity_genes.csv", 'r') as f:
    imm = [x.rstrip().split(',') for x in f.readlines()]

###
#
# Produce lists of data for mapped genes with orthologues
#


dmel_ids = [x['FBgn'] for x in obb]
hmel_ids = [x['locus'] for x in hpg]
mapped_hmel_ids = mapped(hmel_ids)
orthologues_with_data = [x for x in orthologues
                         if x['Dmel_ID'] in dmel_ids
                         and x['Hmel_ID'] in mapped_hmel_ids]
dmel_ids_with_orthologues = list()
for o in orthologues_with_data:
    for l in obb:
        if o['Dmel_ID'] == l['FBgn']:
            dmel_ids_with_orthologues.append(dict(l)) # use dict explicitly to force deep copy
hmel_ids_with_orthologues = list()
for o in orthologues_with_data:
    for l in hpg:
        if o['Hmel_ID'] == l['locus']:
            hmel_ids_with_orthologues.append(dict(l)) # use dict explicitly to force deep copy
#control_gene_names = getControlGenes([x['locus'] for x in hmel_ids_with_orthologues], 20)
control_gene_names = getControlGenes([y[0] for y in imm], 10)
control_genes = [x for x in hpg if x['locus'] in control_gene_names]
immunity_genes = [x for x in hpg if x['locus'] in [y[0] for y in imm]]

###
#
# Add class information to H. mel genes
#

for (g, d) in zip(hmel_ids_with_orthologues, dmel_ids_with_orthologues):
    g['gene_name'] = d['Locus']
    g['Class'] = d['Class']
    g['Cell_Hum'] = d['Cell_Hum']
for g in control_genes:
    g['gene_name'] = ''
    g['Class'] = 'control'
    g['Cell_Hum'] = 'control'
for g in immunity_genes:
    for h in imm:
        if h[0] == g['locus']:
            g['Class'] = h[2]
            g['Cell_Hum'] = h[3]
            g['gene_name'] = h[1]

###
#
# Create MKtest input file for classes or loci
#

classes = {'control' : '0',
           'immune' : '1',
          }

with open("mktest_immuneVScontrol.csv", 'w') as f:
    for x in immunity_genes:
            f.write(','.join([x['D_N'],
                                  x['Ln'],
                                  x['in_P_N'],
                                  x['Ln'],
                                  x['D_S'],
                                  x['Ls'],
                                  x['in_P_S'],
                                  x['Ls'],
                                  x['ingroup_sequences'],
                                  '0',
                                  '1',
                                ]) + '\n')
    for x in control_genes:
            f.write(','.join([x['D_N'],
                                  x['Ln'],
                                  x['in_P_N'],
                                  x['Ln'],
                                  x['D_S'],
                                  x['Ls'],
                                  x['in_P_S'],
                                  x['Ls'],
                                  x['ingroup_sequences'],
                                  '0',
                                  '0',
                                 ]) + '\n')

###
#
# Create MKtest input file for classes or loci
#

classes = {'control' : '0',
           'AMP' : '1',
           'RNAi' : '2',
           'cellular_recognition' : '3',
           'humoral_recognition' : '4',
           'signaling' : '5',
          }

with open("mktest_classes.csv", 'w') as f:
    for x in immunity_genes:
        if x['Cell_Hum'] in classes.keys():
            f.write(','.join([x['D_N'],
                                  x['Ln'],
                                  x['in_P_N'],
                                  x['Ln'],
                                  x['D_S'],
                                  x['Ls'],
                                  x['in_P_S'],
                                  x['Ls'],
                                  x['ingroup_sequences'],
                                  '0',
                                  classes[x['Cell_Hum']],
                                ]) + '\n')
    for x in control_genes:
        if x['Cell_Hum'] in classes.keys():
            f.write(','.join([x['D_N'],
                                  x['Ln'],
                                  x['in_P_N'],
                                  x['Ln'],
                                  x['D_S'],
                                  x['Ls'],
                                  x['in_P_S'],
                                  x['Ls'],
                                  x['ingroup_sequences'],
                                  '0',
                                  classes[x['Cell_Hum']],
                                 ]) + '\n')

###
#
# Create MKtest input file for comparison with drosophila genes
#

classes = {'control' : '0',
           'AMP' : '1',
           'RNAi' : '2',
           'cellular_recognition' : '3',
           'humoral_recognition' : '4',
           'signaling' : '5',
          }

with open("mktest_vsdmel.csv", 'w') as f:
    for x in hmel_ids_with_orthologues:
        if x['Cell_Hum'] in classes.keys():
            f.write(','.join([x['D_N'],
                                  x['Ln'],
                                  x['in_P_N'],
                                  x['Ln'],
                                  x['D_S'],
                                  x['Ls'],
                                  x['in_P_S'],
                                  x['Ls'],
                                  x['ingroup_sequences'],
                                  '0',
                                  classes[x['Cell_Hum']],
                                ]) + '\n')
    for x in control_genes:
        if x['Cell_Hum'] in classes.keys():
            f.write(','.join([x['D_N'],
                                  x['Ln'],
                                  x['in_P_N'],
                                  x['Ln'],
                                  x['D_S'],
                                  x['Ls'],
                                  x['in_P_S'],
                                  x['Ls'],
                                  x['ingroup_sequences'],
                                  '0',
                                  classes[x['Cell_Hum']],
                                 ]) + '\n')

###
#
# Run MKtest
#

#subprocess.call(['../../MKtest-2.0/MKtest', '-a5', '-A', 'mktest_in.csv', '-o', 'loci.mkout'])

###
#
# Scatter plot of Ln vs a for immune and control genes
#

with open('loci.mkout', 'r') as f:
    lines = f.readlines()
    par_names = lines[0].rstrip().split('\t')
    par_vals = [float(x) for x in lines[1].rstrip().split('\t')]

def bootstrap(V, a1, a2, n=1000):
    '''Return the a1 and a2 points of the bootstrapping distribution of the mean of V
    '''
    l = len(V)
    samples = list()
    for i in range(n):
        samples.append(mean([random.choice(V) for j in range(l)]))
    samples.sort()
    i1 = int(round(n*a1)) - 1
    i2 = int(round(n*a2)) - 1
    return (samples[i1], samples[i2])

def bstd(V, n=1000):
    '''Return the bootstrap standard error in the mean of V
    '''
    l = len(V)
    samples = list()
    for i in range(n):
        samples.append(mean([random.choice(V) for j in range(l)]))
    return std(samples)

a = [x for (p,x) in zip(par_names, par_vals) if p.startswith('alpha')] # But these are 'a' estimates not 'alpha'
Lni = [float(x['Ln']) for x in immunity_genes]
Lnc = [float(x['Ln']) for x in control_genes]
ai = a[:len(Lni)]
ac = a[len(Lni):]
names = [x['gene_name'] for x in immunity_genes]

figure()
scatter(Lnc, ac)
scatter(Lni, ai, color='r')
xlim(xmin=-100)
ylabel(r'adaptive substitutions per non-synonymous site, $a$', fontsize=14)
xlabel(r'number of non-synonymous sites, $L_N$', fontsize=14)
legend(('control', 'immune'), loc='lower right', fontsize=14)

figure()
hist(ac, linspace(-0.2, 0.101, 2*len(ac)), normed=True, cumulative=True, histtype='step', color='b')
hist(ai, linspace(-0.2, 0.101, 2*len(ac)), normed=True, cumulative=True, histtype='step', color='r')
xlim(-0.2, 0.1)
ylim(0, 1.02)
xlabel(r'adaptive substitutions per non-synonymous site, $a$', fontsize=14)
ylabel(r'cumulative frequency', fontsize=16)
legend(('control', 'immune'), loc='upper left', fontsize=14)

###
#
# Hmel vs Dmel
#

hmel_a = list()
loci = [x['locus'] for x in immunity_genes]
dmel_a = list()
for (h, d) in zip(hmel_ids_with_orthologues, dmel_ids_with_orthologues):
    for (l, a) in zip(loci, ai):
        if l == h['locus']:
            hmel_a.append(a)
    dmel_a.append(d['a'])



hmel_d = dict()
dmel_d = dict()
for (h, d) in zip(hmel_ids_with_orthologues, dmel_ids_with_orthologues):
    for (l, a) in zip(loci, ai):
        if l == h['locus']:
            K = "nim" if d['Locus'].startswith("nim") else d['Locus']
            if K in dmel_d.keys():
                hmel_d[K].append(a) 
                dmel_d[K].append(float(d['a']))
            else:
                hmel_d[K] = [a,]
                dmel_d[K] = [float(d['a']),]
for (k, v) in dmel_d.items():
    dmel_d[k] = mean(v)
for (k, v) in hmel_d.items():
    hmel_d[k] = mean(v)
hmel_a, dmel_a = zip(*[(hmel_d[k], dmel_d[k]) for k in hmel_d.keys()])

figure()
scatter(hmel_a, dmel_a)
xlabel(r'$a$ estimate from $Heliconius$', fontsize=14)
ylabel(r'$a$ estimate from $Drosophila$', fontsize=16)

for (h, d) in zip(hmel_a, dmel_a):
    if h < -0.01 or d > 0.021 or d < -0.01:
        #print h, d, {v : ('IKKg' if k.startswith('ird5') else k) for (k, v) in hmel_d.items()}[h]
        #scatter([h,], [d,], color='r')
        text(h + 0.001, d - 0.001, {v: (r'IKKg' if k.startswith('ird5') else k) for (k, v) in hmel_d.items()}[h])

with open('loci.mkout', 'r') as f:
    lines = f.readlines()
    par_names = lines[0].rstrip().split('\t')
    par_vals = [float(x) for x in lines[1].rstrip().split('\t')]

a = [x for (p,x) in zip(par_names, par_vals) if p.startswith('alpha')]
Lni = [float(x['Ln']) for x in immunity_genes]
Lnc = [float(x['Ln']) for x in control_genes]
ai = alpha[:len(Lni)]
ac = alpha[len(Lni):]


###
#
# Comparison of alpha values for different genes
#

with open('loci_alpha.mkout', 'r') as f:
    lines = f.readlines()
    par_names = lines[0].rstrip().split('\t')
    par_vals = [float(x) for x in lines[1].rstrip().split('\t')]

alpha = [x for (p,x) in zip(par_names, par_vals) if p.startswith('alpha')]
Lni = [float(x['Ln']) for x in immunity_genes]
Lnc = [float(x['Ln']) for x in control_genes]
alphai = alpha[:len(Lni)]
alphac = alpha[len(Lni):]

names = [x['gene_name'] for x in immunity_genes]
for (i, g) in enumerate(immunity_genes):
    g['a'] = ('%.6f' % ai[i])
for (i, g) in enumerate(control_genes):
    g['a'] = ('%.6f' % ac[i])
for (i, g) in enumerate(immunity_genes):
    g['alpha'] = ('%.3f' % alphai[i])
for (i, g) in enumerate(control_genes):
    g['alpha'] = ('%.3f' % alphac[i])

###
#
# Write lovely tables of genes used in the analysis
#

with open("../results/results_thres0.2.csv", 'r') as f:
    colnames = f.readline().rstrip().split(',') + ['a', 'alpha', 'gene_name', 'Cell_Hum', 'Class']

with open("immune_and_control_a_and_alpha.csv", 'w') as f:
    F = DictWriter(f, colnames)
    F.writeheader()
    for g in immunity_genes:
        F.writerow(g)
    for g in control_genes:
        F.writerow(g)

