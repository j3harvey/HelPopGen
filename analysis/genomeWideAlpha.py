import csv
import subprocess

for th in ["0.00",
           "0.02",
           "0.04",
           "0.06",
           "0.08",
           "0.10",
           "0.12",
           "0.14",
           "0.16",
           "0.18",
           "0.20",
           "0.22",
           "0.24",
           "0.26",
           "0.28",
           "0.30",
           "0.32",
           "0.34",
           "0.36",
           "0.38",
           "0.40",]:
    ###
    #
    # Load HelPopGen output
    #
    
    with open("../results/paper" + th + ".csv", 'r') as f:
        hpg_names = f.readline().rstrip().split(',')
        hpg = [x for x in csv.DictReader(f, hpg_names)]
    
    ###
    #
    # Create MKtest input file for classes or loci
    #
    
    with open("mktest_genomeWide" + th + ".csv", 'w') as f:
        for x in hpg:
            if float(x['Ssites']) > 3 and float(x['NSsites']) > 3:
                f.write(','.join([x['D_N'],
                                  str(float(x['NSsites'])/3),
                                  x['in_P_N'],
                                  str(float(x['NSsites'])/3),
                                  x['D_S'],
                                  str(float(x['Ssites'])/3),
                                  x['in_P_S'],
                                  str(float(x['Ssites'])/3),
                                  x['ingroup_sequences'],
                                 ]) + '\n')
    
    ###
    #
    # Run MKtest
    #
    
    subprocess.call(['../../MKtest-2.0/MKtest', '-f2', 'mktest_genomeWide' + th + '.csv', '-o', 'genomeWideAlpha' + th + '.mkout'])



