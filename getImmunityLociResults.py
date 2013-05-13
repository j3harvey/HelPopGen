#! /usr/bin/env python2.7

with open("data/ImmunityGenes.csv", 'r') as f:
    lines = f.readlines()

ids = [l.split(',')[0].rstrip() for l in lines]
names = [l.split(',')[1].rstrip() for l in lines]
dictNames = {i:n for (i, n) in zip(ids, names)}

print ids
print names

count = 0
with open("results/third_run_00.csv", 'r') as f:
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

