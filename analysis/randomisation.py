with open('randomisation.mkout', 'r') as f:
    LnL = [float(l.rstrip().split('\t')[0]) for l in f.readlines()][1:]

LnL0 = LnL[0]
LnL = LnL[1:]
print sum([x > LnL0 for x in LnL])/float(len(LnL))

