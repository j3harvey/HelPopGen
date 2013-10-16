with open('classes_bootstrap.mkout', 'r') as f:
    lines = f.readlines()
    
lines = [x.rstrip().split('\t') for x in lines]
alpha = [x[-6:] for x in lines]
classes = zip(*sorted(zip(alpha[0], zip(*[map(float, x) for x in alpha[1:]]))))[1]
classes = [sorted(c) for c in classes]

with open('classes_alpha.mkout', 'r') as f:
    lines = f.readlines()
    par_names = lines[0].rstrip().split('\t')
    par_vals = lines[1].rstrip().split('\t')
    alpha = [(x, float(y)) for (x, y) in zip(par_names, par_vals) if x.startswith('alpha')]

upper = array([c[4] for c in classes])
middle = array(zip(*sorted(alpha))[1])
lower = array([c[-5] for c in classes])

labels = ["control",
          "AMPs",
          "RNAi",
          "cellular recognition",
          "humoral recognition",
          "signalling",
        ]

figure()
xlocations = np.array(range(len(middle)))+0.5
width = 0.7
plot1 = plt.bar(xlocations, middle, width=width, color=['b', 'r', 'r', 'r', 'r', 'r'], 
                yerr=(middle - lower, upper - middle), ecolor='k')
subplots_adjust(top=0.95, bottom=0.3)
#yticks(na.array([-0.4, -0.2, 0, 0.2, 0.4]))
xticks(xlocations+ width/2, labels, rotation='vertical', fontsize='16')
xlim(0, xlocations[-1]+width*2)
title(r'$\alpha$ estimates by gene class', fontsize='16')
ylabel(r'$\alpha$', rotation='horizontal', fontsize='20')
gca().get_xaxis().tick_bottom()
gca().get_yaxis().tick_left()


### a

with open('classes_A_bootstrap.mkout', 'r') as f:
    lines = f.readlines()
    
lines = [x.rstrip().split('\t') for x in lines]
alpha = [x[-6:] for x in lines]
classes_a = zip(*sorted(zip(alpha[0], zip(*[map(float, x) for x in alpha[1:]]))))[1]
classes_a = [sorted(c) for c in classes_a]

with open('classes_A.mkout', 'r') as f:
    lines = f.readlines()
    par_names = lines[0].rstrip().split('\t')
    par_vals = lines[1].rstrip().split('\t')
    alpha = [(x, float(y)) for (x, y) in zip(par_names, par_vals) if x.startswith('alpha')]

lower_a = array([c[4] for c in classes_a])
middle_a = array(zip(*sorted(alpha))[1])
upper_a = array([c[-5] for c in classes_a])

labels = ["control",
          "AMPs",
          "RNAi",
          "cellular recognition",
          "humoral recognition",
          "signalling",
        ]

figure()
xlocations = np.array(range(len(middle)))+0.5
width = 0.7
plot1 = plt.bar(xlocations, middle_a, width=width, color=['b', 'r', 'r', 'r', 'r', 'r'],
                yerr=(middle_a - lower_a, upper_a - middle_a), ecolor='k')
subplots_adjust(top=0.95, bottom=0.3)
#yticks(na.array([-0.4, -0.2, 0, 0.2, 0.4]))
xticks(xlocations+ width/2, labels, rotation='vertical', fontsize='16')
xlim(0, xlocations[-1]+width*2)
title(r'$a$ estimates by gene class', fontsize='16')
ylabel(r'a', rotation='horizontal', fontsize='20')
gca().get_xaxis().tick_bottom()
gca().get_yaxis().tick_left()


