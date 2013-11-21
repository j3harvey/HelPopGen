#! /usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt

data = np.array([0.207,
                 0.164,
                 0.562,
                 0.0207,
                 0.585,
                 0.349,
                ])

Lower = np.array([0.185,
                  0.157,
                  0.409,
                 -0.22,
                  0.559,
                  0.034,
        ])

Upper = np.array([0.247,
         0.524,
         0.597,
         0.136,
         0.620,
         0.495,
        ])

Err = (Upper - data,
       data - Lower)
 

labels = ["control",
          "humoral recognition",
          "cellular recognition",
          "signalling",
          "AMPs",
          "RNAi",
         ]

xlocations = np.array(range(len(data)))+0.5
width = 0.7
plot1 = plt.bar(xlocations, data, width=width, color=['b', 'r', 'r', 'r', 'r', 'r'], yerr=Err, ecolor='k')
subplots_adjust(top=0.95, bottom=0.25)
#yticks(na.array([-0.4, -0.2, 0, 0.2, 0.4]))
xticks(xlocations+ width/2, labels, rotation='vertical')
xlim(0, xlocations[-1]+width*2)
title("Alpha estimates by gene class")
ylabel("alpha")
gca().get_xaxis().tick_bottom()
gca().get_yaxis().tick_left()


