import matplotlib.pyplot as plt
from numpy import linspace

alpha = []

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
    with open("genomeWideAlpha" + th + ".mkout", 'r') as f:
        alpha.append(float(f.readlines()[1].split('\t')[-1]))

print alpha

plt.plot( linspace(0, 0.4, 21), alpha )
xlabel("cut-off frequency")
ylabel("alpha")
plt.savefig("alpha.pdf")

