import numpy as np
import matplotlib.pyplot as plt
import math
import os
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as mtick

plt.figure(figsize=(6,6))
plt.rcParams.update({'font.size': 16})

file = open("plots/benchmarks.txt", "r")
lines = file.readlines()
ns = np.zeros(len(lines))
times_dsdv = np.zeros(len(lines))
times_eqmotion = np.zeros(len(lines))

for i, line in enumerate(lines):
    ns[i] = int(line.strip().split(" ")[0])
    times_dsdv[i] = float(line.strip().split(" ")[2])
    times_eqmotion[i] = float(line.strip().split(" ")[3])

# plot system momentum and energy
ns2 = ns**2
ns_log = np.log(ns)
dsdv_log = np.log(times_dsdv)
plt.clf()
plt.ticklabel_format(style='plain', axis='x')
plt.scatter(ns, times_eqmotion, s=0.5, label="Eq. of motion iterations")
plt.plot(ns, times_eqmotion)
plt.scatter(ns, times_dsdv, s=0.5, label="Entropy gradient")
plt.plot(ns, times_dsdv)
legend = plt.legend()
legend.legendHandles[0]._sizes = [30]
legend.legendHandles[1]._sizes = [30]
plt.xlabel("Markers per dimension and per specie")
plt.ylabel("Computational time per timestep [ms]")
plt.yscale("log", base=2)
plt.xscale("log", base=2)
plt.gca().xaxis.set_major_formatter(ScalarFormatter())
plt.gca().yaxis.set_major_formatter(ScalarFormatter())
plt.gcf().subplots_adjust(left=0.2)
plt.savefig("bench.eps")
plt.savefig("bench.png")
