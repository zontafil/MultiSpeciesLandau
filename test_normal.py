from cmath import tau
import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import math
from plot_utils import createDir
from plot_analytic_times import PslowDown, Tisothermal, Tthermal
from plot_charts import MultiSpeciesCharts
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as mtick
from scipy.optimize import curve_fit

# Return a gaussian distribution at an angle alpha from the x-axis
# from astroML for use with curve_fit
def mult_gaussFun_Fit2(z,A,x0,y0,varx,vary):
    X = z[0]
    Y = z[1]
    f = A*np.exp(-((X-x0)**2/(varx)+(Y-y0)**2/(vary)))
    return f
def mult_gaussFun_Fit(z,A,x0,y0,varx,vary,rho,alpha):
    X = z[0]
    Y = z[1]
    assert rho != 1
    a = 1/(2*(1-rho**2))
    f = A*np.exp(-a*((X-x0)**2/(varx)+(Y-y0)**2/(vary)-(2*rho/(np.sqrt(varx*vary)))*(X-x0)*(Y-y0)))
    return f
def gaussFun_Fit(X,A,x0,varx):
    f = A*np.exp(-((X-x0)**2/(varx)))
    return f

plt.figure(figsize=(6,6))
plt.rcParams.update({'font.size': 16})

enableDistributionPlot = False
enableLinePlots = False
enableAnalyticFunctions = True

files = glob.glob("./out/data/dist_*.txt")
files.sort()

fmax = None
nspecies = 0
file = open(files[0], "r")
lines = file.readlines()
nspecies = int(lines[0].strip().split(" ")[0])
nMarkersPerDim = int(math.sqrt(float(lines[0].strip().split(" ")[1])))
n = int(math.sqrt(float(lines[0].strip().split(" ")[2])))
dt = float(lines[0].strip().split(" ")[3])
vx = np.zeros([nspecies, n*n])
vy = np.zeros([nspecies, n*n])
f = np.zeros([nspecies, n*n])
ntimes = len(files)
times = np.zeros(ntimes)
errs = np.zeros(ntimes)

printFirstDistr = True
printSecondDistr = True

for file_i, filename in enumerate(files):
    file = open(filename, "r")
    matches = re.findall(r"dist_(\d+)\.txt", filename)
    t = int(matches[0])
    times[file_i] = t * dt

    print("Analizing file {} of {}".format(file_i, len(files)))

    lines = file.readlines()
    if len(lines) > 2:
        for i in range(1, len(lines)):
            line = lines[i]
            line_s = line.strip().split(" ")
            specie = int(line_s[0])
            id = int(line_s[1])
            vx[specie][id] = float(line_s[2])
            vy[specie][id] = float(line_s[3])
            f[specie][id] = float(line_s[4])

    # record distribution
    if file_i == 0:
        print("nspecies {}".format(nspecies))
        print("n {}".format(n))
        print("n Markers/Dim {}".format(nMarkersPerDim))
    if file_i % 20 == 0:
        print(filename)

    z = np.array((vx[0], vy[0]))
    
    p0 = (f[0].max(), z[0].mean(), z[1].mean(), np.cov(vx[0]), np.cov(vy[1]), 0.5, np.pi/4)
    p0 = (0.6366198/4, 0, 0, 8, 8)

    coeff, var_matrix = curve_fit(mult_gaussFun_Fit2,z,f[0], p0=p0)

    A = coeff[0]
    x0 = coeff[1]
    y0 = coeff[2]
    varx = coeff[3]
    vary = coeff[4]
    f_out = A*np.exp(-((z[0]-x0)**2/(varx)+(z[1]-y0)**2/(vary)))
    # print("v {:.2E} {:.2E} {:.2E} {:.2E}".format(f_out[0], f[0][0], vx[0][0], vy[0][0]))

    df = f[0] - f_out
    errs[file_i] = np.mean(np.abs(df))

    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.scatter(times, errs, s=0.5)
    plt.plot(times, errs)
    plt.xlabel("Time [s]")
    plt.ylabel("|f - f0| mean error")
    plt.gcf().subplots_adjust(left=0.25)
    plt.savefig("err_maxwell.png")
    plt.savefig("err_maxwell.eps")
