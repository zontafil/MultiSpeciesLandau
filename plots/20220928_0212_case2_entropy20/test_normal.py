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
# from pingouin import multivariate_normality
# import pandas as pd
from scipy.optimize import curve_fit

# Return a gaussian distribution at an angle alpha from the x-axis
# from astroML for use with curve_fit
def mult_gaussFun_Fit2(z,A,x0,y0,varx,vary):
    X = z[0]
    Y = z[1]
    # X,Y = np.meshgrid(x,y)
    # f = A*np.exp(-a*((X-x0)**2/(varx)+(Y-y0)**2/(vary)-(2*rho/(np.sqrt(varx*vary)))*(X-x0)*(Y-y0)))
    f = A*np.exp(-((X-x0)**2/(varx)+(Y-y0)**2/(vary)))
    # return Z.ravel()
    return f
def mult_gaussFun_Fit(z,A,x0,y0,varx,vary,rho,alpha):
    X = z[0]
    Y = z[1]
    # X,Y = np.meshgrid(x,y)
    assert rho != 1
    a = 1/(2*(1-rho**2))
    # f = A*np.exp(-a*((X-x0)**2/(varx)+(Y-y0)**2/(vary)-(2*rho/(np.sqrt(varx*vary)))*(X-x0)*(Y-y0)))
    f = A*np.exp(-a*((X-x0)**2/(varx)+(Y-y0)**2/(vary)-(2*rho/(np.sqrt(varx*vary)))*(X-x0)*(Y-y0)))
    # return Z.ravel()
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

if printSecondDistr:
    f2file = open("f0out.txt", "r")
    f2lines = f2file.readlines()
    times2 = np.zeros(len(f2lines))
    errs2 = np.zeros(len(f2lines))
    for i in range(len(f2lines)):
        line = f2lines[i]
        line_s = line.strip().split(" ")
        times2[i] = float(line_s[0])
        errs2[i] = float(line_s[1])

    plt.clf()
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.scatter(times2, errs2, s=0.5)
    plt.plot(times2, errs2)
    plt.xlabel("Time [s]")
    plt.grid()
    plt.ylabel("f0 mean err")
    plt.gcf().subplots_adjust(left=0.25)
    plt.savefig("err_maxwell.png")


if printFirstDistr:
    for file_i, filename in enumerate(files):
        file = open(filename, "r")
        matches = re.findall(r"dist_(\d+)\.txt", filename)
        t = int(matches[0])
        times[file_i] = t * dt

        # if file_i > 100:
        #     continue

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
                # print("dio {:.2E}".format(f[specie][id]))

        # record distribution
        if file_i == 0:
            print("nspecies {}".format(nspecies))
            print("n {}".format(n))
            print("n Markers/Dim {}".format(nMarkersPerDim))
        # else:
        #     break
        if file_i % 20 == 0:
            print(filename)

        z = np.array((vx[0], vy[0]))
        
        p0 = (f[0].max(), z[0].mean(), z[1].mean(), np.cov(vx[0]), np.cov(vy[1]), 0.5, np.pi/4)
        p0 = (0.6366198/4, 0, 0, 8, 8)
        # p0 = (0.25, 0, 8)

        coeff, var_matrix = curve_fit(mult_gaussFun_Fit2,z,f[0], p0=p0)
        # coeff, var_matrix = curve_fit(mult_gaussFun_Fit,z,f[0], p0=p0)
        # coeff, var_matrix = curve_fit(gaussFun_Fit,vx[0],f[0], p0=p0)

        A = coeff[0]
        x0 = coeff[1]
        y0 = coeff[2]
        varx = coeff[3]
        vary = coeff[4]
        # A = 0.6366198/4
        # x0 = 0
        # v0 = 0
        # varx = 8
        # vary = 8
        f_out = A*np.exp(-((z[0]-x0)**2/(varx)+(z[1]-y0)**2/(vary)))
        print("v {:.2E} {:.2E} {:.2E} {:.2E}".format(f_out[0], f[0][0], vx[0][0], vy[0][0]))

        # A = coeff[0]
        # x0 = coeff[1]
        # y0 = coeff[2]
        # varx = coeff[3]
        # vary = coeff[4]
        # rho = coeff[5]
        # alpha = coeff[6]
        # a = 1/(2*(1-rho**2))
        # f_out = A*np.exp(-a*((z[0]-x0)**2/(varx)+(z[1]-y0)**2/(vary)-(2*rho/(np.sqrt(varx*vary)))*(z[0]-x0)*(z[1]-y0)))

        # A = coeff[0]
        # x0 = coeff[1]
        # varx = coeff[2]
        # f_out = A*np.exp(-((vx[0]-x0)**2/(varx)))
        # f_out = 0.25*np.exp(-((vx[0])**2/(8)))

        df = f[0] - f_out
        # print(f_out)
        errs[file_i] = np.mean(np.abs(df))
        # print("ferr {}".format(err))


    # plt.clf()
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.scatter(times, errs, s=0.5)
    plt.plot(times, errs)
    plt.xlabel("Time [s]")
    # plt.grid()
    plt.ylabel("f0 mean err")
    plt.gcf().subplots_adjust(left=0.25)
    plt.savefig("err_maxwell.png")

# fout = open("f0out.txt", "w")
# print(len(times))
# for i in range(len(times)):
#     fout.write("{:.2E} {:.2E}\n".format(times[i], errs[i]))

    # print(coeff)
    # print(var_matrix)
    # # coeff, var_matrix = curve_fit(myfunc,vx[0], f[0])
    

