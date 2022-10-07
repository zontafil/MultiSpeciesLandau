from cmath import tau
import sys
import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import math
# import cv2
from plot_utils import createDir
from plot_analytic_times import PslowDown, Tisothermal, Tthermal
from plot_charts import MultiSpeciesCharts
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as mtick

plt.figure(figsize=(6,6))
plt.rcParams.update({'font.size': 16})

enableDistributionPlot = True
enableLinePlots = False
enableAnalyticFunctions = True

files = glob.glob("./out/data/step_C*.txt")
files.sort()
files = files[:2]

fmax = None
nspecies = 0
file = open(files[0], "r")
lines = file.readlines()
nspecies = int(lines[0].strip().split(" ")[0])
fmax = np.zeros(nspecies)

# read fixed parameters
file = open(files[0], "r")
lines = file.readlines()
nMarkersPerDim = int(math.sqrt(float(lines[0].strip().split(" ")[1])))
nPerDim = int(math.sqrt(float(lines[0].strip().split(" ")[2])))
nspecies = int(lines[0].strip().split(" ")[0])
nmesh = [0] * nspecies
dt = float(lines[0].strip().split(" ")[3])
for s in range(nspecies):
    line = lines[s+2].strip().split(" ")
    if len(line) > 9:
        nmesh[s] = int(math.sqrt(float(line[9])))
        print("nmesh {} {}".format(s, nmesh[s]))
    else:
        nmesh[s] = nPerDim
n = np.max(nmesh)
if len(lines[0].strip().split(" ")) > 4:
    clog00 = float(lines[0].strip().split(" ")[4])
    clog01 = float(lines[0].strip().split(" ")[5])
    clog11 = float(lines[0].strip().split(" ")[6])
else:
    # use some dummy default values for old sim outputs
    clog00 = 14.6
    clog01 = 6.6
    clog11 = 14.1


if enableDistributionPlot == True:
    for file_i, filename in enumerate(files):
        file = open(filename, "r")
        matches = re.findall(r"step_C_(\d+)\.txt", filename)
        t = int(matches[0])

        print("Analizing file {} of {}".format(file_i, len(files)))

        lines = file.readlines()
        if len(lines) > 2+nspecies:
            for i in range(2+nspecies, len(lines)):
                line = lines[i]
                line_s = line.strip().split(" ")
                if (len(line_s)<5):
                    continue
                specie = int(line_s[0])
                index = int(line_s[1])
                fmax[specie] = max(fmax[specie], float(line_s[4]))

        # record distribution
        if file_i == 0:
            print("nspecies {}".format(nspecies))
            print("n {}".format(n))
            print("n Markers/Dim {}".format(nMarkersPerDim))
        if file_i % 20 == 0:
            print(filename)

    vx = np.zeros([nspecies, n*n])
    vy = np.zeros([nspecies, n*n])
    specie = np.zeros([nspecies, n*n])
    f = np.zeros([nspecies, n*n])
    for file_i, filename in enumerate(files):
        matches = re.findall(r"step_C_(\d+)\.txt", filename)
        t = int(matches[0])

        print("Plotting distribution for file {} of {} (time {})".format(file_i, len(files), t))
        file = open(filename, "r")
        lines = file.readlines()
        if len(lines) > 2+nspecies:
            for i in range(2+nspecies, len(lines)):
                line = lines[i]
                line_s = line.strip().split(" ")
                if len(line_s) >= 5:
                    specie = int(line_s[0])
                    index = int(line_s[1])
                    vx[specie][index] = float(line_s[2])
                    vy[specie][index] = float(line_s[3])
                    f[specie][index] = float(line_s[4])

            #plot distribution to file
            plt.cla()

            # plot distribution for current time step
            colourMaps = ["viridis", "inferno", "viridis", "summer"]
            alphas = [0.6, 0.6]
            resolutions = [20, 20]
            vmaxes = [fmax[0]/2, fmax[1]]
            if nspecies > 1:
                vmaxes.append(fmax[1])
            for i in range(nspecies):
                nmeshs = nmesh[i]
                vxs = vx[i][:nmeshs*nmeshs]
                vys = vy[i][:nmeshs*nmeshs]
                fs = f[i][:nmeshs*nmeshs]
                plt.contourf(vxs[:].reshape(nmeshs,nmeshs), vys[:].reshape(nmeshs,nmeshs), fs[:].reshape(nmeshs,nmeshs), resolutions[i], vmax=vmaxes[i], alpha=alphas[i], cmap=colourMaps[i])
            imgname = 'out/data/plot_C_' + str(t+1).zfill(6) + '.png'
            plt.axis('equal')
            plt.ylabel("vy [m/s]")
            plt.xlabel("vx [m/s]")
            plt.draw()
            plt.savefig(imgname)

if enableLinePlots:
    E = np.zeros(len(files))
    Eerr = np.zeros(len(files))
    dS = np.zeros(len(files))
    P = np.zeros((len(files), 2))
    Pnorm = np.zeros(len(files))
    Perr = np.zeros((len(files), 2))
    times = np.zeros(len(files))
    thermalAnalytic = np.zeros(len(files))
    Especies = np.zeros((nspecies, len(files)))
    Tspecies = np.zeros((nspecies, len(files)))
    TspeciesAxes = np.zeros((nspecies, 2, len(files)))
    VFlowSpecies = np.zeros((nspecies, 2, len(files)))
    VFlowSpeciesNorm = np.zeros((nspecies, len(files)))
    ns = np.zeros((nspecies, len(files)))
    ms = np.zeros((nspecies, len(files)))
    specieNames = []

    print("Number of files: {}".format(len(files)))

    for file_i, filename in enumerate(files):
        print(str(file_i))
        matches = re.findall(r"step_C_(\d+)\.txt", filename)
        t = int(matches[0])
        file = open(filename, "r")
        lines = file.readlines()
        idx = np.zeros([nspecies, n*n])

        # store Energy and momentum
        line = lines[1].strip().split(" ")
        E[file_i] = float(line[0])
        Eerr[file_i] = float(line[1])
        P[file_i,0] = float(line[2])
        P[file_i,1] = float(line[3])
        Pnorm[file_i] = np.sqrt(float(line[3])**2 + float(line[2])**2)
        Perr[file_i,0] = float(line[4])
        Perr[file_i,1] = float(line[5])
        thermalAnalytic[file_i] = float(line[6])
        if len(line) > 7:
            dS[file_i] = float(line[7])
        times[file_i] = t
        for s in range(nspecies):
            line = lines[2+s].strip().split(" ")
            Especies[s, file_i] = float(line[0])
            VFlowSpecies[s, 0, file_i] = float(line[1])
            VFlowSpecies[s, 1, file_i] = float(line[2])
            VFlowSpeciesNorm[s, file_i] = np.sqrt(VFlowSpecies[s,0,file_i]**2 + VFlowSpecies[s,1,file_i]**2)
            Tspecies[s, file_i] = float(line[3])
            TspeciesAxes[s, 0, file_i] = float(line[4])
            TspeciesAxes[s, 1, file_i] = float(line[5])
            specieNames.append(line[6])

            if len(line) > 7:
                ns[s, file_i] = float(line[7])
                ms[s, file_i] = float(line[8])
            else:
                # use dummy data for old sim output
                ms[0, file_i] = 9.1093837*10**(-31)
                ms[1, file_i] = 3.3435837724*10**(-27)
                ns[0, file_i] = 10**20
                ns[1, file_i] = 10**20

        file.close()

    print(dS)

    # find final temperature for case 1 in paper.
    # assumes that it's a relaxation example with equal masses and thermalized species
    Te0 = Tspecies[0, 0]
    me = 9.1093837*10**(-31)
    k = 1.602176565*10**(-19)
    ve0 = np.array([VFlowSpecies[0,0,0], VFlowSpecies[0,1,0]])
    vi0 = np.array([VFlowSpecies[1,0,0], VFlowSpecies[1,1,0]])
    Te1 = Te0 + me / (2*k) * np.linalg.norm((vi0-ve0)/2)**2
    
    # compute basic quantities
    Tax = TspeciesAxes[0, 0, 0]
    Tay = TspeciesAxes[0, 1, 0]
    Tbx = TspeciesAxes[1, 0, 0]
    Tby = TspeciesAxes[1, 1, 0]
    Tb = 0.5*(Tbx+Tby)
    Ta = 0.5*(Tax+Tay)
    ma = ms[0, 0]
    mb = ms[1, 0]
    na = ns[0, 0]
    nb = ns[1, 0]

    createDir("./plots/eps")
    createDir("./plots/png")

    # Species names override
    # specieNames = ["Species 1", "Species 2"]

    times = times * dt

    # init plot functions
    charts = MultiSpeciesCharts(plt, times, nspecies, ma, mb, Tax, Tay, Tbx, Tby, na, specieNames, clog00, clog01, clog11)

    charts.plotEnergy(Eerr)
    charts.plotEnergySpecies(Especies)
    charts.plotPnorm(Pnorm)
    charts.plotVAxes(VFlowSpecies, True, enableAnalyticFunctions)
    charts.plotPerr(Perr, Pnorm)
    charts.plotEntropy(dS)
    charts.plotTemperatureSpeciesAxes(TspeciesAxes, True, enableAnalyticFunctions)
    charts.plotTemperatureSpecies(Tspecies, False, enableAnalyticFunctions, Te1)

# VIDEO OUTPUT
# enable this and import cv2
# if enableDistributionPlot:
#     img_array = []
#     size = (0,0)
#     for filename in glob.glob("./out/data/plot_C*.png"):
#         img = cv2.imread(filename)
#         height, width, layers = img.shape
#         size = (width,height)
#         img_array.append(img)

#     print(size)
#     out = cv2.VideoWriter('./out/out.avi',cv2.VideoWriter_fourcc(*'DIVX'), 30, size)

#     for i in range(len(img_array)):
#         out.write(img_array[i])

#     out.release()
