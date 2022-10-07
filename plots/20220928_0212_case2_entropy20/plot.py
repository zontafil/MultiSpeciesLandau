from cmath import tau
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

enableDistributionPlot = False
enableLinePlots = True
enableAnalyticFunctions = True

files = glob.glob("./out/data/step_C*.txt")
files.sort()
# files = files[:100]

fmax = None
nspecies = 0
file = open(files[0], "r")
lines = file.readlines()
nspecies = int(lines[0].strip().split(" ")[0])
fmax = np.zeros(nspecies)

# read fixed parameters
file = open(files[0], "r")
lines = file.readlines()
n = int(math.sqrt(float(lines[0].strip().split(" ")[2])))
nMarkersPerDim = int(math.sqrt(float(lines[0].strip().split(" ")[1])))
nspecies = int(lines[0].strip().split(" ")[0])
dt = float(lines[0].strip().split(" ")[3])
if len(lines[0].strip().split(" ")) > 4:
    clog00 = float(lines[0].strip().split(" ")[4])
    clog01 = float(lines[0].strip().split(" ")[5])
    clog11 = float(lines[0].strip().split(" ")[6])
else:
    clog00 = 14.6
    clog01 = 6.6
    clog11 = 14.1


print("lol")

if enableDistributionPlot == True:
    for file_i, filename in enumerate(files):
        file = open(filename, "r")
        matches = re.findall(r"step_C_(\d+)\.txt", filename)
        t = int(matches[0])

        # if t!=0 and t!=52000:
        #     continue

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

        if t!=0 and t!=52000:
            continue

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
                plt.contourf(vx[i,:].reshape(n,n), vy[i,:].reshape(n,n), f[i,:].reshape(n,n), resolutions[i], vmax=vmaxes[i], alpha=alphas[i], cmap=colourMaps[i])
            imgname = 'out/data/plot_C_' + str(t+1).zfill(6) + '.png'
            plt.axis('equal')
            plt.ylabel("V [m/s]")
            plt.xlabel("V [m/s]")
            plt.draw()
            print("savefig")
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
        # if file_i > 100:
        #     continue
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
                ms[0, file_i] = 9.1093837*10**(-31)
                ms[1, file_i] = 3.3435837724*10**(-27)
                ns[0, file_i] = 10**20
                ns[1, file_i] = 10**20

        file.close()

    print(dS)
    
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
    charts.plotVAxes(VFlowSpecies)
    # charts.plotVAxesLog(VFlowSpecies)
    # charts.plotPerr(Perr, Pnorm)
    charts.plotEntropy(dS)
    if enableAnalyticFunctions:
        charts.plotTemperatureAnalytic(TspeciesAxes)

    plt.clf()
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.scatter(times, Pnorm, s=0.5)
    plt.plot(times, Pnorm)
    plt.xlabel("Time [s]")
    plt.grid()
    plt.ylabel("|P|")
    plt.savefig("out/P.eps")
    plt.savefig("out/P.png")

    # plot single species energy
    plt.clf()
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    for s in range(0, nspecies):
        plt.scatter(times, Especies[s], label=specieNames[s], s=0.5)
        plt.plot(times, Especies[s])
    legend = plt.legend()
    for s in range(0, nspecies):
        legend.legendHandles[s]._sizes = [30]
    plt.xlabel("Time [s]")
    plt.ylabel("Total Energy [J]")
    plt.grid()
    plt.gcf().subplots_adjust(left=0.15)
    plt.xscale("log")
    plt.ticklabel_format(axis='y', style='sci', scilimits=(3,4))
    plt.savefig("out/EnergySpecies_log.eps")
    plt.savefig("out/EnergySpecies_log.png")

    # plot single species energy
    plt.clf()
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    for s in range(0, nspecies):
        plt.scatter(times, Especies[s], label=specieNames[s], s=0.5)
        plt.plot(times, Especies[s])
    legend = plt.legend()
    for s in range(0, nspecies):
        legend.legendHandles[s]._sizes = [30]
    plt.xlabel("Time [s]")
    plt.ylabel("Total Energy [J]")
    plt.grid()
    plt.gcf().subplots_adjust(left=0.15)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(3,4))
    plt.savefig("out/EnergySpecies.eps")
    plt.savefig("out/EnergySpecies.png")

    # plot single species temeperature separated for each axis (log scale)
    plt.clf()
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    for s in range(0, nspecies):
        plt.scatter(times, TspeciesAxes[s, 0], label=specieNames[s]+" x", s=0.5)
        plt.plot(times, TspeciesAxes[s, 0])
        plt.scatter(times, TspeciesAxes[s, 1], label=specieNames[s]+" y", s=0.5)
        plt.plot(times, TspeciesAxes[s, 1])
    legend = plt.legend()
    for s in range(0, 2*nspecies):
        legend.legendHandles[s]._sizes = [30]
    # plt.ylim(bottom=0)
    plt.xscale("log")
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Temperature [eV]")
    plt.gcf().subplots_adjust(left=0.15)
    current_values = plt.gca().get_yticks()
    plt.gca().set_yticklabels(['{:.0f}'.format(x) for x in current_values])
    plt.savefig("out/TemperatureSpeciesAxesLog.eps")
    plt.savefig("out/TemperatureSpeciesAxesLog.png")


    # plot single species temeperature separated for each axis plt.clf()
    if enableAnalyticFunctions:
        plt.clf()
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        Tanalytic = Tisothermal(times, ma, mb, na, Tax, Tay, Tbx, Tby, clog00, clog01, clog11)
        for s in range(0, nspecies):
            plt.scatter(times, TspeciesAxes[s, 0], label=specieNames[s]+" x", s=0.5)
            plt.plot(times, TspeciesAxes[s, 0])
            plt.scatter(times, TspeciesAxes[s, 1], label=specieNames[s]+" y", s=0.5)
            plt.plot(times, TspeciesAxes[s, 1])
        # analytic solution for electrons
        plt.scatter(times, Tanalytic[0], label="Analytic", s=0.5)
        plt.plot(times, Tanalytic[0], color="green", linestyle="dashed")
        plt.scatter(times, Tanalytic[1], s=0.5)
        plt.plot(times, Tanalytic[1], color="green", linestyle="dashed")
        plt.scatter(times, Tanalytic[2], label="Analytic", s=0.5)
        plt.plot(times, Tanalytic[2], color="green", linestyle="dashed")
        plt.scatter(times, Tanalytic[3], s=0.5)
        plt.plot(times, Tanalytic[3], color="green", linestyle="dashed")
        for s in range(0, 2*nspecies):
            legend.legendHandles[s]._sizes = [30]
        plt.xlim(right=10**(-4), left=0)
        plt.grid()
        plt.xlabel("Time [s]")
        plt.ylabel("Temperature [eV]")
        plt.gcf().subplots_adjust(left=0.15)
        current_values = plt.gca().get_yticks()
        plt.gca().set_yticklabels(['{:.0f}'.format(x) for x in current_values])
        plt.savefig("out/TemperatureSpeciesAxesAnalytic_zoom.eps")
        plt.savefig("out/TemperatureSpeciesAxesAnalytic_zoom.png")

    # plot single species temeperature separated for each axis
    plt.clf()
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    for s in range(0, nspecies):
        plt.scatter(times, TspeciesAxes[s, 0], label=specieNames[s]+" x", s=0.5)
        plt.plot(times, TspeciesAxes[s, 0])
        plt.scatter(times, TspeciesAxes[s, 1], label=specieNames[s]+" y", s=0.5)
        plt.plot(times, TspeciesAxes[s, 1])
    legend = plt.legend()
    for s in range(0, 2*nspecies):
        legend.legendHandles[s]._sizes = [30]
    # plt.ylim(bottom=0)
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Temperature [eV]")
    plt.gcf().subplots_adjust(left=0.15)
    current_values = plt.gca().get_yticks()
    plt.gca().set_yticklabels(['{:.0f}'.format(x) for x in current_values])
    plt.savefig("out/TemperatureSpeciesAxes.eps")
    plt.savefig("out/TemperatureSpeciesAxes.png")

    # plot single species temeperature
    plt.clf()
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    for s in range(0, nspecies):
        plt.scatter(times, Tspecies[s], label=specieNames[s], s=0.5)
        # plt.scatter(times, Tspecies[s], label=specieNames[s], s=0.5)
        plt.plot(times, Tspecies[s])
    legend = plt.legend()
    for s in range(0, nspecies):
        legend.legendHandles[s]._sizes = [30]
    plt.xlabel("Time [s]")
    plt.ylabel("Temperature [eV]")
    plt.gcf().subplots_adjust(left=0.17)
    plt.xscale("log")
    current_values = plt.gca().get_yticks()
    plt.grid()
    # plt.gca().set_yticklabels(['{:.4f}'.format(x) for x in current_values])
    plt.savefig("out/TemperatureSpecies_log.eps")
    plt.savefig("out/TemperatureSpecies_log.png")

    # plot single species temeperature
    plt.clf()
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    for s in range(0, nspecies):
        plt.scatter(times, Tspecies[s], label=specieNames[s], s=0.5)
        # plt.scatter(times, Tspecies[s], label=specieNames[s], s=0.5)
        plt.plot(times, Tspecies[s])
    legend = plt.legend()
    for s in range(0, nspecies):
        legend.legendHandles[s]._sizes = [30]
    plt.xlabel("Time [s]")
    plt.ylabel("Temperature [eV]")
    plt.gcf().subplots_adjust(left=0.17)
    current_values = plt.gca().get_yticks()
    plt.grid()
    # plt.gca().set_yticklabels(['{:.4f}'.format(x) for x in current_values])
    plt.savefig("out/TemperatureSpecies.eps")
    plt.savefig("out/TemperatureSpecies.png")

    # # plot single species temeperature
    plt.clf()
    an = Tthermal(times, ma, mb, na, clog00, Tax, Tay, Tbx, Tby)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    for s in range(0, nspecies):
        plt.scatter(times, Tspecies[s], label=specieNames[s], s=0.5)
        plt.plot(times, Tspecies[s])
    plt.scatter(times, an[0], label="Analytic", s=0.5)
    plt.plot(times, an[0])
    legend = plt.legend()
    for s in range(0, nspecies):
        legend.legendHandles[s]._sizes = [30]
    legend.legendHandles[nspecies]._sizes = [30]
    plt.xlabel("Time [s]")
    plt.ylabel("Temperature [eV]")
    plt.gcf().subplots_adjust(left=0.15)
    current_values = plt.gca().get_yticks()
    plt.grid()
    plt.gca().set_yticklabels(['{:.0f}'.format(x) for x in current_values])
    plt.savefig("out/TemperatureSpeciesAnalytic.eps")
    plt.savefig("out/TemperatureSpeciesAnalytic.png")

if enableDistributionPlot:
    img_array = []
    size = (0,0)
    for filename in glob.glob("./out/data/plot_C*.png"):
        img = cv2.imread(filename)
        height, width, layers = img.shape
        size = (width,height)
        img_array.append(img)

    print(size)
    out = cv2.VideoWriter('./out/out.avi',cv2.VideoWriter_fourcc(*'DIVX'), 30, size)

    for i in range(len(img_array)):
        out.write(img_array[i])

    out.release()
