import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import math
import cv2
import os
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as mtick

plt.figure(figsize=(6,6))
plt.rcParams.update({'font.size': 16})

ma = 9.109*10**(-31)
mb = 3.34*10**(-27)
density = 1*10**(20)
l = 6.6
Tax = 400
Tay = 300
Tbx = 250
Tby = 200

def createDir(path):
    try: 
        os.mkdir(path)
    except OSError as error:
        pass

def time_thermalization(t, ma, mb, n, l, Ta, Tb):
    ma = ma * 1000
    mb = mb * 1000
    n = n * 10**(-6)
    nu =  1.8*10**(-19) * np.sqrt(ma*mb) * n * l / (ma*Ta + mb*Tb)**(1.5)
    tau = 1. / (2*nu)
    return tau

def time_isotropization(t, m, n, l, Tx, Ty):
    Tpar = Tx
    Tperp = Ty
    e = 1.602176565*10**(-19) * 3*10**9
    m = m * 1000
    n = n * 10**(-6)
    k = 1.6*10**(-12)
    A = Tperp/Tpar - 1
    if (A>0):
        f = np.arctan(np.sqrt(A))/np.sqrt(A)
    else:
        f = np.arctanh(np.sqrt(-A))/np.sqrt(-A)
    nu = (2*np.sqrt(np.pi)*e**4*n*l)/(np.sqrt(m)*(k*Tpar)**(1.5)) \
        * A**(-2) * (-3 + (A+3)*(f))
    tau = 1 / nu
    ret = [tau/4, tau/4]
    # ret[0] =  (Tpar+Tperp)/2. + (Tpar-Tperp)/2.* np.exp(-4.*t/tau)
    # ret[1] =  (Tperp+Tpar)/2. + (Tperp-Tpar)/2.* np.exp(-2.*t/tau)
    return ret

def Tisothermal(t, ma, mb, n, l, Tax, Tay, Tbx, Tby):
    # compute times
    Ta = np.sqrt(Tax*Tay)
    Tb = np.sqrt(Tbx*Tby)
    Ta = (Tax+Tay)/2
    Tb = (Tbx+Tby)/2
    tau_a_iso = time_isotropization(t, ma, n, l, Tax, Tay)
    tau_b_iso = time_isotropization(t, mb, n, l, Tbx, Tby)
    tau_a_thermal = time_thermalization(t, ma, mb, n, l, Ta, Tb)
    tau_b_thermal = time_thermalization(t, mb, ma, n, l, Tb, Ta)
    
    ret = [0,0,0,0]
    ret[0] = (Ta+Tb)/2. + (Ta-Tb)/2.* np.exp(-t/tau_a_thermal) + (Tax-Tay)/2 * np.exp(-t/tau_a_iso[0])
    ret[1] = (Ta+Tb)/2. + (Ta-Tb)/2.* np.exp(-t/tau_a_thermal) + (Tay-Tax)/2 * np.exp(-t/tau_a_iso[1])
    ret[2] = (Ta+Tb)/2. + (Tb-Ta)/2.* np.exp(-t/tau_b_thermal) + (Tbx-Tby)/2 * np.exp(-t/tau_b_iso[0])
    ret[3] = (Ta+Tb)/2. + (Tb-Ta)/2.* np.exp(-t/tau_b_thermal) + (Tby-Tbx)/2 * np.exp(-t/tau_b_iso[1])

    return ret

def Tthermal(t, ma, mb, n, l, Tax, Tay, Tbx, Tby):
    # compute times
    Ta = np.sqrt(Tax*Tay)
    Ta = (Tax+Tay)/2
    Tb = np.sqrt(Tbx*Tby)
    Tb = (Tbx+Tby)/2
    tau_a_thermal = time_thermalization(t, ma, mb, n, l, Ta, Tb)
    tau_b_thermal = time_thermalization(t, mb, ma, n, l, Tb, Ta)
    
    ret = [0,0]
    ret[0] = (Ta+Tb)/2. + (Ta-Tb)/2.* np.exp(-t/tau_a_thermal)
    ret[1] = (Ta+Tb)/2. + (Tb-Ta)/2.* np.exp(-t/tau_b_thermal)

    return ret


files = glob.glob("./out/data/step_C*.txt")
files.sort()

# files = files[0:100]

fmax = None
nspecies = 0
dt = 0
file = open(files[0], "r")
lines = file.readlines()
nspecies = int(lines[0].strip().split(" ")[0])
fmax = np.zeros(nspecies)
for file_i, filename in enumerate(files):
    print("Analizing file {} of {}".format(file_i, len(files)))
    file = open(filename, "r")
    lines = file.readlines()
    n = int(math.sqrt(float(lines[0].strip().split(" ")[2])))
    nspecies = int(lines[0].strip().split(" ")[0])
    dt = float(lines[0].strip().split(" ")[3])
    if len(lines) > 2+nspecies:
        for i in range(2+nspecies, len(lines)):
            line = lines[i]
            line_s = line.strip().split(" ")
            if (len(line_s)<5):
                continue
            specie = int(line_s[0])
            index = int(line_s[1])
            fmax[specie] = max(fmax[specie], float(line_s[4]))

E = np.zeros(len(files))
Eerr = np.zeros(len(files))
P = np.zeros((len(files), 2))
Pnorm = np.zeros(len(files))
PSpeciesNorm = np.zeros((nspecies, len(files)))
Perr = np.zeros((len(files), 2))
times = np.zeros(len(files))
thermalAnalytic = np.zeros(len(files))
Especies = np.zeros((nspecies, len(files)))
Tspecies = np.zeros((nspecies, len(files)))
TspeciesAxes = np.zeros((nspecies, 2, len(files)))
Pspecies = np.zeros((nspecies, 2, len(files)))
specieNames = []

for file_i, filename in enumerate(files):
    matches = re.findall(r"step_C_(\d+)\.txt", filename)
    t = int(matches[0])
    file = open(filename, "r")
    lines = file.readlines()
    n = int(math.sqrt(float(lines[0].strip().split(" ")[2])))
    nMarkersPerDim = int(math.sqrt(float(lines[0].strip().split(" ")[1])))
    idx = np.zeros([nspecies, n*n])
    vx = np.zeros([nspecies, n*n])
    vy = np.zeros([nspecies, n*n])
    specie = np.zeros([nspecies, n*n])
    f = np.zeros([nspecies, n*n])

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
    times[file_i] = t
    for s in range(nspecies):
        line = lines[2+s].strip().split(" ")
        Especies[s, file_i] = float(line[0])
        Pspecies[s, 0, file_i] = float(line[1])
        Pspecies[s, 1, file_i] = float(line[2])
        PSpeciesNorm[s, file_i] = np.sqrt(Pspecies[s,0,file_i]**2 + Pspecies[s,1,file_i]**2)
        Tspecies[s, file_i] = float(line[3])
        TspeciesAxes[s, 0, file_i] = float(line[4])
        TspeciesAxes[s, 1, file_i] = float(line[5])
        specieNames.append(line[6])

    # debug print
    if file_i == 0:
        print("nspecies {}".format(nspecies))
        print("n {}".format(n))
        print("n Markers/Dim {}".format(nMarkersPerDim))
    if file_i % 20 == 0:
        print(filename)

    # record distribution
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
        vmaxes = [fmax[0]/2]
        if nspecies > 1:
            vmaxes.append(fmax[1])
        for i in range(nspecies):
            plt.contourf(vx[i,:].reshape(n,n), vy[i,:].reshape(n,n), f[i,:].reshape(n,n), resolutions[i], vmax=vmaxes[i], alpha=alphas[i], cmap=colourMaps[i])
        imgname = 'out/data/plot_C_' + str(t+1).zfill(6) + '.png'
        plt.axis('equal')
        plt.ylabel("V [m/s]")
        plt.xlabel("V [m/s]")
        plt.draw()
        plt.savefig(imgname)

createDir("./plots/eps")
createDir("./plots/png")

times = times * dt

# plot system momentum and energy
plt.clf()
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.scatter(times, Eerr, s=0.5)
plt.plot(times, Eerr)
plt.xlabel("Time [s]")
plt.grid()
plt.ylabel("(E-E0)/E0")
plt.gcf().subplots_adjust(left=0.15)
plt.savefig("out/EnergyError.eps")
plt.savefig("out/EnergyError.png")

# plot system momentum and energy
# Perr = Perr / P
PerrNorm = np.sqrt(Perr[:,0]**2 + Perr[:,1]**2)
PerrNorm = PerrNorm / Pnorm
plt.clf()
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.scatter(times, PerrNorm, s=0.5)
plt.plot(times, PerrNorm)
plt.xlabel("Time [s]")
plt.ylabel("|P - P0|/|P0|")
plt.grid()
plt.gcf().subplots_adjust(left=0.15)
plt.savefig("out/PerrNorm.eps")
plt.savefig("out/PerrNorm.png")

plt.clf()
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.scatter(times, Pnorm, s=0.5)
plt.plot(times, Pnorm)
plt.xlabel("Time [s]")
plt.grid()
plt.ylabel("|P|")
plt.savefig("out/P.eps")
plt.savefig("out/P.png")

# plot single species momentum separated for each axis
plt.clf()
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,2))
for s in range(0, nspecies):
    plt.scatter(times, Pspecies[s, 0], label=specieNames[s]+" x", s=0.5)
    plt.plot(times, Pspecies[s, 0])
    plt.scatter(times, Pspecies[s, 1], label=specieNames[s]+" y", s=0.5)
    plt.plot(times, Pspecies[s, 1])
legend = plt.legend()
for s in range(0, 2*nspecies):
    legend.legendHandles[s]._sizes = [30]
# plt.ylim(bottom=0)
plt.xlabel("Time [s]")
plt.ylabel("P [kg m/s]")
plt.grid()
plt.gcf().subplots_adjust(left=0.15)
current_values = plt.gca().get_yticks()
# plt.gca().set_yticklabels(['{:.0f}'.format(x) for x in current_values])
plt.savefig("out/Vaxes.eps")
plt.savefig("out/Vaxes.png")

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
plt.clf()
times_log = times.copy()
# times_log[0] = 10**(-7)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
Tanalytic = Tisothermal(times_log, ma, mb, density, l, Tax, Tay, Tbx, Tby)
for s in range(0, nspecies):
    plt.scatter(times_log, TspeciesAxes[s, 0], label=specieNames[s]+" x", s=0.5)
    plt.plot(times_log, TspeciesAxes[s, 0])
    plt.scatter(times_log, TspeciesAxes[s, 1], label=specieNames[s]+" y", s=0.5)
    plt.plot(times_log, TspeciesAxes[s, 1])
# analytic solution for electrons
plt.scatter(times_log, Tanalytic[0], color="black", label="Analytic", s=0.5)
plt.plot(times_log, Tanalytic[0], color="black", linestyle="dashed")
plt.scatter(times_log, Tanalytic[1], s=0.5)
plt.plot(times_log, Tanalytic[1], color="black", linestyle="dashed")
plt.scatter(times_log, Tanalytic[2], s=0.5)
plt.plot(times_log, Tanalytic[2], color="black", linestyle="dashed")
plt.scatter(times_log, Tanalytic[3], s=0.5)
plt.plot(times_log, Tanalytic[3], color="black", linestyle="dashed")
# legend = plt.legend(loc='upper right')
legend = plt.legend(bbox_to_anchor=(1,1), loc="upper right")

for s in range(0, 5):
    legend.legendHandles[s]._sizes = [30]
# plt.ylim(bottom=0)
plt.grid()
plt.xscale("log")
plt.xlabel("Time [s]")
plt.ylabel("Temperature [eV]")
plt.gcf().subplots_adjust(left=0.15)
current_values = plt.gca().get_yticks()
plt.gca().set_yticklabels(['{:.0f}'.format(x) for x in current_values])
plt.savefig("out/TemperatureSpeciesAxesAnalytic_log.eps")
plt.savefig("out/TemperatureSpeciesAxesAnalytic_log.png")


# plot single species temeperature separated for each axis plt.clf()
plt.clf()
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
Tanalytic = Tisothermal(times, ma, mb, density, l, Tax, Tay, Tbx, Tby)
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
    plt.plot(times, Tspecies[s])
legend = plt.legend()
for s in range(0, nspecies):
    legend.legendHandles[s]._sizes = [30]
plt.xlabel("Time [s]")
plt.ylabel("Temperature [eV]")
plt.gcf().subplots_adjust(left=0.15)
current_values = plt.gca().get_yticks()
plt.grid()
# plt.gca().set_yticklabels(['{:.4f}'.format(x) for x in current_values])
plt.savefig("out/TemperatureSpecies.eps")
plt.savefig("out/TemperatureSpecies.png")

# plot single species temeperature
plt.clf()
an = Tthermal(times, ma, mb, density, l, Tax, Tay, Tbx, Tby)
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
