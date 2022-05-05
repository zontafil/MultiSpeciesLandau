import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import math
import cv2
import os

plt.figure(figsize=(6,6))

def createDir(path):
    try: 
        os.mkdir(path)
    except OSError as error:
        pass


files = glob.glob("./out/data/step_C*.txt")
files.sort()
fmax = None
nspecies = 0
dt = 0
file = open(files[0], "r")
lines = file.readlines()
nspecies = int(lines[0].strip().split(" ")[0])
fmax = np.zeros(nspecies)
for filename in files:
    file = open(filename, "r")
    lines = file.readlines()
    n = int(math.sqrt(float(lines[0].strip().split(" ")[1])))
    nspecies = int(lines[0].strip().split(" ")[0])
    dt = float(lines[0].strip().split(" ")[2])
    for s in range(nspecies):
        for i in range(2+nspecies+s*n*n, s*n*n+n*n):
            fmax[s] = max(fmax[s], float(lines[i].strip().split(" ")[4]))

E = np.zeros(len(files))
Eerr = np.zeros(len(files))
distmin = np.zeros(len(files))
P = np.zeros((len(files), 2))
Pnorm = np.zeros(len(files))
PSpeciesNorm = np.zeros((nspecies, len(files)))
Perr = np.zeros((len(files), 2))
times = np.zeros(len(files))
Especies = np.zeros((nspecies, len(files)))
Tspecies = np.zeros((nspecies, len(files)))
Pspecies = np.zeros((nspecies, 2, len(files)))
specieNames = []

for file_i, filename in enumerate(files):
    matches = re.findall(r"step_C_(\d+)\.txt", filename)
    t = int(matches[0])
    file = open(filename, "r")
    lines = file.readlines()
    n = int(math.sqrt(float(lines[0].strip().split(" ")[1])))
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
    distmin[file_i] = float(line[6])
    times[file_i] = t
    for s in range(nspecies):
        line = lines[2+s].strip().split(" ")
        Especies[s, file_i] = float(line[0])
        Pspecies[s, 0, file_i] = float(line[1])
        Pspecies[s, 1, file_i] = float(line[2])
        PSpeciesNorm[s, file_i] = np.sqrt(Pspecies[s,0,file_i]**2 + Pspecies[s,1,file_i]**2)
        Tspecies[s, file_i] = float(line[3])
        specieNames.append(line[4])

    # debug print
    if file_i == 0:
        print("nspecies {}".format(nspecies))
        print("n {}".format(n))
    if file_i % 20 == 0:
        print(filename)

    # record distribution
    for i in range(2+nspecies, len(lines)):
        line = lines[i]
        line_s = line.strip().split(" ")
        specie = int(line_s[0])
        index = int(line_s[1])
        vx[specie][index] = float(line_s[2])
        vy[specie][index] = float(line_s[3])
        f[specie][index] = float(line_s[4])

    #plot distribution to file
    plt.cla()

    # plot distribution for current time step
    colourMaps = ["summer", "inferno", "viridis", "summer"]
    alphas = [0.7, 0.7]
    resolutions = [20, 20]
    vmaxes = [fmax[0]/10, fmax[1]/2]
    for i in range(nspecies):
        plt.contourf(vx[i,:].reshape(n,n), vy[i,:].reshape(n,n), f[i,:].reshape(n,n), resolutions[i], vmax=vmaxes[i], alpha=alphas[i], cmap=colourMaps[i])
    imgname = 'out/data/plot_C_' + str(t+1).zfill(6) + '.png'
    plt.axis('equal')
    plt.draw()
    plt.savefig(imgname)

createDir("./plots/eps")
createDir("./plots/png")

# plot system momentum and energy
times = times * dt
plt.clf()
plt.scatter(times, Eerr, s=0.5)
plt.plot(times, Eerr)
plt.xlabel("Time [s]")
plt.ylabel("Energy Error")
plt.savefig("out/EnergyError.eps")
plt.savefig("out/EnergyError.png")

plt.clf()
plt.scatter(times, distmin, s=0.5)
plt.plot(times, distmin)
plt.xlabel("Time [s]")
plt.ylabel("Minimum velocity distance")
plt.savefig("out/minimumdist.eps")
plt.savefig("out/minimumdist.png")

plt.clf()
plt.scatter(times, Pnorm, s=0.5)
plt.plot(times, Pnorm)
plt.xlabel("Time [s]")
plt.ylabel("|P|")
plt.savefig("out/P.eps")
plt.savefig("out/P.png")

# plot single species momentum and energy
plt.clf()
for s in range(0, nspecies):
    plt.scatter(times, Especies[s], label=specieNames[s], s=0.5)
    plt.plot(times, Especies[s])
    plt.legend()
plt.xlabel("Time [s]")
plt.ylabel("Energy [J]")
plt.savefig("out/EnergySpecies.eps")
plt.savefig("out/EnergySpecies.png")

# plot single species temeperature
plt.clf()
for s in range(0, nspecies):
    plt.scatter(times, PSpeciesNorm[s], label=specieNames[s], s=0.5)
    plt.plot(times, PSpeciesNorm[s])
    plt.legend()
plt.xlabel("Time [s]")
plt.ylabel("|P|")
plt.savefig("out/PspeciesNorm.eps")
plt.savefig("out/PspeciesNorm.png")

# plot single species temeperature
plt.clf()
for s in range(0, nspecies):
    plt.scatter(times, Tspecies[s], label=specieNames[s], s=0.5)
    plt.plot(times, Tspecies[s])
plt.legend()
plt.ylim(bottom=0)
plt.xlabel("Time [s]")
plt.ylabel("Temperature [eV]")
plt.savefig("out/TemperatureSpecies.eps")
plt.savefig("out/TemperatureSpecies.png")

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
