import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import math
import cv2

files = glob.glob("./out/step_C*.txt")
files.sort()
fmax = 0
for filename in files:
    file = open(filename, "r")
    lines = file.readlines()
    for i in range(1, len(lines)):
        fmax = max(fmax, float(lines[i].strip().split(" ")[4]))
for filename in files:
    print(filename)
    matches = re.findall(r"step_C_(\d+)\.txt", filename)
    t = int(matches[0])
    file = open(filename, "r")
    lines = file.readlines()
    n = int(math.sqrt(float(lines[0].strip().split(" ")[1])))
    nspecies = int(lines[0].strip().split(" ")[0])
    print("nspecies {}".format(nspecies))
    print("n {}".format(n))
    idx = np.zeros([nspecies, n*n])
    vx = np.zeros([nspecies, n*n])
    vy = np.zeros([nspecies, n*n])
    specie = np.zeros([nspecies, n*n])
    f = np.zeros([nspecies, n*n])
    for i in range(1, len(lines)):
        line = lines[i]
        line_s = line.strip().split(" ")
        specie = int(line_s[0])
        index = int(line_s[1])
        vx[specie][index] = float(line_s[2])
        vy[specie][index] = float(line_s[3])
        f[specie][index] = float(line_s[4])

    #plot distribution to file
    plt.cla()

    colourMaps = ["inferno", "viridis", "summer"]
    for i in range(nspecies):
        plt.contourf(vx[i,:].reshape(n,n), vy[i,:].reshape(n,n), f[i,:].reshape(n,n), 20, vmax=fmax, alpha=0.7, cmap=colourMaps[i])
    plt.draw()
    imgname = 'out/plot_C_' + str(t+1).zfill(3) + '.png'
    plt.savefig(imgname)

img_array = []
size = (0,0)
for filename in glob.glob("./out/plot_C*.png"):
    img = cv2.imread(filename)
    height, width, layers = img.shape
    size = (width,height)
    img_array.append(img)

print(size)
out = cv2.VideoWriter('./out/project.avi',cv2.VideoWriter_fourcc(*'DIVX'), 30, size)

for i in range(len(img_array)):
    out.write(img_array[i])

out.release()
