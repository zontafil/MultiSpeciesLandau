import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import math
import cv2

files = glob.glob("./step_C*.txt")
for filename in files:
    print(filename)
    matches = re.findall(r"step_C_(\d+)\.txt", filename)
    t = int(matches[0])
    file = open(filename, "r")
    lines = file.readlines()
    n = int(math.sqrt(len(lines)))
    idx = np.zeros(n*n)
    vx = np.zeros(n*n)
    vy = np.zeros(n*n)
    f = np.zeros(n*n)
    for i in range(len(lines)):
        line = lines[i]
        line_s = line.strip().split(" ")
        idx[i] = int(line_s[0])
        vx[i] = float(line_s[1])
        vy[i] = float(line_s[2])
        f[i] = float(line_s[3])

    #plot distribution to file
    plt.cla()
    plt.contourf(vx.reshape(n,n), vy.reshape(n,n), f.reshape(n,n), 20)
    plt.draw()
    imgname = 'plot_C_' + str(t+1).zfill(3) + '.png'
    plt.savefig(imgname)

img_array = []
size = (0,0)
for filename in glob.glob("./plot_C*.png"):
    img = cv2.imread(filename)
    height, width, layers = img.shape
    size = (width,height)
    img_array.append(img)

print(size)
out = cv2.VideoWriter('project.avi',cv2.VideoWriter_fourcc(*'DIVX'), 15, size)

for i in range(len(img_array)):
    out.write(img_array[i])

out.release()
