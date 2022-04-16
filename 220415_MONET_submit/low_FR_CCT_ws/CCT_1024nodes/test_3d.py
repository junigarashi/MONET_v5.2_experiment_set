import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()

f = open('FR_M1_P1_subgrid_with_time.dat')
raw_data = f.readlines()
f.close()

steps = len(raw_data)

n_sg = 5
sg = []

for i_rd, rd in enumerate(raw_data):
    sg.append([])
    line = rd.split(" ")
    line.remove("\n")
    for word in line[1:24]:
        sg[i_rd].append(float(word))

        
for step in range(steps):
    print step 
    ax = Axes3D(fig)
    rows, columns, heights, rgb, size = [], [], [], [], []
    for i_sg, value in enumerate(sg[step]):
        rows.append(i_sg % n_sg)
        columns.append(i_sg / n_sg)
        heights.append(0)
        rgb.append([value / 1000., 0, 0, 1])
        size.append(500)
        
    ax.set_aspect('equal')
    ax.scatter(rows, columns, heights, s=size, c=rgb, edgecolor=rgb)
    ax.view_init(90.0, 0.0)

    plt.pause(.01)

