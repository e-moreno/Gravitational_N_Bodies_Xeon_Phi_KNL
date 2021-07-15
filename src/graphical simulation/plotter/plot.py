#!/usr/bin/python3
from mpl_toolkits.mplot3d import Axes3D
import random
import time
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation
import pandas as pd
import sys


bodiesXData = []
bodiesYData = []
bodiesZData = []

if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    filename = "bodies.csv"

with open(filename,"r") as csvfile:
    firstline = csvfile.readline()
    params = firstline.split(",")
    N = int(params[0].strip())
    iter = int(params[1].strip())
    for i in range(0,iter):
        x_values = []
        y_values = []
        z_values = []
        for j in range(0,N):
            line = csvfile.readline()
            data = line.split(",")
            x_values.append(float(data[0].strip()))
            y_values.append(float(data[1].strip()))
            z_values.append(float(data[2].strip()))
        temp_x = np.array(x_values)
        temp_y = np.array(y_values)
        temp_z = np.array(z_values)

        bodiesXData.append(temp_x)
        bodiesYData.append(temp_y)
        bodiesZData.append(temp_z)


a = np.random.rand(2000, 3)*10
t = np.array([np.ones(100)*i for i in range(20)]).flatten()
df = pd.DataFrame({"time": t ,"x" : a[:,0], "y" : a[:,1], "z" : a[:,2]})



def update_graph(num):
    graph.set_data (bodiesXData[num], bodiesYData[num])
    graph.set_3d_properties(bodiesZData[num])
    title.set_text('3D Test, time={}'.format(num))
    return title, graph,


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
title = ax.set_title('3D Test')

graph, = ax.plot(bodiesXData[0], bodiesYData[0], bodiesZData[0], linestyle="", marker="o")

ani = matplotlib.animation.FuncAnimation(fig, update_graph, iter,
                               interval=40, blit=True)

plt.show()
