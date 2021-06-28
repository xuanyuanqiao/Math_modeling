#!/usr/bin/env python3
#tissue_size, N=324

import math
import matplotlib.pyplot as plt
import random
import numpy as np
from matplotlib import colors


def sigma(x):
    return (1 + math.tanh(2 * x)) / 2


def signaling_gradient(position, time):
    gaussian_profile = S * sigma(1 - time / Taug) * (math.exp((-position * position) / (2 * L * L)) + math.exp(
        (-(1 - position) * (1 - position)) / (2 * L * L)))
    narrower_profile = sigma(time / Taug - 1) * (math.exp((-position * position) / (2 * l * l)) + math.exp(
        (-(1 - position) * (1 - position)) / (2 * l * l)))
    result = gaussian_profile + narrower_profile
    return result


def signal_to_cell(target, area):
    total_contribute = 0
    for cell in area:
        distance = math.sqrt((target.xposition - cell.xposition) * (target.xposition - cell.xposition) + (
                    target.yposition - cell.yposition) * (target.yposition - cell.yposition))
        if target.xindex == cell.xindex and target.yindex == cell.yindex:
            contribute = 0
        else:
            contribute = math.exp((-distance * distance / (2 * l * l)))
        signal = cell.signal_output()
        total_contribute += signal * contribute
    return total_contribute


class Cell:
    def __init__(self, startstate, xposition, yposition, xindex, yindex):
        self.state = startstate
        self.xposition = xposition
        self.yposition = yposition
        self.xindex = xindex
        self.yindex = yindex
        self.signal_receive = 0
        
    def receive(self, t, neighboringinput):
        x = self.xposition
        signal = signaling_gradient(x, t) + neighboringinput
        self.signal_receive = signal 
        
    def signal_output(self):
        state = self.state
        ligand_level = state
        ligand_activity = A0 + 3 * (state * state * state) / (2 + state * state) * A1
        return ligand_level * ligand_activity
    
    def dynamics(self, interval):
        u = self.state
        s = self.signal_receive
        state_change = ((sigma(2 * (u - s)) - u + random.gauss(0, 2 * D * Tau ** 2)) / Tau) * interval
        self.state = state_change + u
        
        
Tau = 0.5
Taug = 1
S = 2
L = 0.2
width = 18
height = 18
N = width * height
lamda = math.sqrt(1 / N)
l = 1.75 * lamda
A0 = 0.05
A1 = 1 - A0
z = 0.2
D = 5e-5

back = []
for i in range(width):
    for j in range(height):
        x = (i + 0.5) / width
        x += z * (1 / width) * random.gauss(0, 1)
        y = (j + 0.5) / height
        y += z * (1 / height) * random.gauss(0, 1)
        back.append(Cell(0, x, y, i + 1, j + 1))
        
developtime = 10.5
step = 0.1
timelist = np.arange(0.0, developtime, step)

tlist = np.arange(0.0, developtime, 1.0)

color_band_1 = colors.LinearSegmentedColormap.from_list('mylist_1', ["darkgreen", "lime"], N=800)
color_band_2 = colors.LinearSegmentedColormap.from_list('mylist_2', ["red", "magenta"], N=800)

for t in timelist:
    xlist = []
    ylist = []
    valuelist = []
    state_list = []
    for cell in back:
        signal = signal_to_cell(cell, back)
        cell.receive(t, signal)
        
        if cell.state < 0.3:
            color = color_band_1(cell.signal_receive)
        else:
            color = color_band_2(cell.state)
        valuelist.append(color)
        
    for cell in back:
        cell.dynamics(step)
        xlist.append(cell.xposition)
        ylist.append(cell.yposition)
        state_list.append(cell.state)
        
    if t in tlist:
        ax = plt.figure(facecolor='k')
        plt.axis("off")
        ax.set_figwidth(10)
        ax.set_figheight(10)
        plt.scatter(xlist, ylist, c=valuelist, s=300)
        plt.savefig('./notch_perturbation/tissue_size/324/' +str(round(t, 2)) + '.png')
        plt.cla()
        
        
        
        