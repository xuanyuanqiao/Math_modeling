import math
import matplotlib.pyplot as plt
import random

import numpy as np
from matplotlib import colors


def sigma(x):
    return (1 + math.tanh(2 * x)) / 2


def included_angel(v1, v2):
    v1_norm2 = math.sqrt(v1[0]**2 + v1[1]**2)
    v2_norm2 = math.sqrt(v2[0]**2 + v2[1]**2)
    dot_product = v1[0]*v2[0] + v1[1]*v2[1]
    return np.arccos(dot_product / (v1_norm2 * v2_norm2))


def signaling_gradient(position, time):
    gaussian_profile = S_0 * sigma(1 - time / Taug) * (math.exp((-position * position) / (2 * L * L)) + math.exp(
        (-(1 - position) * (1 - position)) / (2 * L * L)))
    narrower_profile = N_0 * sigma(time / Taug - 1) * (math.exp((-position * position) / (2 * l * l)) + math.exp(
        (-(1 - position) * (1 - position)) / (2 * l * l)))
    result = gaussian_profile + narrower_profile
    return result


def signal_to_cell(target, area):
    total_contribute = 0
    target_x = target.xindex
    target_y = target.yindex
    surrounding_cell_1 = []
    for sur_row in range(target_y-1, target_y+2):
        for sur_col in range(target_x-1, target_x+2):
            if sur_row < 0 or sur_col < 0 or sur_row >= len(area) or sur_col >= len(area[sur_row]):
                continue
            if target.xindex == area[sur_row][sur_col].xindex and target.yindex == area[sur_row][sur_col].yindex:
                continue
            surrounding_cell_1.append(area[sur_row][sur_col])
    angle_list = np.array([0.]*len(surrounding_cell_1))
    for cell_index in range(len(surrounding_cell_1)):
        angle_list[cell_index] = included_angel([surrounding_cell_1[cell_index].xposition - target.xposition, surrounding_cell_1[cell_index].yposition - target.yposition],
                                  [surrounding_cell_1[(cell_index+1)%len(surrounding_cell_1)].xposition - target.xposition, surrounding_cell_1[(cell_index+1)%len(surrounding_cell_1)].yposition - target.yposition])
    angle_list = angle_list / angle_list.sum()

    for cell_index in range(len(surrounding_cell_1)):
        angle = (angle_list[cell_index-1] + angle_list[(cell_index+1)%len(angle_list)])/2
        distance = math.sqrt(
            (target.xposition - surrounding_cell_1[cell_index].xposition) * (target.xposition - surrounding_cell_1[cell_index].xposition) +
            (target.yposition - surrounding_cell_1[cell_index].yposition) * (target.yposition - surrounding_cell_1[cell_index].yposition))
        contribute = math.exp(-(distance/math.exp(angle))**2 / (2 * l_1 * l_1))
        signal = surrounding_cell_1[cell_index].signal_output()
        total_contribute += signal * contribute

    surrounding_cell_2 = []
    for sur_row in range(target_y-2, target_y+3):
        for sur_col in range(target_x-2, target_x+3):
            if sur_row < 0 or sur_col < 0 or sur_row >= len(area) or sur_col >= len(area[sur_row]):
                continue
            if abs(target.xindex - area[sur_row][sur_col].xindex) < 2 and abs(target.yindex - area[sur_row][sur_col].yindex) < 2:
                continue
            surrounding_cell_2.append(area[sur_row][sur_col])

    for cell_index in range(len(surrounding_cell_2)):
            distance = math.sqrt(
                (target.xposition - surrounding_cell_2[cell_index].xposition) * (
                            target.xposition - surrounding_cell_2[cell_index].xposition) +
                (target.yposition - surrounding_cell_2[cell_index].yposition) * (
                            target.yposition - surrounding_cell_2[cell_index].yposition))
            contribute = math.exp(-(distance/math.exp(1/len(surrounding_cell_2)))**2 / (2 * l_2 * l_2))
            signal = surrounding_cell_2[cell_index].signal_output()
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
        state_change = ((sigma(2 * (u - s)) - u + random.gauss(0, 2 * D * Tau**2))/Tau) * interval
        self.state = state_change + u


Tau = 0.5
Taug = 1
S_0 = 2
N_0 = 2
L = 0.2
width = 18
height = 18
N = width * height
lamda = math.sqrt(1 / N)
l = 1.75 * lamda
l_1 = l # for surrounding cell 1
l_2 = 2*l # for surrounding cell 2
A0 = 0.05
A1 = 1 - A0
z = 0.2
D = 5 * (10 ** (-5))

back = []
for i in range(height):
    row = []
    for j in range(width):
        x = (j + 0.5) / width
        x += z * (1 / width) * random.gauss(0, 1)
        y = (i + 0.5) / height
        y += z * (1 / height) * random.gauss(0, 1)
        row.append(Cell(0, x, y, j, i))
    back.append(row)

developtime = 10.0
step = 0.1
timelist = np.arange(0.0, developtime, step)

tlist = np.arange(0.0, developtime, 0.1)

color_band_1 = colors.LinearSegmentedColormap.from_list('mylist_1', ["darkgreen", "lime"], N=800)
color_band_2 = colors.LinearSegmentedColormap.from_list('mylist_2', ["red", "magenta"], N=800)

for t in timelist:
    xlist = []
    ylist = []
    valuelist = []
    state_list = []
    for row in back:
        for cell in row:
            signal = signal_to_cell(cell, back)
            cell.receive(t, signal)

            if cell.state < 0.3:
                color = color_band_1(cell.signal_receive)
            else:
                color = color_band_2(cell.state)
            valuelist.append(color)

    for row in back:
        for cell in row:
            if cell.xindex <= 7 and 2 <= cell.yindex <= 8 and (cell.yindex - cell.xindex) <= 6 and (
                    cell.yindex - cell.xindex) >= 4:
                cell.signal_receive = 0
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
        plt.savefig('notch_mutant/' + str(round(t, 2)) + '.png')
        plt.cla()



