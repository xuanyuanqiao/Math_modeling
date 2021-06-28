import math
import matplotlib.pyplot as plt
import random

import numpy as np
from matplotlib import colors


def sigma(x):
    return (1 + math.tanh(2 * x)) / 2

def ReLU(x):
    return max(0, x + 0.5)

def signaling_gradient(position, time):
    # if time >= Taug:
    # return 0
    gaussian_profile = S * sigma(1 - time / Taug) * (math.exp((-position * position) / (2 * L * L)) + math.exp(
        (-(1 - position) * (1 - position)) / (2 * L * L)))
    narrower_profile = sigma(time / Taug - 1) * (math.exp((-position * position) / (2 * l * l)) + math.exp(
        (-(1 - position) * (1 - position)) / (2 * l * l)))
    result = gaussian_profile + narrower_profile
    return result

# compute the signal received by the cell
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

# cell class
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
        self.signal_receive = signal + signal * D * random.gauss(0, 1)

    def signal_output(self):
        state = self.state
        ligand_level = state
        ligand_activity = A0 + 3 * (state * state * state) / (2 + state * state) * A1
        return ligand_level * ligand_activity

    def dynamics(self, interval):
        u = self.state
        s = self.signal_receive
        state_change = ((ReLU(2 * (u - s)) - u + random.gauss(0, 2 * D * Tau**2))/Tau) * interval
        self.state = state_change + u

# parameters
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
D = 5 * (10 ** (-5))

# initialization
back = []
for i in range(width):
    for j in range(height):
        x = (i + 0.5) / width
        x += z * (1 / width) * random.gauss(0, 1)
        y = (j + 0.5) / height
        y += z * (1 / height) * random.gauss(0, 1)
        back.append(Cell(0, x, y, i + 1, j + 1))

developtime = 10.0
step = 0.1
timelist = np.arange(0.0, developtime, step)

tlist = np.arange(0.0, developtime, 0.5)

color_band_1 = colors.LinearSegmentedColormap.from_list('mylist_1', ["darkgreen", "lime"], N=800)
color_band_2 = colors.LinearSegmentedColormap.from_list('mylist_2', ["red", "magenta"], N=800)

#lists used to draw the curve for change of cell state and signal production during the simulation
state_list_7 = []
signal_list_7 = []
state_list_8 = []
signal_list_8 = []
state_list_9 = []
signal_list_9 = []
state_list_10 = []
signal_list_10 = []
state_list_11 = []
signal_list_11 = []

# begin simulation
for t in timelist:
    xlist = []
    ylist = []
    valuelist = []
    state_list = []
    cellcount = 0
    for cell in back:
        signal = signal_to_cell(cell, back)
        cell.receive(t, signal)

        if cell.state < 0.3:
            color = color_band_1(cell.signal_receive)
        else:
            color = color_band_2(cell.state)
        valuelist.append(color)

        if cellcount % 18 == 7:
            state_list_7.append(cell.state)
            signal_list_7.append(signal)
            
        if cellcount % 18 == 8:
            state_list_8.append(cell.state)
            signal_list_8.append(signal)
            
        if cellcount % 18 == 9:
            state_list_9.append(cell.state)
            signal_list_9.append(signal)
            
        if cellcount % 18 == 10:
            state_list_10.append(cell.state)
            signal_list_10.append(signal)
            
        if cellcount % 18 == 11:
            state_list_11.append(cell.state)
            signal_list_11.append(signal)
            

        cellcount += 1

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
        plt.savefig('./figure/improvement_change_function/' + str(round(t, 2)) + '.png')
        plt.cla()
    

#Draw the curve for change of cell state and signal production during the simulation
drawing_state = []
drawing_signal = []
plt.figure(figsize=(16, 16))
plt.ylim([0, 1])
plt.xlim([0, 1])
for i in range(18):
    for j in range(len(state_list_7)):
        if j % 18 == i:
            drawing_state.append(state_list_7[j])
            drawing_signal.append(signal_list_7[j])
    #print(drawing_signal)
    plt.plot(drawing_signal, drawing_state, color='#498EB5',linewidth=0.5,linestyle='-')
    drawing_state = []
    drawing_signal = []
    for j in range(len(state_list_8)):
        if j % 18 == i:
            drawing_state.append(state_list_8[j])
            drawing_signal.append(signal_list_8[j])
    plt.plot(drawing_signal, drawing_state, color='#6250B1',linewidth=0.5,linestyle='-')
    drawing_state = []
    drawing_signal = []
    for j in range(len(state_list_9)):
        if j % 18 == i:
            drawing_state.append(state_list_9[j])
            drawing_signal.append(signal_list_9[j])
    plt.plot(drawing_signal, drawing_state, color='#AC23B0',linewidth=0.5,linestyle='-')
    drawing_state = []
    drawing_signal = []
    for j in range(len(state_list_10)):
        if j % 18 == i:
            drawing_state.append(state_list_10[j])
            drawing_signal.append(signal_list_10[j])
    plt.plot(drawing_signal, drawing_state, color='#6250B1',linewidth=0.5,linestyle='-')
    drawing_state = []
    drawing_signal = []
    for j in range(len(state_list_11)):
        if j % 18 == i:
            drawing_state.append(state_list_11[j])
            drawing_signal.append(signal_list_11[j])
    plt.plot(drawing_signal, drawing_state, color='#498EB5',linewidth=0.5,linestyle='-')
    drawing_state = []
    drawing_signal = []
    
ax=plt.gca();
ax.spines['bottom'].set_linewidth(2);
ax.spines['left'].set_linewidth(2);
ax.spines['right'].set_linewidth(2);
ax.spines['top'].set_linewidth(2);
new_ticks = [0,1]
plt.xticks(new_ticks,fontsize=30)
plt.xlabel("s",fontsize=30)
plt.ylabel("u",fontsize=30)
plt.yticks(new_ticks,fontsize=30 )
plt.savefig('u-s.png')
plt.cla()

drawing_state = []
plt.figure(figsize=(16, 16))
plt.ylim([0, 1])
plt.xlim([0, 10])
for i in range(18):
    for j in range(len(state_list_7)):
        if j % 18 == i:
            drawing_state.append(state_list_7[j])
    #print(drawing_signal)
    plt.plot(timelist, drawing_state, color='#498EB5',linewidth=0.5,linestyle='-')
    drawing_state = []
    drawing_signal = []
    for j in range(len(state_list_8)):
        if j % 18 == i:
            drawing_state.append(state_list_8[j])
    plt.plot(timelist, drawing_state, color='#6250B1',linewidth=0.5,linestyle='-')
    drawing_state = []
    drawing_signal = []
    for j in range(len(state_list_9)):
        if j % 18 == i:
            drawing_state.append(state_list_9[j])
    plt.plot(timelist, drawing_state, color='#AC23B0',linewidth=0.5,linestyle='-')
    drawing_state = []
    drawing_signal = []
    for j in range(len(state_list_10)):
        if j % 18 == i:
            drawing_state.append(state_list_10[j])
    plt.plot(timelist, drawing_state, color='#6250B1',linewidth=0.5,linestyle='-')
    drawing_state = []
    drawing_signal = []
    for j in range(len(state_list_11)):
        if j % 18 == i:
            drawing_state.append(state_list_11[j])
    plt.plot(timelist, drawing_state, color='#498EB5',linewidth=0.5,linestyle='-')
    drawing_state = []
    drawing_signal = []

ax=plt.gca();
ax.spines['bottom'].set_linewidth(2);
ax.spines['left'].set_linewidth(2);
ax.spines['right'].set_linewidth(2);
ax.spines['top'].set_linewidth(2);

new_ticks = [0,1]
plt.xlabel("t",fontsize=30)
plt.ylabel("u",fontsize=30)
plt.xticks(fontsize=30)
plt.yticks(new_ticks,fontsize=30)
plt.savefig('cis/u-t.png')
plt.cla()


