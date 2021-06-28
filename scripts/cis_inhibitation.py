import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
import random

def sigma(x):
    return (1 + math.tanh(2 * x)) / 2

# compute the initial delta production
def delta_initial(position, time):
    L_prime = math.sqrt(L * L - l * l)
    result = S_0 * sigma(1 - time / Taug) * (math.exp((-position * position) / (2 * L_prime * L_prime)) + math.exp(
        (-(1 - position) * (1 - position)) / (2 * L_prime * L_prime)))
    return result

# delta production
def delta_production(position, time, state):
    result = delta_initial(position, time) + state
    return result

def compute_signal(cell, t):
    state = cell.state
    ligand_activity = A0 + 3 * ((state ** 3) / (2 + state ** 2)) * A1
    delta = delta_production(cell.xposition, t, state)
    return delta * ligand_activity

# compute the signal received by the cell
def signal_to_cell(target, area, t):
    total_contribute = 0
    for cell in area:
        distance = math.sqrt((target.xposition - cell.xposition) ** 2 + (
                target.yposition - cell.yposition) ** 2)
        if target.cnum == cell.cnum and target.rnum == cell.rnum:
            contribute = 0
        else:
            contribute = math.exp((-(distance * distance) / (2 * l * l)))
        signal = compute_signal(cell, t)
        total_contribute += signal * contribute
    return total_contribute

def cis_inhibition(cell, t):
    state = cell.state
    x = cell.xposition
    delta = delta_production(x, t, state)
    notch = (beta_n - delta) / (2 * gama_n) - gama_d / (2 * k_c) + math.sqrt(
        ((beta_n - delta) / (2 * gama_n) - gama_d / (2 * k_c)) ** 2 + (gama_d * beta_n) / (gama_n * k_c))
    return notch

class Cell:
    def __init__(self, startstate, xposition, yposition, cnum, rnum):
        self.state = startstate
        self.xposition = xposition
        self.yposition = yposition
        self.cnum = cnum
        self.rnum = rnum
        self.signal_receive = 0
        self.notch = 0

    # Signal received by the cell
    def receive(self, signal):
        self.signal_receive = k_t * signal

    # Change the state of cell
    def dynamics(self, interval):
        x = self.xposition
        u = self.state
        s = self.signal_receive
        N = self.notch
        state_change = (((sigma(1.4 * u - s * N) - u) + random.gauss(0, 2 * D * Tau**2)) / Tau) * interval
        self.state = state_change + u

# parameters
width = 18
height = 18
beta_n = 2
gama_n = 1
gama_d = 1
k_c = 100
k_t = 1
Tau = 1 / 3
L = 0.2
N = width ** 2
S_0 = 2
lamda = math.sqrt(1 / N)
l = 1.75 * lamda
A0 = 0.05
A1 = 1 - A0
Taug = 1
D = 5e-5

# initialization;
back = []
for i in range(width):
    for j in range(height):
        x = (i + 0.5) / width
        x += L * (1 / width) * random.gauss(0, 1)
        y = (j + 0.5) / height
        y += L * (1 / height) * random.gauss(0, 1)
        back.append(Cell(0, x, y, i, j))

developtime = 10.0
step = 0.1
timelist = np.arange(0.0, developtime, step)
tlist = np.arange(0.0, developtime, 0.1)

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
    cellcount = 0
    for cell in back:
        signal = signal_to_cell(cell, back, t)
        cell.receive(signal)

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
        cell.notch = cis_inhibition(cell, t)
        cell.dynamics(step)
        xlist.append(cell.xposition)
        ylist.append(cell.yposition)

    if t in tlist:
        ax = plt.figure(facecolor='k')
        plt.axis("off")
        ax.set_figwidth(10)
        ax.set_figheight(10)
        plt.scatter(xlist, ylist, c=valuelist, s=300)
        plt.savefig('/Users/lizekai/Downloads/20210614_codes/' + str(round(t, 2)) + '.png')
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
