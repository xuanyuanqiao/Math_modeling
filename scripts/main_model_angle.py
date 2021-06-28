import math
import matplotlib.pyplot as plt
import random
import numpy as np
from matplotlib import colors


def sigma(x):
    return (1 + math.tanh(2 * x)) / 2


def included_angle(v1, v2):
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

# compute the signal received by the cell
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
        angle_list[cell_index] = included_angle([surrounding_cell_1[cell_index].xposition - target.xposition,
                                                 surrounding_cell_1[cell_index].yposition - target.yposition], [
                                                    surrounding_cell_1[(cell_index + 1) % len(
                                                        surrounding_cell_1)].xposition - target.xposition,
                                                    surrounding_cell_1[(cell_index + 1) % len(
                                                        surrounding_cell_1)].yposition - target.yposition])
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

# cell class
class Cell:
    def __init__(self, startstate, xposition, yposition, xindex, yindex):
        self.state = startstate
        self.xposition = xposition
        self.yposition = yposition
        self.xindex = xindex
        self.yindex = yindex
        self.signal_receive = 0

    # Signal received by the cell
    def receive(self, t, neighboringinput):
        x = self.xposition
        signal = signaling_gradient(x, t) + neighboringinput
        self.signal_receive = signal

    # Change the state of cell
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

# parameters
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

# initialization
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
    for row in back:
        for cell in row:
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

    for row in back:
        for cell in row:
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
        plt.savefig('./figure/improvement_angle/' + str(round(t, 2)) + '.png')
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

