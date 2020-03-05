import numpy as np
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.ndimage import correlate
from scipy import stats

file = open("EquiLib.txt", 'r')
num_lines = sum(1 for line in open("EquiLib.txt"))

## read first line -- header
file.readline()
count = 0

temp = np.zeros(num_lines - 1)
mag = np.zeros(num_lines - 1)
energy = np.zeros(num_lines - 1)

for line in file:
    ## get contents of line
    string = line.split()

    ## store contents in array
    temp[count] = string[0]
    mag[count] = string[1]
    energy[count] = string[2]

    ## increase count
    count += 1

plt.figure()
plt.scatter(temp, mag)
plt.xlabel("Temperature")
plt.ylabel("Average magnetization per spin")
plt.title("Temperature vs. Average Equilibrium Magnetization")

plt.figure()
plt.scatter(temp, energy)
plt.xlabel("Temperature")
plt.ylabel("Average energy per spin")
plt.title("Temperature vs. Average Equilibrium Energy")