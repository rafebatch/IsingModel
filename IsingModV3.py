import numpy as np
import math

import matplotlib as mpl       # As of July 2017 Bucknell computers use v. 2.x
import matplotlib.pyplot as plt
from scipy.ndimage import correlate
from scipy import stats

## define lattice size
lattice_size = 100

## create variable to store number of spins
L_sq = lattice_size**2


def flip(lattice):
    ## create matrix of random indices for row and col
    ran_num = stats.randint.rvs(0, lattice_size, size=2 * L_sq)

    ## calculate array of random probabilities to determine flip or not when energy gain
    prob = stats.uniform.rvs(0, 1, size=L_sq)

    for i in range(L_sq):
        ## set row, col using random indices in array
        col_r = ran_num[i]
        row_r = ran_num[i + L_sq]

        ## calculate initial energy (sum of neighbors)
        e_init = calc_E(lattice, row_r, col_r)

        ## final energy is simply the flipped spin * sum of neigbhors

        deltE = 2 * e_init * lattice[row_r, col_r]
        lattice[row_r, col_r] *= -1

        ## use negative deltE since negative means energy gain
        if (deltE) > 0:
            if deltE == 4 and prob[i] > p1:
                lattice[row_r, col_r] *= -1
            if deltE == 8 and prob[i] > p2:
                lattice[row_r, col_r] *= -1

def calc_E(lattice, row, col):
    ## sum of neighbors
    return lattice[row-1,col] + lattice[(row+1)%lattice_size, col] + lattice[row, col-1] + lattice[row, (col+1)%lattice_size]

def avg_mag(lattice):
    return (np.sum(lattice) / L_sq)    ## avg magnetization per spin


def energy_sys(lattice):
    ## create shifted matrices
    up = np.roll(lattice, 1, axis=0)
    left = np.roll(lattice, -1, axis=1)

    ## calculate multiplication between
    tot_E_up = lattice * up
    tot_E_left = lattice * left

    tot_E = np.sum(tot_E_up + tot_E_left)

    return -tot_E / lattice_size ** 2


## initialize fresh aligned lattice
lattice = np.ones((lattice_size, lattice_size))

## initialize trial count and temp range
num_trials = 500
temp_range = np.linspace(0.1, 4, 400)
temp_range = np.round(temp_range, 2)
avg_num = 100

## create arrays to store avg mag and avg energy
magn = np.zeros(num_trials)
ener = np.zeros(num_trials)

## open/prepare file to store equilibrium results
file_equi = open("EquiLib.dat", "w")
for i in temp_range:
    ## initialize fresh lattice for each temp
    lattice = np.ones((lattice_size, lattice_size))

    ## set temperature
    T = i
    T = round(T, 1)

    ## calculate probabilities to be used
    p1 = np.exp(-4 / T)
    p2 = np.exp(-8 / T)

    ## open file and write first line
    file = open("IM_T" + str(T) + ".dat", "w")
    file.write("Trial\tMag\t\tEnergy\t\tTEMP=" + str(T) + "\n")
    if T >= 2:
        num_trials = 5000
        avg_num = 2000
        magn = np.zeros(num_trials)
        ener = np.zeros(num_trials)
    if T >= 3:
        num_trials = 500
        avg_num = 100
        magn = np.zeros(num_trials)
        ener = np.zeros(num_trials)

    for j in range(num_trials):
        ## call flip and calculate magnetization after 100 flips
        flip(lattice)
        magn[j] = avg_mag(lattice)
        ener[j] = energy_sys(lattice)

        ## write to file -- still need to calculate energy
        file.write(str(j) + "\t\t" + str(magn[j]) + "\t\t" + str(ener[j]) + "\n")

    ## take average magnetization and energy at end of each trial to get equilib values
    equi_mag = np.sum(magn[-avg_num:]) / avg_num
    equi_energy = np.sum(ener[-avg_num:]) / avg_num
    file_equi.write(str(T) + "\t" + '%.6f' % (equi_mag) + "\t\t" + '%.6f' % (equi_energy) + "\n")

    file.close()

file_equi.close()