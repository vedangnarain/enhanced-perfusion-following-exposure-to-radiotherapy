#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 16:47:21 2022

@author: narain

Tested in Python 3.7.4.

This script calculates the changes required to the average diameter of a 
synthetic network to have it match the range of the biological networks.
"""

# Initialise libraries
#import errno
import matplotlib.pyplot as plt
import numpy as np
#import os, os.path
#import scipy.stats as stats

# Import tools for Paraview data
from get_paraview_data import *

# Set LaTex-style font
#from pathlib import Path
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 11})

# =============================================================================
# FUNCTIONS
# =============================================================================

# Define a function to calculate the offset to a network to match the biological diameter range
def calculate_required_offset(reference_network_path):
    
    # Get the synthetic network's diameter data
    segment_diameters_um, mean_diameter, std_dev = get_diameter_distribution(reference_network_path)
#    print(min(segment_diameters_um))
    
    # Calculate the offsets
    min_offset = min_bio_diameter - mean_diameter
    mean_offset = mean_bio_diameter - mean_diameter
    max_offset = max_bio_diameter - mean_diameter
    
    # Print the offsets    
    print('Offsets: ', round(min_offset,2), round(mean_offset,2), round(max_offset,2))
    
    # Print the diameter of the smallest vessel
#    print('Min Diameter: ', min(segment_diameters_um)+min_offset)
    
#    print(np.diff(np.unique(segment_diameters_um)).min())

# =============================================================================
# OFFSETS
# =============================================================================

# Enter the inlet diameter and number of generations
d_inlet = 75
n_gen = 6

# Enter the biological diameters
biological_network_mean_diameters = [22.76,24.62,26.1,31.52,32.36,33.64]

# Calculate the minimum, average, and maximum biological diameter averages
min_bio_diameter = np.min(biological_network_mean_diameters)
mean_bio_diameter = np.mean(biological_network_mean_diameters)
max_bio_diameter = np.max(biological_network_mean_diameters)

parent_network_folder = '/home/narain/Desktop/Scripts/reference_networks/forking_networks/inlet_diameter_' + str(d_inlet) +'_um_generations_' + str(n_gen) + '/Alpha'

# Calculate the required offsets
calculate_required_offset(parent_network_folder+'1.00/FinalHaematocrit.vtp')
calculate_required_offset(parent_network_folder+'1.10/FinalHaematocrit.vtp')
calculate_required_offset(parent_network_folder+'1.20/FinalHaematocrit.vtp')
calculate_required_offset(parent_network_folder+'1.30/FinalHaematocrit.vtp')
calculate_required_offset(parent_network_folder+'1.40/FinalHaematocrit.vtp')

# =============================================================================
# UNUSED FUNCTIONALITY
# =============================================================================

'''
d_P = 100
mean_1 = 8
mean_2 = 18

d_A = np.arange(1,101)

# Obeys parent vessel diameter
plt.plot(d_A,d_P/d_A,label='dP='+str(d_P))

# Obeys the overall mean
#d_B = (2*mean_1)-d_A
plt.plot(d_A,(2*mean_1)-d_A,label='u='+str(mean_1))
plt.plot(d_A,(2*mean_2)-d_A,label='u='+str(mean_2))

# Obeys the standard deviation

plt.xlabel('dA')
plt.ylabel('dB')
plt.legend()
'''

# Plot the effect of diameter on resistance
'''
plt.figure()
x = np.arange(3)
alpha_1 = [0.00167483648429,0.000411644303887,0.000168230172615]
alpha_2 = [0.04283926010914,0.001102676806587,0.00030192085884]
alpha_3 = [1.15873515177815, 0.001757657727125, 0.000393424690973]
width = 0.20
  
# plot data in grouped manner of bar type
plt.bar(x-0.2, alpha_1, width, color='cyan')
plt.bar(x, alpha_2, width)
plt.bar(x+0.2, alpha_3, width)
plt.xticks(x, ['22.75', '28.5', '33.65'])
plt.xlabel("mean diameter (μm)")
plt.ylabel(r'$\overline{R}^{geom}$')
plt.legend(["α = 1.1", "α = 1.3", "α = 1.4"])
plt.yscale('log')
plt.show()
'''

'''
# Set the figure layout
fig, axs = plt.subplots(3, 1, figsize=(7.5, 16), tight_layout = {'pad': 2})
fig.subplots_adjust(hspace = .5, wspace=.25)
#plt.suptitle(solver_name + ' haematocrit solver in the heterogeneous hexagonal vessel network with stochastic pruning')

# Plot the distribution stats for a solver
axs = axs.ravel()
#axs.spines['top'].set_visible(False)
axs[0].bar(tumour_labels, tumour_pq_change, color=colors, label='well-responding')
axs[0].set_ylabel('change in PQ (%)') 
axs[0].spines['right'].set_visible(False)
axs[0].spines['top'].set_visible(False)
axs[0].legend(handles=legend_elements, loc='best')
#    axs[row_index].set_xlabel('radius threshold (μm)')    

axs[1].bar(tumour_labels, tumour_mean_diameters, color=colors)
axs[1].set_ylabel('mean diameter (μm)') 
axs[1].spines['right'].set_visible(False)
axs[1].spines['top'].set_visible(False)

axs[2].bar(tumour_labels, tumour_loops_per_size, color=colors)
axs[2].set_ylabel(r'loops per vessel $(\times 10^{-2})$') 
axs[2].set_xlabel('tumours') 
axs[2].spines['right'].set_visible(False)
axs[2].spines['top'].set_visible(False)
## Plot the PQ for a solver
#for i in range(0,1:
#    axs[row_index].set_ylim([0,10])  # set PQ limits
##    axs[row_index].set_ylim([0,1.1])  # set PQ limits
#    axs[row_index].plot(mean_array[:,1], (pq_composite[i-offset*row_index,:]-pq_composite[i-offset*row_index,2])*100/pq_composite[i-offset*row_index,2], ls = linestyles[i-offset*row_index], label = '${µ}$ = ' + mean_list[i-offset*row_index])
##    axs[row_index].plot(mean_array[:,1], pq_composite[i-offset*row_index,:], ls = linestyles[i-offset*row_index], label = '${µ}$ = ' + mean_list[i-offset*row_index])
##    axs[row_index].fill_between(mean_array[:,1], pq_composite[:,(3*i)-10]+sd_pq_composite[:, (3*i)-10], pq_composite[:,(3*i)-10]-sd_pq_composite[:, (3*i)-10], color='grey', alpha=0.5, label='mean ± SD')
#    axs[row_index].set_xlim(0)
##    axs[row_index].set_xlabel('radius threshold (μm)')    
##    axs[row_index].set_ylabel('PQ') 
#    axs[row_index].set_ylabel('Relative % PQ') 
#    axs[row_index].grid(True)
#    axs[row_index].title.set_text('SD = ' + sd_folder)
#    axs[row_index].legend()
#row_index+=1

# Show plots
plt.show()
'''

'''
# Get the .vtk data
point_coordinates_data_array, point_data_array, cell_data_array, polydata = get_vtk_data(parent_network_folder+'1.00/FinalHaematocrit.vtp')
    
# Convert the radii in metres into the right units
segment_diameters_um = cell_data_array['Vessel Radius m']*2*1000000  # convert radii (m) to diameters (um)
mean_diameter = np.mean(segment_diameters_um)
std_dev = np.std(segment_diameters_um)
    
# Set the figure layout
fig, axs = plt.subplots(nrows=1, ncols=1, sharex=True)
fig.subplots_adjust(hspace = 0.75, wspace=.25)

# Calculate the offsets
min_offset = min_bio_diameter - mean_diameter
mean_offset = mean_bio_diameter - mean_diameter
max_offset = max_bio_diameter - mean_diameter

# Plot the distribution
n, bins, patches = axs.hist(x=segment_diameters_um, bins='auto', alpha=0.9, label='natural', color='black')
n, bins, patches = axs.hist(x=segment_diameters_um+min_offset, bins='auto', alpha=0.5, label='${µ}$ = ' + str(round(np.mean(segment_diameters_um+min_offset),2)) + ' μm')
n, bins, patches = axs.hist(x=segment_diameters_um+mean_offset, bins='auto', alpha=0.5, label='${µ}$ = ' + str(round(np.mean(segment_diameters_um+mean_offset),2)) + ' μm')
n, bins, patches = axs.hist(x=segment_diameters_um+max_offset, bins='auto', alpha=0.5, label='${µ}$ = ' + str(round(np.mean(segment_diameters_um+max_offset),2)) + ' μm')
axs.set_ylabel('number of vessels') 
axs.set_xlabel('vessel diameter (μm)')    
axs.title.set_text('${σ}$ = ' + str(round(std_dev,2)) + ' μm')
axs.tick_params(labelbottom=True)
#axs.axvline(mean_diameter, c='black', ls='--', label='${µ}$ = ' + str(round(mean_diameter,2)) + ' μm')
axs.set_ylim(0, 75)
axs.legend()
'''
