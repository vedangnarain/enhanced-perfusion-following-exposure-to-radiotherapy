#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 15:31:25 2022

@author: narain
  
Tested in Python 3.7.4.

Plots the biological networks.
"""

# =============================================================================
# LIBRARIES & INITIALISATION
# =============================================================================

# Initialise libraries
#import matplotlib.colors as colors
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import time
from pathlib import Path
from matplotlib.patches import Patch

# Starts stopwatch to clock execution time
start_time = time.time()

# Set LaTex-style font
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 22})

# Set poster format
# plt.style.use('seaborn-poster')
# plt.rcParams['font.family'] = 'Arial'
# plt.rcParams.update({'font.size': 25})

# =============================================================================
# DATA
# =============================================================================

tumor_18_2D = np.array([[720, 522, 446],[0.498, 0.076, 0.449]])
tumor_18_3A = np.array([[1000, 709],[0.143,	0.408]])
tumor_24_2B = np.array([[888, 603, 548],[0.38, 0.18, 0.18]])
tumor_33_2B = np.array([[4087, 1374, 1366],[0.26,0.56,0.55]]) 
tumor_35_2C = np.array([[1038, 982],[0.36, 0.55]])
tumor_47_1I = np.array([[2175, 990, 784],[0.615, 0.447, 0.475]])

# =============================================================================
# PLOTS
# =============================================================================
'''
# Set the figure layout
fig, axs = plt.subplots(2, 1, figsize=(10, 10), tight_layout = {'pad': 2})
fig.subplots_adjust(hspace = .5, wspace=.25)
linestyles = ['solid', 'dashed', 'dotted', 'dashdot', (0, (3, 1, 1, 1, 1, 1)), (0,(9,10))]
linecolours = ['g','r']
#plt.suptitle(solver_name + ' haematocrit solver in the heterogeneous hexagonal vessel network with stochastic pruning')

# Plot the distribution stats for a solver
axs = axs.ravel()

axs[0].plot(abs(tumor_18_2D[0] - tumor_18_2D[0,0]), tumor_18_2D[1], 'o', label='18_2D', ls=linestyles[0], c=linecolours[1])
axs[0].plot(abs(tumor_18_3A[0] - tumor_18_3A[0,0]), tumor_18_3A[1], 'o', label='18_3A', ls=linestyles[1], c=linecolours[0])
axs[0].plot(abs(tumor_24_2B[0] - tumor_24_2B[0,0]), tumor_24_2B[1], 'o', label='24_2B', ls=linestyles[2], c=linecolours[1])
axs[0].plot(abs(tumor_33_2B[0] - tumor_33_2B[0,0]), tumor_33_2B[1], 'o', label='33_2B', ls=linestyles[3], c=linecolours[0])
axs[0].plot(abs(tumor_35_2C[0] - tumor_35_2C[0,0]), tumor_35_2C[1], 'o', label='35_2C', ls=linestyles[4], c=linecolours[0])
axs[0].plot(abs(tumor_47_1I[0] - tumor_47_1I[0,0]), tumor_47_1I[1], 'o', label='47_1I', ls=linestyles[5], c=linecolours[1])
axs[0].set_ylabel('PQ') 
axs[0].set_xlabel('kills (vessels)') 
axs[0].legend()
axs[0].grid()

axs[1].plot(abs(tumor_18_2D[0] - tumor_18_2D[0,0])*100/tumor_18_2D[0,0], tumor_18_2D[1], 'o', label='18_2D', ls=linestyles[0], c=linecolours[1])
axs[1].plot(abs(tumor_18_3A[0] - tumor_18_3A[0,0])*100/tumor_18_3A[0,0], tumor_18_3A[1], 'o', label='18_3A', ls=linestyles[1], c=linecolours[0])
axs[1].plot(abs(tumor_24_2B[0] - tumor_24_2B[0,0])*100/tumor_24_2B[0,0], tumor_24_2B[1], 'o', label='24_2B', ls=linestyles[2], c=linecolours[1])
axs[1].plot(abs(tumor_33_2B[0] - tumor_33_2B[0,0])*100/tumor_33_2B[0,0], tumor_33_2B[1], 'o', label='33_2B', ls=linestyles[3], c=linecolours[0])
axs[1].plot(abs(tumor_35_2C[0] - tumor_35_2C[0,0])*100/tumor_35_2C[0,0], tumor_35_2C[1], 'o', label='35_2C', ls=linestyles[4], c=linecolours[0])
axs[1].plot(abs(tumor_47_1I[0] - tumor_47_1I[0,0])*100/tumor_47_1I[0,0], tumor_47_1I[1], 'o', label='47_1I', ls=linestyles[5], c=linecolours[1])
axs[1].set_ylabel('PQ') 
axs[1].set_xlabel('kills (% vessels)') 
#axs[1].legend()
axs[1].grid()
    
# Show plots
plt.show()
'''
# =============================================================================
# DATA
# =============================================================================

tumour_labels = ['18_3A', '33_2B', '35_2C', '47_1I', '24_2B', '18_2D']
tumour_pq_change = [185.31, 115.38, 52.78, -22.76, -37.93, -53.85]
tumour_mean_diameters = [22.75, 24.63, 26.1, 32.35, 33.65, 31.51]
tumour_loops_per_size = [2.2, 6.04, 4.62, 6.16, 7.45, 6.82]
colors = ['royalblue','royalblue','royalblue','crimson','crimson','crimson']
legend_elements = [Patch(facecolor='royalblue', label='well-responding tumour'),
                   Patch(facecolor='crimson', label='poorly-responding tumour')]

# Create the figure
plt.show()
# Set the figure layout
fig, axs = plt.subplots(3, 1, figsize=(10, 16), tight_layout = {'pad': 2})
fig.subplots_adjust(hspace = .5, wspace=.25)
#plt.suptitle(solver_name + ' haematocrit solver in the heterogeneous hexagonal vessel network with stochastic pruning')

# Plot the distribution stats for a solver
axs = axs.ravel()
#axs.spines['top'].set_visible(False)
axs[0].bar(tumour_labels, tumour_pq_change, color=colors, label='well-responding')
axs[0].set_ylabel(r'% change in $\mathcal{P}$') 
axs[0].set_xlabel('tumours') 
axs[0].spines['right'].set_visible(False)
axs[0].spines['top'].set_visible(False)
axs[0].legend(handles=legend_elements, loc='best')
#    axs[row_index].set_xlabel('radius threshold (μm)')    

axs[1].bar(tumour_labels, tumour_mean_diameters, color=colors)
axs[1].set_ylabel('mean diameter (μm)') 
axs[1].set_xlabel('tumours') 
axs[1].spines['right'].set_visible(False)
axs[1].spines['top'].set_visible(False)
axs[1].legend(handles=legend_elements, loc='best')

axs[2].bar(tumour_labels, tumour_loops_per_size, color=colors)
axs[2].set_ylabel(r'loops per vessel $(\times 10^{-2})$') 
axs[2].set_xlabel('tumours') 
axs[2].spines['right'].set_visible(False)
axs[2].spines['top'].set_visible(False)
axs[2].legend(handles=legend_elements, loc='best')

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

# Save image
figure_title = 'summary'
file_path = Path('~/Desktop/Final Figures/' + figure_title + '.svg').expanduser()
fig.savefig(file_path, dpi=500, bbox_inches = 'tight')
file_path = Path('~/Desktop/Final Figures/' + figure_title + '.png').expanduser()
fig.savefig(file_path, dpi=500, bbox_inches = 'tight')