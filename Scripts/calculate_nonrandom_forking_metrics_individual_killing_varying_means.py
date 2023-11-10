#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 21:01:37 2022

@author: narain

Tested in Python 3.7.4.

Calculate metrics for the forking network with individual vessel pruning and varying means.
"""

# =============================================================================
# LIBRARIES & INITIALISATION
# =============================================================================

# Initialise libraries
#import matplotlib.colors as colors
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import time

# Import tools for Paraview data
from get_paraview_data import *

# Starts stopwatch to clock execution time
start_time = time.time()

# Set LaTex-style font
from pathlib import Path
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 22})

# Set poster format
#plt.style.use('seaborn-poster')
#plt.rcParams['font.family'] = 'Arial'
#plt.rcParams.update({'font.size': 20})

# =============================================================================
# FUNCTIONS
# =============================================================================

# Define a function to read an O2 distribution .vti file and return the basic stats
def get_distribution_stats(solver_name, alpha_value, mean_value, kills, hypoxic_threshold_list, diameter_threshold_list, flow_rate_threshold_list, plot=0, read=0):    
    
    # Set the file path
    folder_path = main_folder_path + solver_name + 'Haematocrit/Lambda4/Alpha' + alpha_value + '/Mean' + mean_value + '/Kills' + kills
    # field_path =  folder_path + '/oxygen_solution_0.vti'    
    network_path = folder_path + '/FinalHaematocrit.vtp'
    
    # Print status update
    # print(field_path)

    # Save file or read from existing file
    if read==0:    
        
        # Compute the architectural metrics
        n_vessels, n_perfused_vessels, n_unperfused_vessels, mean_diameter, mean_geometric_resistance, diameter_binned, flow_rate_binned, n_connected_components, n_cycles, size_largest_connected_component, n_cycles_largest_connected_component = get_forking_predictors(network_path, reference_rank_lengths, reference_node_coordinates, pq_threshold, diameter_threshold_list, flow_rate_threshold_list)

    return 0, 0, 0, 0, 0, 0, n_vessels, n_perfused_vessels, n_unperfused_vessels, mean_diameter, mean_geometric_resistance, diameter_binned, flow_rate_binned, n_connected_components, n_cycles, size_largest_connected_component, n_cycles_largest_connected_component


# Define a function to return statistics for all the heterogeneities in the data
def get_solver_stats(solver_name, alpha_value, mean_list, kills_list, hypoxic_threshold_list, plot, read):
#    table = np.array([])
    alpha_table = np.array([])
    for mean_value in mean_list:
        for kills in kills_list:    
            average_kills_data = get_distribution_stats(solver_name, alpha_value, mean_value, kills, hypoxic_threshold_list, diameter_threshold_list, flow_rate_threshold_list, plot, read)
            table_entry = np.hstack([float(mean_value), float(kills), average_kills_data])
            alpha_table = np.vstack([alpha_table, table_entry]) if alpha_table.size else table_entry
    return alpha_table

# =============================================================================
# DISTRIBUTION STATS & HYPOXIC FRACTIONS
# =============================================================================

# Enter details to allow looping over folders
n_gens = 6
d_inlet = 75    
max_kills = 249
# max_kills = 508
# max_kills = 124
pq_threshold = 3.e-12
which_mean = 0
which_alpha = 1

# Set folder path
main_folder_path = '/home/narain/Desktop/Results/Forking/Nonrandom/Without Oxygen/Individual Pruning in Nonrandom Forking Network with Varying Means (' + str(n_gens) + ' Gens ' + str(d_inlet) + ' Inlet)/TestDichotomousNetwork/'
#main_folder_path = '/tmp/narain/testoutput/TestDichotomousNetwork/'
# main_folder_path = '/tmp/narain/testoutput/TestDichotomousNetworkHigherH/'

# Set simulation specifics
solver_list = ['Constant', 'Pries', 'Memory', 'Fung']
alpha_list = ['1.00', '1.10', '1.20', '1.30', '1.40']
graph_alpha_list = ['1.0', '1.1', '1.2', '1.3', '1.4']  # what we want displayed on the graph
mean_list = ['0', '1', '2']
graph_mean_list = ['22.76 μm', '28.50 μm', '33.65 μm']  # what we want displayed on the graph
kills_list = [str(x) for x in range(max_kills + 1)]
hypoxic_threshold_list = [2195, 10000, 15000, 20000, 25000, 27441] 
#hypoxic_threshold_list = [2195, 5488, 10976, 16465, 21953, 27441] 
#hypoxic_threshold_list = [1,2,3,4] 
diameter_threshold_list = [19, 27, 36]  # in um
flow_rate_threshold_list = [1e-12, 2e-12, 8e-12]
#hypoxic_threshold_list_pp = [0.8, 2, 4, 6, 8, 10] 

# Get reference forking network (all node coordinates are the same, regardless of heterogeneity)
generation_coordinates, reference_rank_lengths, reference_node_coordinates = get_reference_forking_network(generations=n_gens, inlet_diameter=d_inlet)

# Set solver and alpha
solver_name = solver_list[0]
alpha_value  = alpha_list[which_alpha]

# Set file name for saving plots
#description = '_nonrandom_forking_individual_pruning_geometric_metrics_with_alpha_' + alpha_value
description = 'forking_individual_pruning_alpha_' + alpha_value + '_mean_' + str(which_mean)

# Get the stats for all solvers (change to read=1 to extract from .vti files directly)
solver_stats = get_solver_stats(solver_name, alpha_value, mean_list, kills_list, hypoxic_threshold_list, plot=0, read=0)

# Save array
#np.save(main_folder_path + solver_name + 'Haematocrit/python_solver_data.npy', solver_stats)
#solver_stats = np.load(main_folder_path + solver_name + 'Haematocrit/python_solver_data.npy')

# Filter by mean
mean_array, hypoxic_fraction_composite, mean_composite, min_composite, half_composite, max_composite, sd_composite, n_vessels_composite, n_perfused_vessels_composite, n_unperfused_vessels_composite, mean_diameter_composite, mean_geometric_resistance_composite, diameter_binned_composite, flow_rate_binned_composite, n_connected_components_composite, n_cycles_composite, size_largest_connected_component_composite, n_cycles_largest_connected_component_composite = filter_by_alpha(mean_list, solver_stats)

# =============================================================================
# PERFUSION QUOTIENTS
# =============================================================================

# Read PQ file
filename = main_folder_path + 'forking_nonrandom_individual_pruning_perfusion_quotients.txt'
pq_df = pd.read_csv(filename, delim_whitespace=True, names=["network_name", "solver_name", "lambda", "alpha", "mean", "kills",  "PQ"])
#pq_df = pd.read_csv(filename, delim_whitespace=True, names=["alpha", "beta", "PQ"], skiprows=1)

# Filter PQ data for multiple solvers and alphas
solver_filter = solver_name + 'Haematocrit/'
pq_df = pq_df.loc[(pq_df["solver_name"] == solver_filter)]
pq_df = pq_df.loc[(pq_df["alpha"] == float(alpha_value))]

# Drop extra data
pq_df = pq_df.loc[(pq_df["kills"] <= max_kills)]

# Separate by alpha 
mean_grouped = pq_df.groupby('mean')
mean_0 = mean_grouped.get_group(0)
mean_1 = mean_grouped.get_group(1)
mean_2 = mean_grouped.get_group(2)
#alpha_3 = alpha_grouped.get_group(1.3)
#alpha_4 = alpha_grouped.get_group(1.4)

## Compute average of all selections for PQ
#line_0 = get_alpha_line(alpha_0, max_beta)
#line_1 = get_alpha_line(alpha_1, max_beta)
#line_2 = get_alpha_line(alpha_2, max_beta)
#line_3 = get_alpha_line(alpha_3, max_beta)
#line_4 = get_alpha_line(alpha_4, max_beta)

# Combine the PQs
pq_composite = np.vstack([mean_0['PQ'], mean_1['PQ'], mean_2['PQ']])

# =============================================================================
# FIGURE 8: PQ + NETWORK COMPOSITION (MEANS)
# =============================================================================
# '''
# Create a function to plot and save graphs for all the means for a single alpha value
def plot_a_single_mean(which_mean):
    
    # Set the figure design
    fig, axs = plt.subplots(1, 1, figsize=(6.25, 5.5))
    # fig.subplots_adjust(hspace = .75, wspace = 0.5)
    linestyles = ['solid', 'dashed', 'dotted', 'dashdot', (0,(5,10)), (0, (3, 1, 1, 1, 1, 1))]
    linecolours = ['#0072B2', '#009E73','#D55E00']
    width = 5

    # Initialise the axes
    # axs = axs.ravel()
    offset = len(mean_list)
    row_index = 0
    x = which_mean
    i = x
    
    # Rectify flow rate error by using second value of array as the initial value
    for mean_i in range(0,3):
        n_perfused_vessels_composite[mean_i,0]=n_perfused_vessels_composite[x,1]
        n_unperfused_vessels_composite[mean_i,0]=n_unperfused_vessels_composite[x,1]
    
    '''
    # Print the relative PF value at certain indices
    kills_of_interest_list = [25, 75, 125, 175, 225]
    for kills_of_interest in kills_of_interest_list:
        print(relative_pf[kills_of_interest])
    '''
    
    # =============================================================================
    # FIGURE 7
    # =============================================================================
    
    '''
    #  Plot the PF
    pq_composite[:,0]=pq_composite[:,1]  # fix flow rate error
    axs.plot(mean_array[:,1], (pq_composite[which_mean,:]-pq_composite[which_mean,1])*100/pq_composite[which_mean,1], label = r'$\overline{d}$ = ' + graph_mean_list[which_mean], c=linecolours[0], lw=width)
    axs.set_xlim(0)
    axs.set_ylabel(r'$\Delta_{\%} \mathcal{P}$') 
    axs.set_ylim(-100,100)  
    axs.set_xlim(0, max_kills)
    axs.grid()
    axs.set_xlabel('dosage (vessels pruned)') 
    row_index+=1    
    '''
    
    '''
    # Plot the number of vessels
    axs.plot(mean_array[:,1], n_vessels_composite[x], ls=linestyles[2], label='total', c=linecolours[0])
    axs.plot(mean_array[:,1], n_perfused_vessels_composite[x], ls=linestyles[1], label='perfused', c=linecolours[1], lw=width)
    axs.plot(mean_array[:,1], n_unperfused_vessels_composite[x], ls=linestyles[0], label='hypoperfused', c=linecolours[2], lw=width)
    axs.set_xlim(0, max_kills)
    axs.set_ylim(0, max_kills)
    axs.set_ylabel('composition (number of vessels)') 
    axs.set_xlabel('dosage (vessels pruned)') 
    axs.legend()
    axs.grid()
    row_index+=1
    '''
    
    # =============================================================================
    # FIGURE 9 (PART)
    # =============================================================================
        
    '''
    # Plot the PQ and network composition for a single mean and SD
    ax2 = axs.twinx()
    #axs.plot(mean_array[:,1], n_vessels_composite[x], ls='dotted', label='total')
    axs.plot(mean_array[:,1], n_perfused_vessels_composite[x], ls=linestyles[1], label='perfused', c=linecolours[1], lw=width)
    axs.plot(mean_array[:,1], n_unperfused_vessels_composite[x], ls=linestyles[0], label='hypoperfused', c=linecolours[2], lw=width)
    axs.legend()
    axs.set_xlim(0,max_kills)
    axs.set_ylim(0,max_kills)
    # axs.set_xlabel(r'dosage (vessels pruned) $({N}_P)$')    
    axs.set_xlabel('dosage (vessels pruned)')    
    axs.set_ylabel('composition (number of vessels)') 
    axs.grid()
    pq_composite[:,0]=pq_composite[:,1]  # fix flow rate error
    vessels_axis = mean_array[:,1]
    relative_pf = (pq_composite[which_mean,:]-pq_composite[which_mean,0])*100/pq_composite[which_mean,0]
    ax2.plot(vessels_axis, relative_pf, label='$\Delta_{\%} \mathcal{P}$', c=linecolours[0], lw=width)
    # ax2.plot(vessels_axis, relative_pf, label = graph_mean_list[i-offset*row_index], label='$\Delta_{\%} \mathcal{P}$', c='C0')
    # markers_on = [20, 45, 95, 118, 147, 176]
    # ax2.plot(vessels_axis[markers_on], relative_pf[markers_on], 'D')
    # ax2.legend(loc=0)
    ax2.set_ylim(-100,200)  
    ax2.set_ylabel(r'$\Delta_{\%} \mathcal{P}$')  #, c='C0')
    # axs.title.set_text(r'$\overline{d}$ = ' + graph_mean_list[which_mean])  
    fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=axs.transAxes)
    # '''
    
    # =============================================================================
    # FIGURE 11 
    # =============================================================================
    
    '''
    # Plot the mean geometric resistance for a solver
    for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
    #    axs[i].set_ylim([0,1.1])  # set PQ limits
        axs.plot(mean_array[:,1], mean_geometric_resistance_composite[i-offset*row_index]*10000, c=linecolours[i-offset*row_index], label = graph_mean_list[i-offset*row_index], lw=width)
        axs.set_xlim(0,max_kills)
        axs.set_ylim(0,17.5)
        axs.set_xlabel('dosage (vessels pruned)')    
        if i==len(mean_list*row_index):
            axs.title.set_text('${α}$ = ' + graph_alpha_list[which_alpha])  
        if which_alpha==1:
            axs.legend()
            axs.set_ylabel(r'$\overline{R}^{geom} (\times 10^{-4})\  (\mathrm{μm}^{-3})$') 
        axs.grid()
    row_index+=1
    '''
    
    # =============================================================================
    # FIGURE 12 (PART)
    # =============================================================================

    '''
    # Plot the raw PQ for a solver
    for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
    #    pq_composite[:,0]=pq_composite[:,1]  # fix flow rate error
        pq_composite[i-offset*row_index,0]=pq_composite[i-offset*row_index,1]  # fix flow rate error
        axs.plot(mean_array[:,1], pq_composite[i-offset*row_index,:], ls = linestyles[i-offset*row_index], c=linecolours[i-offset*row_index], label = graph_mean_list[i-offset*row_index], lw=width)
        axs.set_xlim(0,max_kills)
        axs.set_ylim(0,1)    
        axs.set_xlabel('dosage (vessels pruned)')    
    #    axs.legend(title='mean diameter', title_fontsize='14', fontsize='14')
        # axs.legend()
        if which_alpha==1:
            axs.legend()
            axs.set_ylabel(r'$\mathcal{P}$') 
        axs.grid()
        axs.title.set_text('${α}$ = ' + graph_alpha_list[which_alpha]) 
    '''

    # =============================================================================
    # FIGURE 13
    # =============================================================================

    '''
    # Plot the loops per vessel (for a single mean value since it doesn't make a difference)
    ax2 = axs.twinx()
    loops_per_vessel = n_cycles_composite/n_vessels_composite
    axs.plot(mean_array[:,1], loops_per_vessel[which_mean]*100, label = r'$\overline{\beta_1}$', c=linecolours[2], lw=width, ls=linestyles[1])
    axs.set_xlim(0,max_kills)
    axs.set_ylim([0,30])
    axs.set_xlabel('dosage (vessels pruned)')    
    axs.grid()    
    pq_composite[:,0]=pq_composite[:,1]  # fix flow rate error
    ax2.plot(mean_array[:,1], (pq_composite[which_mean,:]-pq_composite[which_mean,1])*100/pq_composite[which_mean,1], label='$\Delta_{\%} \mathcal{P}$', c=linecolours[0], lw=width, ls=linestyles[0])
    ax2.set_ylim(-100,100)  
    axs.title.set_text('${α}$ = ' + graph_alpha_list[which_alpha]) 
    if which_alpha==1:
        fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=axs.transAxes)
        axs.set_ylabel(r'$\overline{\beta_1}$ $(\times 10^{-2})$') 
    if which_alpha==3:
        ax2.set_ylabel(r'$\Delta_{\%} \mathcal{P}$')  #, c='C0
    '''

    # =============================================================================
    # FIGURE 14
    # =============================================================================

    '''
    # Plot the resistance per loop
    for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
        a = mean_geometric_resistance_composite[i-offset*row_index]
        b = n_cycles_composite[i-offset*row_index]/n_vessels_composite[i-offset*row_index]
        resistance_per_loop = np.divide(a, b, out=np.zeros_like(a), where=b!=0)    
        axs.plot(mean_array[:,1], resistance_per_loop*100, ls = linestyles[i-offset*row_index], c=linecolours[i-offset*row_index], label = graph_mean_list[i-offset*row_index], lw=width)
        axs.set_xlim(0,max_kills)
        axs.set_ylim(pow(10,-2),pow(10,1))
        axs.set_xlabel('dosage (vessels pruned)')    
        if which_alpha==1:
            axs.legend()
            axs.set_ylabel(r'$\overline{R}_β^{geom} (\times 10^{-2})\  (\mathrm{μm}^{-3})$') 
        axs.grid()
        axs.set_yscale('log')
        if i==len(mean_list*row_index):
            axs.title.set_text('${α}$ = ' + graph_alpha_list[which_alpha])  
    row_index+=1
    '''

    # =============================================================================
    # FIGURE 15
    # =============================================================================

    '''
    # Plot the resistance per loop
    for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
        a = mean_geometric_resistance_composite[i-offset*row_index]
        b = n_cycles_composite[i-offset*row_index]/n_vessels_composite[i-offset*row_index]
        resistance_per_loop = np.divide(a, b, out=np.zeros_like(a), where=b!=0)    
        axs.plot(mean_array[:,1], resistance_per_loop*100, ls = linestyles[i-offset*row_index], c=linecolours[i-offset*row_index], label = graph_mean_list[i-offset*row_index], lw=width)
        axs.set_xlim(0,max_kills)
        axs.set_ylim(pow(10,-2),pow(10,1))
        axs.set_xlabel('dosage (vessels pruned)')    
        if which_alpha==1:
            axs.legend()
            axs.set_ylabel(r'$\overline{R}_β^{geom} (\times 10^{-2})\  (\mathrm{μm}^{-3})$') 
        axs.grid()
        axs.set_yscale('log')
        if i==len(mean_list*row_index):
            axs.title.set_text('${α}$ = ' + graph_alpha_list[which_alpha])  
    row_index+=1
    '''
    
    # =============================================================================
    # SAVE IMAGE
    # =============================================================================

    # Show plots
    plt.show()
    
    # Save image
    description = 'forking_individual_pruning_solver_' + solver_name + '_alpha_' + alpha_value + '_mean_' + str(which_mean)
    file_path = Path('~/Desktop/Final Figures/' + description + '.tif').expanduser()
    fig.savefig(file_path, dpi=300, bbox_inches = 'tight')
    # file_path = Path('~/Desktop/Final Figures/' + description + '.svg').expanduser()
    # fig.savefig(file_path, dpi=300, bbox_inches = 'tight')


# plot_a_single_mean(0)
# plot_a_single_mean(1)
plot_a_single_mean(2)

# Prints execution time
print("\n--- Execution Time: %s seconds ---" % (time.time() - start_time))


















# =============================================================================
# UNUSED FUNCTIONALITY
# =============================================================================
# =============================================================================
# FIGURE 2
# =============================================================================

'''
# Set the figure design
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.5), tight_layout = {'pad': 2.5})
fig.subplots_adjust(hspace = .75, wspace = 0.5)
linestyles = ['solid', 'dashed', 'dotted', 'dashdot', (0,(5,10)), (0, (3, 1, 1, 1, 1, 1))]
linecolours = ['C0', 'g','r']

# Initialise the axes
axs = axs.ravel()
offset = len(mean_list)
row_index = 0
x = which_mean
i = x

# Rectify flow rate error by using second value of array as the initial value
for mean_i in range(0,3):
    n_perfused_vessels_composite[mean_i,0]=n_perfused_vessels_composite[x,1]
    n_unperfused_vessels_composite[mean_i,0]=n_unperfused_vessels_composite[x,1]
    '''
'''
# Plot the PQ for a single mean and SD
axs[row_index].plot(mean_array[:,1], (pq_composite[i-offset*row_index,:]-pq_composite[i-offset*row_index,1])*100/pq_composite[i-offset*row_index,1], label = r'$\overline{d}$ = ' + graph_mean_list[i-offset*row_index])
axs[row_index].set_xlim(0,max_kills)
axs[row_index].set_ylim(-100,100)
axs[row_index].set_xlabel('dosage (vessels pruned)')    
axs[row_index].set_ylabel(r'$\Delta_{\%} \mathcal{P}$') 
#axs[row_index].legend()
axs[row_index].grid()
#axs[row_index].title.set_text('${α}$ = ' + alpha_value)
row_index+=1
'''

# =============================================================================
# FIGURE 3: % CHANGE IN PQ
# =============================================================================

'''
# Plot the % change in PQ for a solver
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
    pq_composite[i-offset*row_index,0]=pq_composite[i-offset*row_index,1]  # fix flow rate error
    axs[row_index].plot(mean_array[:,1], (pq_composite[i-offset*row_index,:]-pq_composite[i-offset*row_index,1])*100/pq_composite[i-offset*row_index,1], ls = linestyles[i-offset*row_index], c=linecolours[i-offset*row_index], label = graph_mean_list[i-offset*row_index])
    axs[row_index].set_xlim(0,max_kills)
    axs[row_index].set_ylim(-100,200)    
    axs[row_index].set_xlabel('dosage (vessels pruned)')    
    axs[row_index].set_ylabel(r'$\Delta_{\%} \mathcal{P}$') 
#    axs[row_index].legend(title='mean diameter', title_fontsize='14', fontsize='14')
    axs[row_index].legend()
    axs[row_index].grid()
    axs[row_index].title.set_text('${α}$ = ' + graph_alpha_list[which_alpha])  
'''

# =============================================================================
# FIGURE 4: GEOMETRIC RESISTANCE
# =============================================================================

'''
# Plot the mean geometric resistance for a solver
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
#    axs[row_index].set_ylim([0,1.1])  # set PQ limits
    axs[row_index].plot(mean_array[:,1], mean_geometric_resistance_composite[i-offset*row_index]*10000, ls = linestyles[i-offset*row_index], c=linecolours[i-offset*row_index], label = graph_mean_list[i-offset*row_index], lw=width)
    axs[row_index].set_xlim(0,max_kills)
    axs[row_index].set_ylim([0,17])
    axs[row_index].set_xlabel('dosage (vessels pruned)')    
    axs[row_index].set_ylabel(r'$\overline{\mathcal{R}}^{geom}$ $(\times 10^{-4})$') 
    axs[row_index].legend()
    axs[row_index].grid()
#    axs[row_index].set_yscale('log')
    axs[row_index].title.set_text('${α}$ = ' + graph_alpha_list[which_alpha])  

row_index+=1
'''

# =============================================================================
# FIGURE 5: LOOPS + PQ
# =============================================================================
'''

# Set the figure design
fig, axs = plt.subplots(1, 1, figsize=(6.25, 5.5))
# fig.subplots_adjust(hspace = .75, wspace = 0.5)
linestyles = ['solid', 'dashed', 'dotted', 'dashdot', (0,(5,10)), (0, (3, 1, 1, 1, 1, 1))]
linecolours = ['C0', 'g','r']

# Initialise the axes
# axs = axs.ravel()
offset = len(mean_list)
row_index = 0
x = which_mean
i = x

# Rectify flow rate error by using second value of array as the initial value
for mean_i in range(0,3):
    n_perfused_vessels_composite[mean_i,0]=n_perfused_vessels_composite[x,1]
    n_unperfused_vessels_composite[mean_i,0]=n_unperfused_vessels_composite[x,1]

# Plot the loops per vessel (for a single mean value since it doesn't make a difference)
# for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
ax2 = axs.twinx()
loops_per_vessel = n_cycles_composite/n_vessels_composite
axs.plot(mean_array[:,1], loops_per_vessel[which_mean]*100, label = r'$\overline{\beta_1}$', c='C3', lw=width)
axs.set_xlim(0,max_kills)
axs.set_ylim([0,30])
axs.set_xlabel('dosage (vessels pruned)')    
axs.set_ylabel(r'$\overline{\beta_1}$ $(\times 10^{-2})$') 
axs.grid()    
pq_composite[:,0]=pq_composite[:,1]  # fix flow rate error
ax2.plot(mean_array[:,1], (pq_composite[which_mean,:]-pq_composite[which_mean,1])*100/pq_composite[which_mean,1], label='$\Delta_{\%} \mathcal{P}$', c='C0', lw=width)
ax2.set_ylim(-100,100)  
ax2.set_ylabel(r'$\Delta_{\%} \mathcal{P}$') 

#    axs.legend()which_mean
# fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=axs.transAxes)
axs.title.set_text('${α}$ = ' + graph_alpha_list[which_alpha])  
'''
# =============================================================================
# FIGURE 6: PQ + NETWORK COMPOSITION (ALPHAS)
# =============================================================================

'''
# Plot the PQ and network composition for a single mean and SD
ax2 = axs[row_index].twinx()
#axs[row_index].plot(mean_array[:,1], n_vessels_composite[x], ls='dotted', label='total')
axs[row_index].plot(mean_array[:,1], n_perfused_vessels_composite[x], ls='dashed', label='perfused', c='g')
axs[row_index].plot(mean_array[:,1], n_unperfused_vessels_composite[x], ls='dotted', label='hypoperfused', c='r')
axs[row_index].legend()
axs[row_index].set_xlim(0,max_kills)
axs[row_index].set_ylim(0,max_kills)
axs[row_index].set_xlabel('dosage (vessels pruned)')    
axs[row_index].set_ylabel('composition (number of vessels)') 
axs[row_index].grid()
pq_composite[:,0]=pq_composite[:,1]  # fix flow rate error
ax2.plot(mean_array[:,1], (pq_composite[which_mean,:]-pq_composite[which_mean,1])*100/pq_composite[which_mean,1], label = graph_mean_list[i-offset*row_index], c='C0')
ax2.set_ylim(-100,200)  
ax2.set_ylabel(r'$\Delta_{\%} \mathcal{P}$', c='C0')
axs[row_index].title.set_text('${α}$ = ' + graph_alpha_list[which_alpha])  
'''

'''
# Set the figure layout
fig, axs = plt.subplots(6, len(mean_list), figsize=(20, 20), tight_layout = {'pad': 2})
fig.subplots_adjust(hspace = .5, wspace=.25)
linestyles = ['solid', 'dashed', 'dotted', 'dashdot', (0,(5,10)), (0, (3, 1, 1, 1, 1, 1))]
linecolours = ['#1f77b4','r']
#plt.suptitle(solver_name + ' haematocrit solver in the heterogeneous hexagonal vessel network with stochastic pruning')

# Plot the distribution stats for a solver
axs = axs.ravel()
offset = len(mean_list)
row_index = 0

# Plot the number of vessels
for i in range(len(mean_list)):
    axs[i].plot(mean_array[:,1], n_vessels_composite[i], label='total')
    axs[i].plot(mean_array[:,1], n_perfused_vessels_composite[i], ls='dashed', label='perfused')
    axs[i].plot(mean_array[:,1], n_unperfused_vessels_composite[i], ls='dotted', label='unperfused')
    axs[i].set_xlim(0)
    axs[i].set_ylim(0,300)
#    axs[i].ticklabel_format(axis="y", scilimits=(0,0))
#    axs[i].set_xlabel('radius threshold (μm)')    
    if i==0:
        axs[i].set_ylabel('${N}_V$') 
        axs[i].legend()
    axs[i].grid()
    axs[i].title.set_text(r'$\overline{d}$ = ' + graph_mean_list[i])
row_index+=1

# Plot the PQ for a solver
for i in range(len(mean_list),len(mean_list)*(row_index+1)):
    axs[i].set_ylim([0,1.1])  # set PQ limits
    axs[i].plot(mean_array[:,1], pq_composite[i-offset,:])
    axs[i].set_xlim(0)
#    axs[i].set_xlabel('radius threshold (μm)')    
    if i==len(mean_list):
        axs[i].set_ylabel('PQ') 
#    axs[i].legend()
    axs[i].grid()
row_index+=1

# Plot the mean diameter for a solver
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
#    axs[i].set_ylim([0,1.1])  # set PQ limits
    axs[i].plot(mean_array[:,1], mean_diameter_composite[i-offset*row_index])
    axs[i].set_xlim(0)
    axs[i].set_ylim([0,100])
#    axs[i].set_xlabel('radius threshold (μm)')    
    if i==len(mean_list*row_index):
        axs[i].set_ylabel('$\overline{D}$ (μm)') 
#    axs[i].legend()
    axs[i].grid()
row_index+=1

# Plot the mean geometric resistance for a solver
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
#    axs[i].set_ylim([0,1.1])  # set PQ limits
    axs[i].plot(mean_array[:,1], mean_geometric_resistance_composite[i-offset*row_index]*10000)
    axs[i].set_xlim(0)
#    axs[i].set_ylim([0,1])
    axs[i].set_xlabel('kills (vessels)')    
    if i==len(mean_list*row_index):
        axs[i].set_ylabel(r'$\overline{R}^{geom} (\times 10^{-4})$') 
#    axs[i].legend()
    axs[i].grid()
row_index+=1

'''
'''
# Plot the diameter composition
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
    diameter_bin_labels = ['< '+str(diameter_threshold_list[0])+' μm', str(diameter_threshold_list[0])+' μm – '+str(diameter_threshold_list[1])+' μm', str(diameter_threshold_list[1])+' μm – '+str(diameter_threshold_list[2])+' μm', '> '+str(diameter_threshold_list[2])+' μm']
    highest_bracket = np.array(n_vessels_composite[i-offset*(row_index+1)] - [diameter_bins[0] for diameter_bins in diameter_binned_composite[i-offset*(row_index+1),:]] - [diameter_bins[1] for diameter_bins in diameter_binned_composite[i-offset*(row_index+1),:]] - [diameter_bins[2] for diameter_bins in diameter_binned_composite[i-offset*(row_index+1),:]], dtype=int)
    axs[i].stackplot(np.array(mean_array[:,1], dtype=int), np.array([diameter_bins[0] for diameter_bins in diameter_binned_composite[i-offset*row_index,:]], dtype=int), np.array([diameter_bins[1] for diameter_bins in diameter_binned_composite[i-offset*row_index,:]], dtype=int), np.array([diameter_bins[2] for diameter_bins in diameter_binned_composite[i-offset*row_index,:]], dtype=int), highest_bracket, labels=[diameter_bin_labels[0],diameter_bin_labels[1],diameter_bin_labels[2],diameter_bin_labels[3]])
#    axs[i].plot(mean_array[:,1], [diameter_bins[0] for diameter_bins in diameter_binned_composite[i-offset*row_index,:]], label=str(diameter_threshold_list[0])+' μm')
#    axs[i].plot(mean_array[:,1], [diameter_bins[1] for diameter_bins in diameter_binned_composite[i-offset*row_index,:]], ls='dashed', label=str(diameter_threshold_list[1])+' μm')
#    axs[i].plot(mean_array[:,1], [diameter_bins[2] for diameter_bins in diameter_binned_composite[i-offset*row_index,:]], ls='dashed', label=str(diameter_threshold_list[2])+' μm')
    axs[i].set_xlim(0)
    axs[i].set_ylim([0,300])
    axs[i].set_xlabel('kills (vessels)')    
    if i==len(mean_list*row_index):
        axs[i].set_ylabel('${N}_V$') 
        handles, labels = axs[i].get_legend_handles_labels()
        axs[i].legend(reversed(handles), reversed(labels), title='Diameter')    
    axs[i].grid()
row_index+=1
    
# Plot the flow rate composition
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
    flow_rate_bin_labels = ['< '+str(flow_rate_threshold_list[0])+' $m^{3}$/s', str(flow_rate_threshold_list[0])+' $m^{3}$/s – '+str(flow_rate_threshold_list[1])+' $m^{3}$/s', str(flow_rate_threshold_list[1])+' $m^{3}$/s – '+str(flow_rate_threshold_list[2])+' $m^{3}$/s', '> '+str(flow_rate_threshold_list[2])+' $m^{3}$/s']
    highest_bracket = np.array(n_vessels_composite[i-offset*row_index] - [flow_rate_bins[0] for flow_rate_bins in flow_rate_binned_composite[i-offset*row_index,:]] - [flow_rate_bins[1] for flow_rate_bins in flow_rate_binned_composite[i-offset*row_index,:]] - [flow_rate_bins[2] for flow_rate_bins in flow_rate_binned_composite[i-offset*row_index,:]], dtype=int)
    axs[i].stackplot(np.array(mean_array[:,1], dtype=int), np.array([flow_rate_bins[0] for flow_rate_bins in flow_rate_binned_composite[i-offset*row_index,:]], dtype=int), np.array([flow_rate_bins[1] for flow_rate_bins in flow_rate_binned_composite[i-offset*row_index,:]], dtype=int), np.array([flow_rate_bins[2] for flow_rate_bins in flow_rate_binned_composite[i-offset*row_index,:]], dtype=int), highest_bracket, labels=[flow_rate_bin_labels[0],flow_rate_bin_labels[1],flow_rate_bin_labels[2],flow_rate_bin_labels[3]])
#    axs[i].plot(mean_array[:,1], [flow_rate_bins[0] for flow_rate_bins in flow_rate_binned_composite[i-offset*row_index,:]], label=str(flow_rate_threshold_list[0])+' $m^{3}$/s')
#    axs[i].plot(mean_array[:,1], [flow_rate_bins[1] for flow_rate_bins in flow_rate_binned_composite[i-offset*row_index,:]], ls='dashed', label=str(flow_rate_threshold_list[1])+' $m^{3}$/s')
#    axs[i].plot(mean_array[:,1], [flow_rate_bins[2] for flow_rate_bins in flow_rate_binned_composite[i-offset*row_index,:]], ls='dashed', label=str(flow_rate_threshold_list[2])+' $m^{3}$/s')
    axs[i].set_xlim(0)
    axs[i].set_ylim([0,300])
    axs[i].set_xlabel('kills (vessels)')    
    if i==len(mean_list*row_index):
        axs[i].set_ylabel('${N}_V$') 
        handles, labels = axs[i].get_legend_handles_labels()
        axs[i].legend(reversed(handles), reversed(labels), title='Flow Rate')    
    axs[i].grid()
row_index+=1

# Plot the number of connected components
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
    axs[i].plot(mean_array[:,1], n_connected_components_composite[i-offset*row_index])
    axs[i].set_xlim(0)
#    axs[i].set_ylim([0,300])
    axs[i].set_xlabel('kills (vessels)')    
    if i==len(mean_list*row_index):
        axs[i].set_ylabel('$β_0$') 
    axs[i].grid()
row_index+=1
    
# Plot the number of loops
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
    axs[i].plot(mean_array[:,1], n_cycles_composite[i-offset*row_index])
    axs[i].set_xlim(0)
#    axs[i].set_ylim([0,300])
    axs[i].set_xlabel('kills (vessels)')    
    if i==len(mean_list*row_index):
        axs[i].set_ylabel('$β_1$') 
    axs[i].grid()
row_index+=1
'''
'''
# Plot the loops per vessel
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
    loops_per_vessel = n_cycles_composite[i-offset*row_index]/n_vessels_composite[i-offset*row_index]
    axs[i].plot(mean_array[:,1], loops_per_vessel*100)
    axs[i].set_xlim(0)
#    axs[i].set_ylim([0,.1])
    axs[i].set_xlabel('kills (vessels)')    
    if i==len(mean_list*row_index):
        axs[i].set_ylabel(r'$\overline{β_1} (\times 10^{-2})$') 
    axs[i].grid()
row_index+=1

# Plot the resistance per loop
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
    
    a = mean_geometric_resistance_composite[i-offset*row_index]
    b = n_cycles_composite[i-offset*row_index]/n_vessels_composite[i-offset*row_index]
    resistance_per_loop = np.divide(a, b, out=np.zeros_like(a), where=b!=0)
    axs[i].plot(mean_array[:,1], resistance_per_loop*100)
    axs[i].set_xlim(0)
#    axs[i].set_ylim([0,.1])
    axs[i].set_xlabel('kills (vessels)')    
    if i==len(mean_list*row_index):
        axs[i].set_ylabel(r'$\overline{R}_β^{geom} (\times 10^{-2})$') 
    axs[i].grid()

'''

# =============================================================================
# POSTER PLOTS
# =============================================================================

'''
mean_list = ['0']#, '1', '2']

# Set the figure layout
fig, axs = plt.subplots(3, 1, figsize=(7.5, 12), tight_layout = {'pad': 2})
fig.subplots_adjust(hspace = .5, wspace=.25)
linestyles = ['dotted', 'dashdot', 'solid', 'dashed', (0,(5,10)), (0, (3, 1, 1, 1, 1, 1))]
colors = ['royalblue', 'crimson','forestgreen']
#plt.suptitle(solver_name + ' haematocrit solver in the heterogeneous hexagonal vessel network with stochastic pruning')
line_changes = np.diff(n_perfused_vessels_composite)

# Plot the distribution stats for a solver
axs = axs.ravel()
offset = len(mean_list)
row_index = 0

# Plot the loops per vessel
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
    loops_per_vessel = n_cycles_composite/n_vessels_composite
    axs[row_index].plot(mean_array[:,1], loops_per_vessel*100, label = r'$\overline{d}$ = ' + graph_mean_list[i-offset*row_index], c=colors[0])
    axs[row_index].set_xlim(0)
    axs[row_index].axvline(157, c='grey', alpha=0.5)  # highlight loops per size
#    axs[row_index].set_ylim([0,.1])
    axs[row_index].set_xlabel('vessels killed')    
    axs[row_index].set_ylabel(r'loops per vessel $(\times 10^{-2})$') 
#    axs[row_index].grid()    
    axs[row_index].spines['right'].set_visible(False)
    axs[row_index].spines['top'].set_visible(False)
row_index+=1

#
##plt.axvspan(3, 6, color='red', alpha=0.5)
#plt.axvspan(3, 6, color='red', alpha=0.5)
#


# Show plots
plt.show()
'''