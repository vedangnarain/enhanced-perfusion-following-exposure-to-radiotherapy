#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 18:39:03 2022

@author: narain

Tested in Python 3.7.4.

Calculate metrics of nonrandom forking network with threshold vessel pruning and varying means.
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
#matplotlib.rcParams['mathtext.fontset'] = 'stix'
#matplotlib.rcParams['font.family'] = 'STIXGeneral'
#matplotlib.rcParams.update({'font.size': 22})

# Set poster format
plt.style.use('seaborn-poster')
plt.rcParams['font.family'] = 'Arial'
plt.rcParams.update({'font.size': 20})

# =============================================================================
# FUNCTIONS
# =============================================================================

# Define a function to read an O2 distribution .vti file and return the basic stats
def get_distribution_stats(solver_name, alpha_value, mean_value, radius_threshold, hypoxic_threshold_list, diameter_threshold_list, flow_rate_threshold_list, plot=0, read=0):    
    
    # Set the file path
    folder_path = main_folder_path + solver_name + 'Haematocrit/Lambda4/Alpha' + alpha_value + '/Mean' + mean_value + '/RadiusThreshold' + radius_threshold
    field_path =  folder_path + '/oxygen_solution_0.vti'    
    network_path = folder_path + '/FinalHaematocrit.vtp'
    
    # Print status update
    print(field_path)

    # Save file or read from existing file
    if read==0:    
        
        # Compute the architectural metrics
        n_vessels, n_perfused_vessels, n_unperfused_vessels, mean_diameter, mean_geometric_resistance, diameter_binned, flow_rate_binned, n_connected_components, n_cycles, size_largest_connected_component, n_cycles_largest_connected_component = get_forking_predictors(network_path, reference_rank_lengths, reference_node_coordinates, pq_threshold, diameter_threshold_list, flow_rate_threshold_list)
        '''
        # Compute the O2 distribution
        field_data, field_spacing = get_vti_data(field_path)  # import the O2 distribution for the network
        middle_field = get_forking_domain(field_data, 6, generation_coordinates)  # get the O2 distribution in the middle of the field for summary statistic and HF
        
        # Plot the O2 distribution if desired
        if plot==1:  
            O2_field, _ = get_plottable_field(middle_field)  # get the O2 mesh for plots and stats
            fig = plt.figure()
            ax = plt.axes()
            colour_map = plt.cm.get_cmap('jet')
            plt.suptitle('O$_2$ distribution generated by the ' + solver_name + ' solver in the forking vessel network with with α = ' + alpha_value + 'and radius threshold = ' + radius_threshold + ' segments (1 unit = ' + str(field_spacing [0]) + ' μm)')
            ref_map = ax.imshow(O2_field, cmap=colour_map, origin='lower')
            fig.colorbar(ref_map, ax=ax, label='nM')
            plt.show()
        flat_field = middle_field['oxygen'].to_numpy().flatten()
#        np.save('/scratch/narain/Stochastic Pruning in Hexagonal Network with 100 Selections/TestHexagonalNetwork/' + solver_name + 'Haematocrit_Selection' + layout_selection + '_Sigma' + sigma + '_Beta' + beta + '_Distribution.npy', flat_field)
#    else:
#        flat_field = np.load('/scratch/narain/Stochastic Pruning in Hexagonal Network with 100 Selections/TestHexagonalNetwork/' + solver_name + 'Haematocrit_Selection' + layout_selection + '_Sigma' + sigma + '_Beta' + beta + '_Distribution.npy')

    # Get the basic stats 
#    middle_O2_stats = pd.DataFrame(flat_field).describe()
            
    # Get the number of total points
    number_of_points = 1
    for dim in np.shape(flat_field): number_of_points *= dim

    # Calculate the number of points below the given hypoxic thresholds
    hypoxic_fraction_list = []
    for hypoxic_threshold in hypoxic_threshold_list:
        hypoxic_points = (flat_field < hypoxic_threshold).sum()
        hypoxic_fraction = hypoxic_points/number_of_points  # calculate the hypoxic fraction
#        print(hypoxic_threshold, hypoxic_points)
        hypoxic_fraction_list.append(hypoxic_fraction)

    # Return the stats
    return hypoxic_fraction_list, np.mean(flat_field), np.amin(flat_field), np.percentile(flat_field, 50), np.amax(flat_field), np.std(flat_field), n_vessels, mean_diameter, mean_geometric_resistance, diameter_adjacency_matrix, length_adjacency_matrix
    '''
    
    return 0, 0, 0, 0, 0, 0, n_vessels, n_perfused_vessels, n_unperfused_vessels, mean_diameter, mean_geometric_resistance, diameter_binned, flow_rate_binned, n_connected_components, n_cycles, size_largest_connected_component, n_cycles_largest_connected_component


# Define a function to return statistics for all the heterogeneities in the data
def get_solver_stats(solver_name, alpha_value, mean_list, radius_threshold_list, hypoxic_threshold_list, plot, read):
#    table = np.array([])
    alpha_table = np.array([])
    for mean_value in mean_list:
        for radius_threshold in radius_threshold_list:    
            average_radius_threshold_data = get_distribution_stats(solver_name, alpha_value, mean_value, radius_threshold, hypoxic_threshold_list, diameter_threshold_list, flow_rate_threshold_list, plot, read)
            table_entry = np.hstack([float(mean_value), float(radius_threshold), average_radius_threshold_data])
            alpha_table = np.vstack([alpha_table, table_entry]) if alpha_table.size else table_entry
    return alpha_table

# =============================================================================
# DISTRIBUTION STATS & HYPOXIC FRACTIONS
# =============================================================================

# Enter details to allow looping over folders
n_gens = 6
d_inlet = 75    
pq_threshold = 1.e-12
which_mean = 0

# Set radius threshold
max_radius_threshold = 25
step_size = 0.02
#radius_threshold_list = [str(x) for x in range(max_radius_threshold + 1)]
radius_threshold_list = [str(format(x, '.2f')) for x in np.arange(0, max_radius_threshold+step_size, step_size)]

# Set folder path
main_folder_path = '/home/narain/Desktop/Results/Forking/Nonrandom/Radius Threshold Pruning in Nonrandom Forking Network with Varying Means (' + str(n_gens) + ' Gens ' + str(d_inlet) + ' Inlet)/TestDichotomousNetwork/'
#main_folder_path = '/tmp/narain/testoutput/TestDichotomousNetwork/'

# Set simulation specifics
solver_list = ['Constant', 'Pries', 'Memory', 'Fung']
alpha_list = ['1.00', '1.10', '1.20', '1.30', '1.40']
graph_alpha_list = ['1.0', '1.1', '1.2', '1.3', '1.4']  # what we want displayed on the graph
mean_list = ['0', '1', '2']
graph_mean_list = ['min', 'mean', 'max']  # what we want displayed on the graph
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
alpha_value  = alpha_list[1]

# Set file name for saving plots
description = '_nonrandom_forking_threshold_pruning_geometric_metrics_with_alpha_' + alpha_value

# Get the stats for all solvers (change to read=1 to extract from .vti files directly)
solver_stats = get_solver_stats(solver_name, alpha_value, mean_list, radius_threshold_list, hypoxic_threshold_list, plot=0, read=0)

# Save array
np.save(main_folder_path + solver_name + 'Haematocrit/python_solver_data.npy', solver_stats)
#solver_stats = np.load(main_folder_path + solver_name + 'Haematocrit/python_solver_data.npy')

# Filter by mean
mean_array, hypoxic_fraction_composite, mean_composite, min_composite, half_composite, max_composite, sd_composite, n_vessels_composite, n_perfused_vessels_composite, n_unperfused_vessels_composite, mean_diameter_composite, mean_geometric_resistance_composite, diameter_binned_composite, flow_rate_binned_composite, n_connected_components_composite, n_cycles_composite, size_largest_connected_component_composite, n_cycles_largest_connected_component_composite = filter_by_alpha(mean_list, solver_stats)

# =============================================================================
# PERFUSION QUOTIENTS
# =============================================================================

# Read PQ file
filename = main_folder_path + 'forking_nonrandom_threshold_pruning_perfusion_quotients.txt'
pq_df = pd.read_csv(filename, delim_whitespace=True, names=["network_name", "solver_name", "lambda", "alpha", "mean", "radius_threshold",  "PQ"])
#pq_df = pd.read_csv(filename, delim_whitespace=True, names=["alpha", "beta", "PQ"], skiprows=1)

# Filter PQ data for multiple solvers and alphas
solver_filter = solver_name + 'Haematocrit/'
pq_df = pq_df.loc[(pq_df["solver_name"] == solver_filter)]
pq_df = pq_df.loc[(pq_df["alpha"] == float(alpha_value))]

# Drop extra data
pq_df = pq_df.loc[(pq_df["radius_threshold"] <= max_radius_threshold)]

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
# PLOTS FOR GEOMETRIC METRICS
# =============================================================================
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
    axs[i].title.set_text('${µ}$ = ' + graph_mean_list[i])
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
    axs[i].set_xlabel('radius threshold (μm)')    
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
    axs[i].set_xlabel('radius threshold (μm)')    
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
    axs[i].set_xlabel('radius threshold (μm)')    
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
    axs[i].set_xlabel('radius threshold (μm)')    
    if i==len(mean_list*row_index):
        axs[i].set_ylabel('$β_0$') 
    axs[i].grid()
row_index+=1
    
# Plot the number of loops
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
    axs[i].plot(mean_array[:,1], n_cycles_composite[i-offset*row_index])
    axs[i].set_xlim(0)
#    axs[i].set_ylim([0,300])
    axs[i].set_xlabel('radius threshold (μm)')    
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
    axs[i].set_xlabel('radius threshold (μm)')    
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
    axs[i].set_xlabel('radius threshold (μm)')    
    if i==len(mean_list*row_index):
        axs[i].set_ylabel(r'$\overline{R}_β^{geom} (\times 10^{-2})$') 
    axs[i].grid()

'''
# =============================================================================
# COMBINED PLOTS
# =============================================================================
#'''
# Set the figure layout
fig, axs = plt.subplots(3, 1, figsize=(5, 15), tight_layout = {'pad': 2})
fig.subplots_adjust(hspace = .5, wspace=.25)
linestyles = ['solid', 'dashed', 'dotted', 'dashdot', (0,(5,10)), (0, (3, 1, 1, 1, 1, 1))]
linecolours = ['#1f77b4','r']
#plt.suptitle(solver_name + ' haematocrit solver in the heterogeneous hexagonal vessel network with stochastic pruning')

# Plot the distribution stats for a solver
axs = axs.ravel()
offset = len(mean_list)
row_index = 0

# Plot the PQ for a solver
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
#for i in range(len(mean_list),len(mean_list)*(row_index+1)):
#    axs[row_index].set_ylim([0,100])  # set PQ limits
    axs[row_index].plot(mean_array[:,1], (pq_composite[i-offset*row_index,:]-pq_composite[i-offset*row_index,1])*100/pq_composite[i-offset*row_index,1], ls = linestyles[i-offset*row_index], label = '${µ}$ = ' + graph_mean_list[i-offset*row_index])
#    axs[row_index].plot(mean_array[:,1], (pq_composite[int(mean_list[0]),:]-pq_composite[int(mean_list[0]),1])*100/pq_composite[int(mean_list[0]),1], label = '${µ}$ = ' + graph_mean_list[int(mean_list[0])], c=colors[0])

    axs[row_index].set_xlim(0)
#    axs[row_index].set_xlabel('radius threshold (μm)')    
    axs[row_index].set_ylabel('change in PQ (%)') 
    axs[row_index].legend()
    axs[row_index].grid()
    axs[row_index].title.set_text('${α}$ = ' + alpha_value)
row_index+=1
#'''
'''
# Plot the mean diameter for a solver
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
#    axs[row_index].set_ylim([0,1.1])  # set PQ limits
#    axs[row_index].plot(mean_array[:,1], mean_diameter_composite[i-offset*row_index], ls = linestyles[i-offset*row_index], label = '${µ}$ = ' + graph_mean_list[i-offset*row_index])
#    axs[row_index].set_xlim(0)
#    axs[row_index].set_ylim([0,100])
##    axs[row_index].set_xlabel('radius threshold (μm)')    
#    axs[row_index].set_ylabel('$\overline{D}$ (μm)') 
##    axs[row_index].legend()
#    axs[row_index].grid()
row_index+=1
'''
'''
# Plot the vessel composition for a solver
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
    axs[row_index].plot(mean_array[:,1], n_vessels_composite[i-offset*row_index], label='total')
    axs[row_index].plot(mean_array[:,1], n_perfused_vessels_composite[i-offset*row_index], ls='dashed', label='perfused')
    axs[row_index].plot(mean_array[:,1], n_unperfused_vessels_composite[i-offset*row_index], ls='dotted', label='unperfused')
row_index+=1
'''
#'''
# Plot the vessel composition for a solver
x=which_mean
axs[row_index].plot(mean_array[:,1], n_vessels_composite[x], label='total')
axs[row_index].plot(mean_array[:,1], n_perfused_vessels_composite[x], ls='dashed', label='perfused')
axs[row_index].plot(mean_array[:,1], n_unperfused_vessels_composite[x], ls='dotted', label='unperfused')
axs[row_index].legend()
axs[row_index].set_xlim(0)
axs[row_index].grid()
row_index+=1
#'''
'''
# Plot the mean geometric resistance for a solver
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
#    axs[row_index].set_ylim([0,1.1])  # set PQ limits
    axs[row_index].plot(mean_array[:,1], mean_geometric_resistance_composite[i-offset*row_index]*10000, ls = linestyles[i-offset*row_index], label = '${µ}$ = ' + graph_mean_list[i-offset*row_index])
    axs[row_index].set_xlim(0)
#    axs[row_index].set_ylim([0,1])
#    axs[row_index].set_xlabel('radius threshold (μm)')    
    axs[row_index].set_ylabel(r'$\overline{R}^{geom} (\times 10^{-4})$') 
#    axs[row_index].legend()
    axs[row_index].grid()
    axs[row_index].set_yscale('log')
row_index+=1
'''
'''
# Plot the loops per vessel
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
    loops_per_vessel = n_cycles_composite[i-offset*row_index]/n_vessels_composite[i-offset*row_index]
    axs[row_index].plot(mean_array[:,1], loops_per_vessel*100, ls = linestyles[i-offset*row_index], label = '${µ}$ = ' + graph_mean_list[i-offset*row_index])
    axs[row_index].set_xlim(0)
#    axs[row_index].set_ylim([0,.1])
#    axs[row_index].set_xlabel('radius threshold (μm)')    
    axs[row_index].set_ylabel(r'$\overline{β_1} (\times 10^{-2})$') 
    axs[row_index].grid()
row_index+=1
'''
# Plot the resistance per loop
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
    
    a = mean_geometric_resistance_composite[i-offset*row_index]
    b = n_cycles_composite[i-offset*row_index]/n_vessels_composite[i-offset*row_index]
    resistance_per_loop = np.divide(a, b, out=np.zeros_like(a), where=b!=0)
    axs[row_index].plot(mean_array[:,1], resistance_per_loop*100, ls = linestyles[i-offset*row_index], label = '${µ}$ = ' + graph_mean_list[i-offset*row_index])
    axs[row_index].set_xlim(0)
#    axs[row_index].set_ylim([0,.1])
    axs[row_index].set_xlabel('radius threshold (μm)')    
    axs[row_index].set_ylabel(r'$\overline{R}_β^{geom} (\times 10^{-2})$') 
    axs[row_index].grid()
    axs[row_index].set_yscale('log')

# Show plots
plt.show()

# Save image
file_path = Path('~/Desktop/Final Figures/' + solver_name + description + '.png').expanduser()
fig.savefig(file_path, dpi=250, bbox_inches = 'tight')
file_path = Path('~/Desktop/Final Figures/' + solver_name + description + '.png').expanduser()
fig.savefig(file_path, dpi=250, bbox_inches = 'tight')
#'''
# =============================================================================
#%% POSTER PLOTS
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

# Plot the PQ for a solver
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
#for i in range(len(mean_list),len(mean_list)*(row_index+1)):
#    axs[row_index].set_ylim([0,100])  # set PQ limits
    axs[row_index].plot(mean_array[:,1], (pq_composite[int(mean_list[0]),:]-pq_composite[int(mean_list[0]),1])*100/pq_composite[int(mean_list[0]),1], label = '${µ}$ = ' + graph_mean_list[int(mean_list[0])], c=colors[0])
    axs[row_index].set_xlim(0)
#    axs[row_index].set_xlabel('radius threshold (μm)')    
    axs[row_index].set_ylabel('change in PQ (%)') 
#    for index in range(0,len(line_changes)):
#        if line_changes[index]>0:
    axs[row_index].axvline(61, c='grey', alpha=0.5)  # highlight flow-rerouting
#    axs[row_index].text(60.1,0,'A',rotation=90)
    axs[row_index].axvspan(91, 129, facecolor='grey', alpha=0.5)  # highlight selective pruning
    axs[row_index].axvline(157, c='grey', alpha=0.5)  # highlight loops per size
#    axs[row_index].legend()
#    axs[row_index].grid()
    axs[row_index].spines['right'].set_visible(False)
    axs[row_index].spines['top'].set_visible(False)
#    axs[row_index].title.set_text('${α}$ = ' + alpha_value)
row_index+=1

# Plot the number of vessels
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
    axs[row_index].plot(mean_array[:,1], n_vessels_composite, ls='dashed', label='total', c=colors[0])
    axs[row_index].plot(mean_array[:,1], n_perfused_vessels_composite, label='perfused', c=colors[2])
    axs[row_index].plot(mean_array[:,1], n_unperfused_vessels_composite, label='unperfused', c=colors[1])
    axs[row_index].axvline(61, c='grey', alpha=0.5)  # highlight flow-rerouting
    axs[row_index].axvspan(91, 129, facecolor='grey', alpha=0.5)  # highlight selective pruning
    axs[row_index].set_xlim(0)
    axs[row_index].set_ylim(0,300)
#    axs[i].ticklabel_format(axis="y", scilimits=(0,0))
#    axs[i].set_xlabel('radius threshold (μm)')    
    axs[row_index].set_ylabel('network composition (vessels)') 
    axs[row_index].legend()
#    axs[row_index].grid()
    axs[row_index].spines['right'].set_visible(False)
    axs[row_index].spines['top'].set_visible(False)
#    axs[row_index].title.set_text('${µ}$ = ' + graph_mean_list[row_index])
row_index+=1

# Plot the loops per vessel
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
    loops_per_vessel = n_cycles_composite/n_vessels_composite
    axs[row_index].plot(mean_array[:,1], loops_per_vessel*100, label = '${µ}$ = ' + graph_mean_list[i-offset*row_index], c=colors[0])
    axs[row_index].set_xlim(0)
    axs[row_index].axvline(157, c='grey', alpha=0.5)  # highlight loops per size
#    axs[row_index].set_ylim([0,.1])
    axs[row_index].set_xlabel('radius threshold (μm)')    
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

# =============================================================================
# PAPER PLOTS
# =============================================================================

# Plot the mean diameter vs the resistance



# Save image
figure_title = 'summary_forking'
file_path = Path('~/Desktop/Final Figures/' + figure_title + '.svg').expanduser()
fig.savefig(file_path, dpi=500, bbox_inches = 'tight')
file_path = Path('~/Desktop/Final Figures/' + figure_title + '.png').expanduser()
fig.savefig(file_path, dpi=500, bbox_inches = 'tight')
'''

# Prints execution time
print("\n--- Execution Time: %s seconds ---" % (time.time() - start_time))
