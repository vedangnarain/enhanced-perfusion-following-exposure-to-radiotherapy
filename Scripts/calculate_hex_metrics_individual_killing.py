#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 16:58:31 2022

@author: narain
  
Tested in Python 3.7.4.

Calculate metrics for the hexagonal network with log-normally distributed radii and individual vessel pruning.
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
import scipy.io
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

# =============================================================================
# FUNCTIONS
# =============================================================================

# Define a function to read an O2 distribution .vti file and return the basic stats
def get_distribution_stats(solver_name, layout_selection, sigma, kills, hypoxic_threshold_list, diameter_threshold_list, flow_rate_threshold_list, plot=0, read=0):    
    
    # Set the file path
#    folder_path = '/home/narain/Desktop/Stochastic Pruning with 100 Trials/' + solver_name + 'Haematocrit/Lambda4/Alpha' + alpha_value + '/Beta' + beta + '/Trial' + trial
    folder_path = main_folder_path + solver_name + 'Haematocrit/Mu' + sigma + '/Selection' + layout_selection + '/Kills' + kills
    field_path =  folder_path + '/oxygen_solution_0.vti'    
    network_path = folder_path + '/FinalHaematocrit.vtp'

    # Print status update
    print(field_path)

    # Save file or read from existing file
    if read==0:    
        
        # Compute the architectural metrics
        n_vessels, n_perfused_vessels, n_unperfused_vessels, mean_diameter, mean_geometric_resistance, diameter_binned, flow_rate_binned, n_connected_components, n_cycles, size_largest_connected_component, n_cycles_largest_connected_component = get_hex_predictors(network_path, vessel_length_m, reference_node_coordinates, inlet_radius_m, pq_threshold, diameter_threshold_list, flow_rate_threshold_list)
    '''    
    # Write the adjacency matrices to .mat files
    matrix_name = 'mu' + sigma + '_selection' + layout_selection + '_kills' + kills  
    scipy.io.savemat(main_folder_path + 'DiameterAdjacencyMatrices/' + matrix_name + '_diameter_adjacency_matrix' + '.mat', {matrix_name+ '_diameter_adjacency_matrix': diameter_adjacency_matrix})
    scipy.io.savemat(main_folder_path + 'LengthAdjacencyMatrices/' + matrix_name + '_length_adjacency_matrix' + '.mat', {matrix_name+ '_length_adjacency_matrix': length_adjacency_matrix})
    '''
    return 0, 0, 0, 0, 0, 0, n_vessels, n_perfused_vessels, n_unperfused_vessels, mean_diameter, mean_geometric_resistance, diameter_binned, flow_rate_binned, n_connected_components, n_cycles, size_largest_connected_component, n_cycles_largest_connected_component
     
'''
        # Compute the O2 distribution
        field_data, field_spacing = get_vti_data(field_path)  # import the data for the network
        middle_field = get_hex_domain(field_data, field_spacing)  # get the distribution in the middle of the field (replace with designated file)
        
        # Plot the O2 distribution if desired
        if plot==1:  # plot the O2 distribution
            O2_field, _ = get_plottable_field(middle_field)  # get the O2 mesh for plots and stats
            fig = plt.figure()
            ax = plt.axes()
            colour_map = plt.cm.get_cmap('jet')
            plt.suptitle('O$_2$ distribution generated by the ' + solver_name + ' solver in the hexagonal vessel network with radii selection = ' + layout_selection + ' (σ = ' + sigma + ' and kills  = ' + kills + ' vessels (1 unit = ' + str(field_spacing [0]) + ' μm)')
            ref_map = ax.imshow(O2_field, cmap=colour_map, origin='lower')
            fig.colorbar(ref_map, ax=ax, label='nM')
            plt.show()
        flat_field = middle_field['oxygen'].to_numpy().flatten()
#        np.save('/home/narain/Temporary Python Files/Radius Threshold Pruning/' + solver_name + 'Haematocrit_Selection' + layout_selection + '_Sigma' + sigma + '_RadiusThreshold' + threshold + '_Distribution.npy', flat_field)
#    else:
#        flat_field = np.load('/home/narain/Temporary Python Files/Radius Threshold Pruning/' + solver_name + 'Haematocrit_Selection' + layout_selection + '_Sigma' + sigma + '_RadiusThreshold' + threshold + '_Distribution.npy')

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
    return hypoxic_fraction_list, np.mean(flat_field), np.amin(flat_field), np.percentile(flat_field, 50), np.amax(flat_field), np.std(flat_field)
'''

# Define function to compute the average of all layouts in a kills selection
def compute_average_kills(solver_name, sigma, kills, max_layouts, hypoxic_threshold_list, diameter_threshold_list, flow_rate_threshold_list, plot, read):
    
    # Create table to store all the kills data trials in an mean group
    kills_table = np.array([])

    # Extract metrics from all trials and store in table
    for layout_selection in range(1, max_layouts+1):    
#        print(layout_selection)
        hypoxic_fraction_list, mean_value, min_value, half_value, max_value, std_value, n_vessels_value, n_perfused_vessels_value, n_unperfused_vessels_value, mean_diameter_value, mean_geometric_resistance_value, diameter_binned_value, flow_rate_binned_value, n_connected_components_value, n_cycles_value, size_largest_connected_component_value, n_cycles_largest_connected_component_value = get_distribution_stats(solver_name, str(layout_selection), sigma, str(kills), hypoxic_threshold_list, diameter_threshold_list, flow_rate_threshold_list, plot=0, read=0)
        table_entry = np.hstack([hypoxic_fraction_list, mean_value, min_value, half_value, max_value, std_value, n_vessels_value, n_perfused_vessels_value, n_unperfused_vessels_value, mean_diameter_value, mean_geometric_resistance_value, diameter_binned_value, flow_rate_binned_value, n_connected_components_value, n_cycles_value, size_largest_connected_component_value, n_cycles_largest_connected_component_value])
        kills_table = np.vstack([kills_table, table_entry]) if kills_table.size else table_entry

    # Return the hypoxic fractions, mean, min, 50%, max, and std for a kills averaged across all layouts
    return np.average(kills_table, axis=0)
#    return np.average(threshold_table[0]), np.average(threshold_table[1]), np.average(threshold_table[2]), np.average(threshold_table[3]), np.average(threshold_table[4])

# Define a function to return statistics for all the heterogeneities in the data
def get_solver_stats(solver_name, mean_list, kills_list, max_layouts, hypoxic_threshold_list, plot, read):
#    table = np.array([])
    mean_table = np.array([])
    for sigma in mean_list:
        for kills in kills_list:    
            average_kills_data = compute_average_kills(solver_name, sigma, kills, max_layouts, hypoxic_threshold_list, diameter_threshold_list, flow_rate_threshold_list, plot, read)
            table_entry = np.hstack([float(sigma), float(kills), average_kills_data])
            mean_table = np.vstack([mean_table, table_entry]) if mean_table.size else table_entry
    return mean_table

# =============================================================================
# DISTRIBUTION STATS & HYPOXIC FRACTIONS
# =============================================================================

# Enter details to allow looping over folders
sd_list = ['8.68', '13.23', '17.49']
solver_list = ['Constant', 'Pries', 'Memory', 'Fung']
mean_list = ['22.76', '28.5', '33.64']
max_kills = 200
max_layouts = 100
kills_list = [str(x) for x in range(0, max_kills + 1)]
#hypoxic_threshold_list = [2195, 10000, 15000, 20000, 25000, 27441] 
#hypoxic_threshold_list = [2195, 5488, 10976, 16465, 21953, 27441] 
hypoxic_threshold_list = [2195, 27441] 
#hypoxic_threshold_list_pp = [0.8, 2, 4, 6, 8, 10] 
graph_mean_list = ['22.76 μm', '28.50 μm', '33.65 μm']  # what we want displayed on the graph

# Set network details (in metres)
vessel_length_m = 100*(10**-6)
inlet_radius_m = 7.500000e-05
#pq_threshold = 1.e-13
pq_threshold = 3.e-13
diameter_threshold_list = [22, 35, 50]  # in um
flow_rate_threshold_list = [1e-13, 2e-12, 5e-12]

# Set the folder path
solver_name = solver_list[0]
which_sd = 0
sd_folder = sd_list[which_sd]
which_mean = 2
# main_folder_path = '/scratch/narain/Hexagonal/Log Normal Distribution/Individual Pruning in Hexagonal Network with 100 Selections and Varying Means with SD of ' + sd_folder + '/TestHexagonalNetwork/'
main_folder_path = '/home/narain/Desktop/Results/Hexagonal/Without Oxygen/Individual Pruning in Hexagonal Network with 100 Selections and Varying Means with SD of ' + sd_folder + '/TestHexagonalNetwork/'
# main_folder_path = '/tmp/narain/testoutput/TestHexagonalNetwork/'

# Get reference hexagonal network (all node coordinates are the same, regardless of heterogeneity)
reference_node_coordinates = get_reference_hexagonal_network('/home/narain/Desktop/Scripts/reference_networks/square_hexagonal_network_100_um/Selection1/RadiusThreshold0/FinalHaematocrit.vtp')
#reference_node_coordinates = get_reference_hexagonal_network('/home/narain/Desktop/Scripts/reference_networks/single_feed_hexagonal_network_100_um/FinalHaematocrit.vtp')

# Set file name for saving plots
description = '_lognormal_hexagonal_individual_pruning_geometric_metrics_sigma_' + sd_folder

# Get the stats for all solvers (change to read=1 to extract from .vti files directly)
# solver_stats = get_solver_stats(solver_name, mean_list, kills_list, max_layouts, hypoxic_threshold_list, plot=0, read=0)

# Save array
# np.save(main_folder_path + solver_name + 'Haematocrit/python_solver_data.npy', solver_stats)
solver_stats = np.load(main_folder_path + solver_name + 'Haematocrit/python_solver_data.npy')

# Filter by mean
mean_array, hypoxic_fraction_composite, mean_composite, min_composite, half_composite, max_composite, sd_composite, n_vessels_composite, n_perfused_vessels_composite, n_unperfused_vessels_composite, mean_diameter_composite, mean_geometric_resistance_composite, diameter_binned_composite, flow_rate_binned_composite, n_connected_components_composite, n_cycles_composite, size_largest_connected_component_composite, n_cycles_largest_connected_component_composite = filter_by_mean_hex(mean_list, solver_stats)

# =============================================================================
# PERFUSION QUOTIENTS
# =============================================================================

# Define a function to generate the data for a single mean value
def get_mean_line(mean_group, max_kills):
    
    # Compute the averages for the kills in the mean group
    mean_pq_table = np.array([])
    sd_pq_table = np.array([])

    mean_kills_grouped = mean_group.groupby(mean_group.kills)
    for kills in range(0, max_kills+1):
        mean_pq_table_entry = np.array([])
        sd_pq_table_entry = np.array([])
        mean_kills_group = mean_kills_grouped.get_group(kills)
        mean_pq_table_entry = np.array([mean_kills_group["mean"].mean(), mean_kills_group["kills"].mean(), mean_kills_group["PQ"].mean()])
        sd_pq_table_entry = np.array([mean_kills_group["mean"].mean(), mean_kills_group["kills"].mean(), mean_kills_group["PQ"].std()])
        mean_pq_table = np.vstack([mean_pq_table, mean_pq_table_entry]) if mean_pq_table.size else mean_pq_table_entry
        sd_pq_table = np.vstack([sd_pq_table, sd_pq_table_entry]) if sd_pq_table.size else sd_pq_table_entry
    
    # Return the table for reference
    return mean_pq_table[:,2], sd_pq_table[:,2]

# Read PQ file
filename = main_folder_path + 'hex_lognormal_individual_perfusion_quotients.txt'
pq_df = pd.read_csv(filename, delim_whitespace=True, names=["network_name", "solver_name", "mean", "selection", "kills", "PQ"])
#pq_df = pd.read_csv(filename, delim_whitespace=True, names=["mean", "threshold", "PQ"], skiprows=1)

# Filter PQ data for multiple solvers
solver_filter = solver_name + 'Haematocrit'
pq_df = pq_df.loc[(pq_df["solver_name"] == solver_filter)]
#pq_df = pq_df.loc[(pq_df["selection"] == 1)]

# Drop extra data
#max_beta = 35
#pq_df = pq_df.loc[(pq_df["beta"] <= max_beta)]

'''
# Separate by mean 
mean_grouped = pq_df.groupby(pq_df.mean)
#mean_0 = mean_grouped.get_group(0)
mean_1 = mean_grouped.get_group(1)
mean_2 = mean_grouped.get_group(2)
mean_3 = mean_grouped.get_group(3)
mean_4 = mean_grouped.get_group(4)
'''

# Separate by mean 
mean_grouped = pq_df.groupby(['mean'])
#mean_0 = mean_grouped.get_group(0)
mean_1 = mean_grouped.get_group(float(mean_list[0]))
mean_2 = mean_grouped.get_group(float(mean_list[1]))
mean_3 = mean_grouped.get_group(float(mean_list[2]))
#mean_4 = mean_grouped.get_group(float(mean_list[3]))

# Compute average of all selections for PQ
line_1, line_1_sd = get_mean_line(mean_1, max_kills)
line_2, line_2_sd = get_mean_line(mean_2, max_kills)
line_3, line_3_sd = get_mean_line(mean_3, max_kills)
#line_4, line_4_sd = get_mean_line(mean_4, max_kills)

# Combine the PQs
#pq_composite = np.vstack([line_1, line_2, line_3 ,line_4])
#sd_pq_composite = np.vstack([line_1_sd, line_2_sd, line_3_sd, line_4_sd])
pq_composite = np.vstack([line_1, line_2, line_3])
sd_pq_composite = np.vstack([line_1_sd, line_2_sd, line_3_sd])

# =============================================================================
# COMBINED PLOTS
# =============================================================================

linestyles = ['solid', 'dashed', 'dotted', 'dashdot', (0,(5,10)), (0, (3, 1, 1, 1, 1, 1))]
linecolours = ['#0072B2', '#009E73','#D55E00']
width = 5

'''
# Set the figure layout
fig, axs = plt.subplots(5, 1, figsize=(7, 20), tight_layout = {'pad': 2})
fig.subplots_adjust(hspace = .5, wspace=.25)
linestyles = ['solid', 'dashed', 'dotted', 'dashdot', (0,(5,10)), (0, (3, 1, 1, 1, 1, 1))]
linecolours = ['#1f77b4','r']
#plt.suptitle(solver_name + ' haematocrit solver in the heterogeneous hexagonal vessel network with stochastic pruning')

# Plot the distribution stats for a solver
axs = axs.ravel()
offset = len(mean_list)
row_index = 0

# Plot the PQ for a solver
for i in range(offset*row_index,offset*(row_index+1)):
#    axs[row_index].set_ylim([0,25])  # set PQ limits
#    axs[row_index].set_ylim([0,1.1])  # set PQ limits
    axs[row_index].plot(mean_array[:,1], (pq_composite[i-offset*row_index,:]-pq_composite[i-offset*row_index,2])*100/pq_composite[i-offset*row_index,2], ls = linestyles[i-offset*row_index], label = r'$\overline{d}$ = ' + mean_list[i-offset*row_index])
#    axs[row_index].plot(mean_array[:,1], pq_composite[i-offset*row_index,:], ls = linestyles[i-offset*row_index], label = r'$\overline{d}$ = ' + mean_list[i-offset*row_index])
#    axs[row_index].fill_between(mean_array[:,1], pq_composite[:,(3*i)-10]+sd_pq_composite[:, (3*i)-10], pq_composite[:,(3*i)-10]-sd_pq_composite[:, (3*i)-10], color='grey', alpha=0.5, label='mean ± SD')
    axs[row_index].set_xlim(0)
#    axs[row_index].set_xlabel('radius threshold (μm)')    
#    axs[row_index].set_ylabel('PQ') 
    axs[row_index].set_ylabel('Relative % PQ') 
    axs[row_index].grid(True)
    axs[row_index].title.set_text('SD = ' + sd_folder)
#    axs[row_index].legend()
row_index+=1
    
# Plot the mean diameter for a solver
for i in range(offset*row_index,offset*(row_index+1)):
#    axs[row_index].set_ylim([0,1.1])  # set PQ limits
    axs[row_index].plot(mean_array[:,1], mean_diameter_composite[i-offset*row_index], ls = linestyles[i-offset*row_index], label = r'$\overline{d}$ = ' + mean_list[i-offset*row_index])
    axs[row_index].set_xlim(0)
#    axs[row_index].set_ylim([0,100])
#    axs[row_index].set_xlabel('radius threshold (μm)')    
    axs[row_index].set_ylabel('$\overline{D}$ (μm)') 
#    axs[row_index].legend()
    axs[row_index].grid(True)
row_index+=1
    
# Plot the mean geometric resistance for a solver
for i in range(offset*row_index,offset*(row_index+1)):
#    axs[row_index].set_ylim([0,1.1])  # set PQ limits
#    print(i,offset,i-offset)
    axs[row_index].plot(mean_array[:,1], mean_geometric_resistance_composite[i-offset*row_index]*10000, ls = linestyles[i-offset*row_index], label = r'$\overline{d}$ = ' + mean_list[i-offset*row_index])
    axs[row_index].set_xlim(0)
#    axs[row_index].set_ylim([0,4])
#    axs[row_index].set_xlabel('kills (vessels)')    
    axs[row_index].set_ylabel(r'$\overline{R}^{geom} (\times 10^{-4})$') 
#    axs[row_index].legend()
    axs[row_index].grid(True)
row_index+=1

# Plot the loops per vessel
for i in range(offset*row_index,offset*(row_index+1)):
    loops_per_vessel = n_cycles_composite[i-offset*row_index]/n_vessels_composite[i-offset*row_index]
    axs[row_index].plot(mean_array[:,1], loops_per_vessel*100, ls = linestyles[i-offset*row_index], label = r'$\overline{d}$ = ' + mean_list[i-offset*row_index])
    axs[row_index].set_xlim(0)
#    axs[row_index].set_ylim([0,.1])
#    axs[row_index].set_xlabel('kills (vessels)')    
    axs[row_index].set_ylabel(r'$\overline{β_1} (\times 10^{-2})$') 
    axs[row_index].grid(True)
row_index+=1

# Plot the resistance per loop
for i in range(offset*row_index,offset*(row_index+1)):
    
    a = mean_geometric_resistance_composite[i-offset*row_index]
    b = n_cycles_composite[i-offset*row_index]/n_vessels_composite[i-offset*row_index]
    resistance_per_loop = np.divide(a, b, out=np.zeros_like(a), where=b!=0)
    axs[row_index].plot(mean_array[:,1], resistance_per_loop*100, ls = linestyles[i-offset*row_index], label = r'$\overline{d}$ = ' + mean_list[i-offset*row_index])
    axs[row_index].set_xlim(0)
#    axs[row_index].set_ylim([0,.1])
    axs[row_index].set_xlabel('kills (vessels)')    
    axs[row_index].set_ylabel(r'$\overline{R}_β^{geom} (\times 10^{-2})$') 
    axs[row_index].grid(True)
    
# Show plots
plt.show()
'''

# =============================================================================
# RAW + % CHANGE PQ PLOTS
# =============================================================================

mean_list = ['0', '1', '2']

'''
# Set the figure layout
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.5), tight_layout = {'pad': 2})
fig.subplots_adjust(hspace = .5, wspace=.25)
linestyles = ['dotted', 'dashdot', 'solid', 'dashed', (0,(5,10)), (0, (3, 1, 1, 1, 1, 1))]
linecolours = ['C0', 'g','r']
#plt.suptitle(solver_name + ' haematocrit solver in the heterogeneous hexagonal vessel network with stochastic pruning')
line_changes = np.diff(n_perfused_vessels_composite)

# Plot the distribution stats for a solver
axs = axs.ravel()
offset = len(mean_list)
row_index = 0


# Plot the raw PQ for a solver
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
    axs[row_index].set_ylim([0,1])  # set PQ limits
#    axs[row_index].plot(mean_array[:,1], (pq_composite[i-offset*row_index,:]-pq_composite[i-offset*row_index,1])*100/pq_composite[i-offset*row_index,1], ls = linestyles[i-offset*row_index], label = r'$\overline{d}$ = ' + graph_mean_list[i-offset*row_index])
    pq_composite[i-offset*row_index,0:2]=pq_composite[i-offset*row_index,2]  # fix flow rate error
    axs[row_index].plot(mean_array[:,1], (pq_composite[i-offset*row_index,:]), ls = linestyles[i-offset*row_index], c = linecolours[i-offset*row_index], label = graph_mean_list[i-offset*row_index])
    axs[row_index].set_xlim(0)
#    axs[row_index].set_ylabel('change in PQ (%)') 
    axs[row_index].set_ylabel('PQ') 
    axs[row_index].set_xlabel('dosage (vessels pruned)')    
#    for index in range(0,len(line_changes)):
#        if line_changes[index]>0:
#    axs[row_index].axvline(61, c='grey', alpha=0.5)  # highlight flow-rerouting
#    axs[row_index].text(60.1,0,'A',rotation=90)
#    axs[row_index].axvspan(91, 129, facecolor='grey', alpha=0.5)  # highlight selective pruning
#    axs[row_index].axvline(157, c='grey', alpha=0.5)  # highlight loops per size
    axs[row_index].legend()
    axs[row_index].grid()
#    axs[row_index].spines['right'].set_visible(False)
#    axs[row_index].spines['top'].set_visible(False)
    axs[row_index].title.set_text('${σ}$ = ' + sd_folder + ' µm')
row_index+=1

# '''
'''
# Plot the % PQ change for a solver
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
    axs[row_index].set_ylim([-100,10])  # set PQ limits
    pq_composite[i-offset*row_index,0:2]=pq_composite[i-offset*row_index,2]  # fix flow rate error
    axs[row_index].plot(mean_array[:,1], (pq_composite[i-offset*row_index,:]-pq_composite[i-offset*row_index,1])*100/pq_composite[i-offset*row_index,1], ls = linestyles[i-offset*row_index], c = linecolours[i-offset*row_index], label = graph_mean_list[i-offset*row_index])
#    axs[row_index].plot(mean_array[:,1], (pq_composite[i-offset*row_index,:]), ls = linestyles[i-offset*row_index], c = linecolours[i-offset*row_index], label = graph_mean_list[i-offset*row_index])
    axs[row_index].set_xlim(0)
#    axs[row_index].set_ylabel('change in PQ (%)') 
    axs[row_index].set_ylabel(r'$\Delta_{\%} \mathcal{P}$') 
    axs[row_index].set_xlabel('dosage (vessels pruned)')    
#    for index in range(0,len(line_changes)):
#        if line_changes[index]>0:
#    axs[row_index].axvline(61, c='grey', alpha=0.5)  # highlight flow-rerouting
#    axs[row_index].text(60.1,0,'A',rotation=90)
#    axs[row_index].axvspan(91, 129, facecolor='grey', alpha=0.5)  # highlight selective pruning
#    axs[row_index].axvline(157, c='grey', alpha=0.5)  # highlight loops per size
    axs[row_index].legend()
    axs[row_index].grid()
#    axs[row_index].spines['right'].set_visible(False)
#    axs[row_index].spines['top'].set_visible(False)
    axs[row_index].title.set_text('${σ}$ = ' + sd_folder + ' µm')
row_index+=1

# '''

'''
# Plot the loops per vessel
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
#    loops_per_vessel = n_cycles_composite/n_vessels_composite
    loops_per_vessel = n_cycles_composite[i-offset*row_index]/n_vessels_composite[i-offset*row_index]
    axs[row_index].plot(mean_array[:,1], loops_per_vessel*100, label = r'$\overline{d}$ = ' + graph_mean_list[i-offset*row_index], ls = linestyles[i-offset*row_index])
    axs[row_index].set_xlim(0)
#    axs[row_index].axvline(157, c='grey', alpha=0.5)  # highlight loops per size
#    axs[row_index].set_ylim([0,.1])
#    axs[row_index].set_xlabel('vessels killed')    
    axs[row_index].set_ylabel(r'$\overline{β_1} (\times 10^{-2})$') 
    axs[row_index].grid()    
#    axs[row_index].spines['right'].set_visible(False)
#    axs[row_index].spines['top'].set_visible(False)
row_index+=1
'''

'''
# Plot the resistance per loop
for i in range(len(mean_list)*row_index,len(mean_list)*(row_index+1)):
    
    a = mean_geometric_resistance_composite[i-offset*row_index]
    b = n_cycles_composite[i-offset*row_index]/n_vessels_composite[i-offset*row_index]
    resistance_per_loop = np.divide(a, b, out=np.zeros_like(a), where=b!=0)
    axs[row_index].plot(mean_array[:,1], resistance_per_loop*100, ls = linestyles[i-offset*row_index], label = r'$\overline{d}$ = ' + graph_mean_list[i-offset*row_index])
    axs[row_index].set_xlim(0)
#    axs[row_index].set_ylim([0,.1])
    axs[row_index].set_xlabel('vessels pruned')    
    axs[row_index].set_ylabel(r'$\overline{R}_β^{geom} (\times 10^{-2})$') 
    axs[row_index].legend()
    axs[row_index].grid(True)
   
# Show plots
plt.show()
'''

# =============================================================================
# NETWORK COMPOSITION
# =============================================================================
# '''
# Set the figure layout
fig, axs = plt.subplots(1, 1, figsize=(6.25, 5.5))
# fig.subplots_adjust(hspace = .5, wspace=.25)
linestyles = ['dotted', 'dashdot', 'solid', 'dashed', (0,(5,10)), (0, (3, 1, 1, 1, 1, 1))]
#plt.suptitle(solver_name + ' haematocrit solver in the heterogeneous hexagonal vessel network with stochastic pruning')

# Plot the distribution stats for a solver
# axs = axs.ravel()
offset = len(mean_list)
row_index = 0
#x = which_mean
#i = x

# Rectify flow rate error by using second value of array as the initial value
for mean_i in range(0,3):
    n_perfused_vessels_composite[mean_i,0:2]=n_perfused_vessels_composite[which_mean,2]
    n_unperfused_vessels_composite[mean_i,0:2]=n_unperfused_vessels_composite[which_mean,2]
    pq_composite[mean_i,0:2]=pq_composite[which_mean,2]  # fix flow rate error

ax2 = axs.twinx()
#axs.plot(mean_array[:,1], n_vessels_composite[x], ls='dotted', label='total')
axs.plot(mean_array[:,1], n_perfused_vessels_composite[which_mean], ls='dashed', label='perfused', lw=width, c=linecolours[1])
axs.plot(mean_array[:,1], n_unperfused_vessels_composite[which_mean], ls='dotted', label='hypoperfused', lw=width, c=linecolours[2])
axs.set_xlim(0)
axs.set_ylim(0,400)
axs.set_xlabel('dosage (vessels pruned)')    
axs.grid()
# pq_composite[:,0]=pq_composite[:,1]  # fix flow rate error
ax2.plot(mean_array[:,1], (pq_composite[which_mean,:]-pq_composite[which_mean,1])*100/pq_composite[which_mean,1], lw=width, c=linecolours[0], label='$\Delta_{\%} \mathcal{P}$')
ax2.set_ylim(-100,20)  
axs.title.set_text(r'$\overline{d}$ = ' + graph_mean_list[which_mean])  
if which_mean==0:
    axs.set_ylabel('composition (number of vessels)') 
    # axs.text(-0.5, 0.5, '${σ}$ = ' + sd_folder + ' µm', va='center', ha='right', transform=axs[which_sd, which_mean].transAxes)
    if which_sd==0:
        fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=axs.transAxes)
if which_mean==2:
    ax2.set_ylabel(r'$\Delta_{\%} \mathcal{P}$')

# Show plots
plt.show()

# '''
# Save image
#file_path = Path('~/Desktop/Final Figures/' + solver_name + description + '.svg').expanduser()
#fig.savefig(file_path, dpi=500, bbox_inches = 'tight')
#file_path = Path('~/Desktop/Final Figures/' + solver_name + description + '.png').expanduser()
file_path = Path('~/Desktop/Final Figures/' + 'sigma_' + sd_folder + '_mean_' + graph_mean_list[which_mean] + '.tif').expanduser()
fig.savefig(file_path, dpi=300, bbox_inches = 'tight')
#'''
# Prints execution time
print("\n--- Execution Time: %s seconds ---" % (time.time() - start_time))

# =============================================================================
# UNUSED FUNCTIONALITY
# =============================================================================

# =============================================================================
# PLOTS FOR O2 STATS, PQ, AND HF
# =============================================================================
'''
# Set the figure layout
fig, axs = plt.subplots(3, len(alpha_list), figsize=(20, 12), tight_layout = {'pad': 2})
fig.subplots_adjust(hspace = .5, wspace=.25)
linestyles = ['solid', 'dashed', 'dotted', 'dashdot', (0,(5,10)), (0, (3, 1, 1, 1, 1, 1))]
linecolours = ['#1f77b4','b','b','b','b', 'r']
linelegends = ['anoxia', 2, 4, 6, 8, 'hypoxia']
#plt.suptitle(solver_name + ' haematocrit solver in the heterogeneous hexagonal vessel network with radius threshold pruning')

# Plot the distribution stats for a solver
axs = axs.ravel()
for i in range(len(alpha_list)):
    axs[i].set_ylim([0,30000])
#    axs[i].plot(alpha_array[:,1], mean_composite[i], ls='dashed', label='mean')
#    axs[i].plot(alpha_array[:,1], min_composite[i], ls='dotted', label='min')
##    axs[i].plot(alpha_array[:,1], half_composite[i], ls=(0, (3, 5, 1, 5)), label='50%')
#    axs[i].plot(alpha_array[:,1], max_composite[i], ls='dashdot', label='max')
#    axs[i].plot(alpha_array[:,1], sd_composite[i], ls='solid', label='SD')
#    axs[i].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
#    axs[i].set_xlabel('radius threshold (μm)')    
    if i==0:
        axs[i].set_ylabel('oxygen (nM)') 
        axs[i].legend(loc="best", prop={'size': 15})
    axs[i].set_xlim(0)
#    axs[i].set_ylim(0)
    axs[i].grid()
    axs[i].title.set_text('${σ}$ = ' + alpha_list[i])

# Plot the PQ for a solver
for i in range(len(alpha_list),len(alpha_list)*2):
    axs[i].set_ylim([0,1.1])  # set PQ limits
    axs[i].plot(line_1[:,1], pq_composite[:, (3*i)-10], label='PQ')
#    axs[i].set_xlabel('radius threshold (μm)')    
    if i==len(alpha_list):
        axs[i].set_ylabel('PQ') 
    axs[i].set_xlim(0)
    axs[i].set_ylim(0)
#    axs[i].legend()
    axs[i].grid()

# Plot the HF for a solver
for i in range(len(alpha_list*2),len(alpha_list)*3):
    offset = i-len(alpha_list*2)
    for threshold_index in range(len(hypoxic_threshold_list)):  # plot a line for each threshold
        if threshold_index==0 or threshold_index==5:
            axs[i].plot(alpha_array[:,1], hypoxic_fraction_composite[offset*len(kills_list):(offset+1)*len(kills_list), threshold_index], label=linelegends[threshold_index], ls=linestyles[threshold_index], c=linecolours[threshold_index])
        if i==len(alpha_list*2):
            axs[i].set_ylabel('HF/AF') 
            axs[i].legend(loc="best", prop={'size': 15})
        axs[i].set_xlabel('kills (vessels)')    
        axs[i].set_xlim(0)
        axs[i].set_ylim([0,1.1])  # set HF limits
        axs[i].grid(b=True)

# Show plots
plt.show()

# Save image
file_path = Path('~/Desktop/Final Figures/' + solver_name + '_lognormal_hexagonal_individual_pruning.svg').expanduser()
fig.savefig(file_path, dpi=500, bbox_inches = 'tight')
file_path = Path('~/Desktop/Final Figures/' + solver_name + '_lognormal_hexagonal_individual_pruning.png').expanduser()
fig.savefig(file_path, dpi=500, bbox_inches = 'tight')
'''

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
for i in range(offset):
    axs[i].plot(mean_array[:,1], n_vessels_composite[i], label='total')
    axs[i].plot(mean_array[:,1], n_perfused_vessels_composite[i], ls='dashed', label='perfused')
    axs[i].plot(mean_array[:,1], n_unperfused_vessels_composite[i], ls='dotted', label='unperfused')
    axs[i].set_xlim(0)
    axs[i].set_ylim(0,500)
#    axs[i].ticklabel_format(axis="y", scilimits=(0,0))
#    axs[i].set_xlabel('radius threshold (μm)')    
    if i==0:
        axs[i].set_ylabel('${N}_V$') 
        axs[i].legend()
    axs[i].grid()
    axs[i].title.set_text(r'$\overline{d}$ = ' + mean_list[i])
row_index+=1

# Plot the PQ for a solver
for i in range(offset,offset*(row_index+1)):
    axs[i].set_ylim([0,1.1])  # set PQ limits
    axs[i].plot(mean_array[:,1], pq_composite[:,(3*i)-10])
    axs[i].fill_between(mean_array[:,1], pq_composite[:,(3*i)-10]+sd_pq_composite[:, (3*i)-10], pq_composite[:,(3*i)-10]-sd_pq_composite[:, (3*i)-10], color='grey', alpha=0.5, label='mean ± SD')
    axs[i].set_xlim(0)
#    axs[i].set_xlabel('radius threshold (μm)')    
    if i==offset:
        axs[i].set_ylabel('PQ') 
        axs[i].legend()
    axs[i].grid()
row_index+=1
    
# Plot the mean diameter for a solver
for i in range(offset*row_index,offset*(row_index+1)):
#    axs[i].set_ylim([0,1.1])  # set PQ limits
    print(i,offset,i-offset)
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
for i in range(offset*row_index,offset*(row_index+1)):
#    axs[i].set_ylim([0,1.1])  # set PQ limits
#    print(i,offset,i-offset)
    axs[i].plot(mean_array[:,1], mean_geometric_resistance_composite[i-offset*row_index]*10000)
    axs[i].set_xlim(0)
#    axs[i].set_ylim([0,4])
    axs[i].set_xlabel('kills (vessels)')    
    if i==len(mean_list*row_index):
        axs[i].set_ylabel(r'$\overline{\mathcal{R}}^{geom}$ (\times 10^{-4})$') 
#    axs[i].legend()
    axs[i].grid()
row_index+=1
'''
'''
# Plot the diameter composition
for i in range(offset*row_index,offset*(row_index+1)):
#    offset_2 = i-len(mean_list*row_index)
    diameter_bin_labels = ['< '+str(diameter_threshold_list[0])+' μm', str(diameter_threshold_list[0])+' μm – '+str(diameter_threshold_list[1])+' μm', str(diameter_threshold_list[1])+' μm – '+str(diameter_threshold_list[2])+' μm', '> '+str(diameter_threshold_list[2])+' μm']
#    for mean_set_index in range(offset):  # plot a line for each threshold
    offset_2 = (i-len(mean_list*row_index))
    mean_set_1 = [diameter_binned_composite[offset_2*len(kills_list):(offset_2+1)*len(kills_list),0]]
    mean_set_2 = [diameter_binned_composite[offset_2*len(kills_list):(offset_2+1)*len(kills_list),1]]
    mean_set_3 = [diameter_binned_composite[offset_2*len(kills_list):(offset_2+1)*len(kills_list),2]]
#        mean_set_4 = [diameter_binned_composite[mean_set_index*len(kills_list):(mean_set_index+1)*len(kills_list),3]]
    highest_bracket = np.array(n_vessels_composite[i-offset*row_index] - mean_set_1 - mean_set_2 - mean_set_3)
    axs[i].stackplot(np.array(mean_array[:,1], dtype=int), mean_set_1, mean_set_2, mean_set_3, highest_bracket, labels=[diameter_bin_labels[0],diameter_bin_labels[1],diameter_bin_labels[2],diameter_bin_labels[3]])
    axs[i].set_xlim(0)
    axs[i].set_ylim([0,500])
    axs[i].set_xlabel('kills (vessels)')    
    if i==len(mean_list*row_index):
        axs[i].set_ylabel('${N}_V$') 
        handles, labels = axs[i].get_legend_handles_labels()
        axs[i].legend(reversed(handles), reversed(labels), title='Diameter')    
    axs[i].grid()
row_index+=1
    
# Plot the flow rate composition
for i in range(offset*row_index,offset*(row_index+1)):
#    offset_2 = i-len(mean_list*4)
    flow_rate_bin_labels = ['< '+str(flow_rate_threshold_list[0])+' $m^{3}$/s', str(flow_rate_threshold_list[0])+' $m^{3}$/s – '+str(flow_rate_threshold_list[1])+' $m^{3}$/s', str(flow_rate_threshold_list[1])+' $m^{3}$/s – '+str(flow_rate_threshold_list[2])+' $m^{3}$/s', '> '+str(flow_rate_threshold_list[2])+' $m^{3}$/s']
#    for mean_set_index in range(offset):  # plot a line for each threshold
    offset_2 = (i-len(mean_list*row_index))
    mean_set_1 = [flow_rate_binned_composite[offset_2*len(kills_list):(offset_2+1)*len(kills_list),0]]
    mean_set_2 = [flow_rate_binned_composite[offset_2*len(kills_list):(offset_2+1)*len(kills_list),1]]
    mean_set_3 = [flow_rate_binned_composite[offset_2*len(kills_list):(offset_2+1)*len(kills_list),2]]
#        mean_set_4 = [diameter_binned_composite[mean_set_index*len(kills_list):(mean_set_index+1)*len(kills_list),3]]
    highest_bracket = np.array(n_vessels_composite[i-offset*row_index] - mean_set_1 - mean_set_2 - mean_set_3)
    axs[i].stackplot(np.array(mean_array[:,1], dtype=int), mean_set_1, mean_set_2, mean_set_3, highest_bracket, labels=[flow_rate_bin_labels[0],flow_rate_bin_labels[1],flow_rate_bin_labels[2],flow_rate_bin_labels[3]])
    axs[i].set_xlim(0)
    axs[i].set_ylim([0,500])
    axs[i].set_xlabel('kills (vessels)')    
    if i==len(mean_list*row_index):
        axs[i].set_ylabel('${N}_V$') 
        handles, labels = axs[i].get_legend_handles_labels()
        axs[i].legend(reversed(handles), reversed(labels), title='Flow Rate')    
    axs[i].grid()
row_index+=1

# Plot the number of connected components
for i in range(offset*row_index,offset*(row_index+1)):
    axs[i].plot(mean_array[:,1], n_connected_components_composite[i-offset*row_index])
    axs[i].set_xlim(0)
#    axs[i].set_ylim([0,300])
    axs[i].set_xlabel('kills (vessels)')    
    if i==len(mean_list*row_index):
        axs[i].set_ylabel('$β_0$') 
    axs[i].grid()
row_index+=1
    
# Plot the number of loops
for i in range(offset*row_index,offset*(row_index+1)):
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
for i in range(offset*row_index,offset*(row_index+1)):
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
for i in range(offset*row_index,offset*(row_index+1)):
    
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
    
# Show plots
plt.show()
'''
