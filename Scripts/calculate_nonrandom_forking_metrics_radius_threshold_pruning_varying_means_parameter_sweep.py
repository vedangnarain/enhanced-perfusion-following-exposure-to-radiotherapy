#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 18:24:17 2023

@author: narain

Tested in Python 3.7.4.

Generate a 3x3 figure for a new set of parameters for the forking network with 
diameter threshold pruning and varying means.
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

# =============================================================================
# FUNCTIONS
# =============================================================================

# Define a function to read an O2 distribution .vti file and return the basic stats
def get_distribution_stats(solver_name, alpha_value, mean_value, radius_threshold, hypoxic_threshold_list, diameter_threshold_list, flow_rate_threshold_list, plot=0, read=0):    
    
    # Set the file path
    folder_path = main_folder_path + solver_name + 'Haematocrit/Lambda4/Alpha' + alpha_value + '/Mean' + mean_value + '/RadiusThreshold' + radius_threshold
    # field_path =  folder_path + '/oxygen_solution_0.vti'    
    network_path = folder_path + '/FinalHaematocrit.vtp'
    
    # Print status update
    # print(field_path)

    # Save file or read from existing file
    if read==0:    
        
        # Compute the architectural metrics
        n_vessels, n_perfused_vessels, n_unperfused_vessels, mean_diameter, mean_geometric_resistance, diameter_binned, flow_rate_binned, n_connected_components, n_cycles, size_largest_connected_component, n_cycles_largest_connected_component = get_forking_predictors(network_path, reference_rank_lengths, reference_node_coordinates, pf_threshold, diameter_threshold_list, flow_rate_threshold_list)

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

# Define a function to return statistics for all the means in the data
def get_solver_stats(solver_name, alpha_value, mean_list, radius_threshold_list, hypoxic_threshold_list, plot, read):
#    table = np.array([])
    alpha_table = np.array([])
    for mean_value in mean_list:
        for radius_threshold in radius_threshold_list:    
            average_radius_threshold_data = get_distribution_stats(solver_name, alpha_value, mean_value, radius_threshold, hypoxic_threshold_list, diameter_threshold_list, flow_rate_threshold_list, plot, read)
            table_entry = np.hstack([alpha_value, int(mean_value), float(radius_threshold), average_radius_threshold_data])
            alpha_table = np.vstack([alpha_table, table_entry]) if alpha_table.size else table_entry
    return alpha_table

# =============================================================================
# SIMULATION CHOICE
# =============================================================================

# Enter details to allow looping over folders
n_gens = 6
d_inlet = 75    
which_solver = 0
pf_threshold = 3.e-12
# which_mean = 2
# which_alpha = 2
max_radius_threshold = 25
# step_size = 1.00
step_size = 0.1

# Set folder path
main_folder_path = '/home/narain/Desktop/Results/Forking/Nonrandom/Without Oxygen/Radius Threshold Pruning in Nonrandom Forking Network with Varying Means (' + str(n_gens) + ' Gens ' + str(d_inlet) + ' Inlet)/TestDichotomousNetwork/'
# main_folder_path = '/home/narain/Desktop/Results/Forking/Nonrandom/Without Oxygen/Individual Pruning in Nonrandom Forking Network with Varying Means (' + str(n_gens) + ' Gens ' + str(d_inlet) + ' Inlet) and Higher Haematocrit/TestDichotomousNetwork/'
# main_folder_path = '/home/narain/Desktop/Results/Forking/Nonrandom/Without Oxygen/Individual Pruning in Nonrandom Forking Network with Varying Means (' + str(n_gens) + ' Gens ' + str(d_inlet) + ' Inlet) and Lower Haematocrit/TestDichotomousNetwork/'
#main_folder_path = '/tmp/narain/testoutput/TestDichotomousNetwork/'
# main_folder_path = '/tmp/narain/testoutput/TestDichotomousNetworkHigherH/'

# Set file name for saving plots
description = 'forking_radius_threshold_pruning_' + 'parameter_sweep_'

# =============================================================================
# SIMULATION PARAMETERS
# =============================================================================

# Set simulation specifics
solver_list = ['Constant', 'Pries', 'Memory', 'Fung']
# alpha_list = ['1.00', '1.10', '1.20', '1.30', '1.40']
alpha_list = ['1.10', '1.20', '1.30']
# graph_alpha_list = ['1.0', '1.1', '1.2', '1.3', '1.4']  # what we want displayed on the graph
graph_alpha_list = ['1.1', '1.2', '1.3']  # what we want displayed on the graph
mean_list = ['0', '1', '2']
graph_mean_list = ['22.76 μm', '28.50 μm', '33.65 μm']  # what we want displayed on the graph
# radius_threshold_list = [str(x) for x in range(max_kills + 1)]
radius_threshold_list = [str(format(x, '.2f')) for x in np.arange(0, max_radius_threshold+step_size, step_size)]
hypoxic_threshold_list = [2195, 10000, 15000, 20000, 25000, 27441] 
#hypoxic_threshold_list = [2195, 5488, 10976, 16465, 21953, 27441] 
#hypoxic_threshold_list = [1,2,3,4] 
diameter_threshold_list = [19, 27, 36]  # in um
flow_rate_threshold_list = [1e-12, 2e-12, 8e-12]
#hypoxic_threshold_list_pp = [0.8, 2, 4, 6, 8, 10] 
solver_name = solver_list[which_solver]

# Get reference forking network (all node coordinates are the same, regardless of heterogeneity)
generation_coordinates, reference_rank_lengths, reference_node_coordinates = get_reference_forking_network(generations=n_gens, inlet_diameter=d_inlet)

# =============================================================================
# DISTRIBUTION STATS & HYPOXIC FRACTIONS
# =============================================================================
'''
# Get the stats for all solvers, varying the value of alpha (change to read=1 to extract from .vti files directly)
solver_stats_1 = get_solver_stats(solver_name, alpha_list[0], mean_list, radius_threshold_list, hypoxic_threshold_list, plot=0, read=0)
solver_stats_2 = get_solver_stats(solver_name, alpha_list[1], mean_list, radius_threshold_list, hypoxic_threshold_list, plot=0, read=0)
solver_stats_3 = get_solver_stats(solver_name, alpha_list[2], mean_list, radius_threshold_list, hypoxic_threshold_list, plot=0, read=0)

# Save array
#np.save(main_folder_path + solver_name + 'Haematocrit/python_solver_data.npy', solver_stats)
#solver_stats = np.load(main_folder_path + solver_name + 'Haematocrit/python_solver_data.npy')

# Combine all the data into one dataframe
# mean_array, hypoxic_fraction_composite, mean_composite, min_composite, half_composite, max_composite, sd_composite, n_vessels_composite, n_perfused_vessels_composite, n_unperfused_vessels_composite, mean_diameter_composite, mean_geometric_resistance_composite, diameter_binned_composite, flow_rate_binned_composite, n_connected_components_composite, n_cycles_composite, size_largest_connected_component_composite, n_cycles_largest_connected_component_composite = filter_by_alpha(mean_list, solver_stats)
# mean_array, _, _, _, _, _, _, n_vessels_composite, n_perfused_vessels_composite, n_unperfused_vessels_composite, mean_diameter_composite, mean_geometric_resistance_composite, diameter_binned_composite, _, n_connected_components_composite, n_cycles_composite, size_largest_connected_component_composite, n_cycles_largest_connected_component_composite = filter_by_alpha(mean_list, solver_stats_1)
complete_alpha_df_1 = pd.DataFrame(solver_stats_1, columns=['alpha', 'mean', 'radius_threshold', '_', '_', '_', '_', '_', '_', 'n_vessels_composite', 'n_perfused_vessels_composite', 'n_unperfused_vessels_composite', 'mean_diameter_composite', 'mean_geometric_resistance_composite', '_', '_', '_', '_', '_', '_'])
# relevant_alpha_df_1 = complete_alpha_df_1.drop(columns='_')
complete_alpha_df_2 = pd.DataFrame(solver_stats_2, columns=['alpha', 'mean', 'radius_threshold', '_', '_', '_', '_', '_', '_', 'n_vessels_composite', 'n_perfused_vessels_composite', 'n_unperfused_vessels_composite', 'mean_diameter_composite', 'mean_geometric_resistance_composite', '_', '_', '_', '_', '_', '_'])
# relevant_alpha_df_2 = complete_alpha_df_2.drop(columns='_')
complete_alpha_df_3 = pd.DataFrame(solver_stats_3, columns=['alpha', 'mean', 'radius_threshold', '_', '_', '_', '_', '_', '_', 'n_vessels_composite', 'n_perfused_vessels_composite', 'n_unperfused_vessels_composite', 'mean_diameter_composite', 'mean_geometric_resistance_composite', '_', '_', '_', '_', '_', '_'])
frames = [complete_alpha_df_1, complete_alpha_df_2, complete_alpha_df_3]
extra_big_df = pd.concat(frames)

# Remove irrelevant information
big_df = extra_big_df.drop(columns='_')

# =============================================================================
# PERFUSION QUOTIENTS
# =============================================================================

# Read PF file
filename = main_folder_path + 'forking_nonrandom_threshold_pruning_perfusion_quotients.txt'
pf_df = pd.read_csv(filename, delim_whitespace=True, names=["network_name", "solver_name", "lambda", "alpha", "mean", "radius_threshold", "PF"])
#pf_df = pd.read_csv(filename, delim_whitespace=True, names=["alpha", "beta", "PF"], skiprows=1)

# Filter PF data for multiple solvers and alphas
solver_filter = solver_name + 'Haematocrit/'
pf_df = pf_df.loc[(pf_df["solver_name"] == solver_filter)]

# Drop extra radius_threshold data
pf_df = pf_df.loc[(pf_df["radius_threshold"] <= max_radius_threshold)]

# Merge with big_df
merger_pf_df = pf_df.drop(columns=["network_name", "solver_name", "lambda"])  # remove irrelevant data 
merger_pf_df['alpha'] = merger_pf_df['alpha'].apply(lambda x: '{:.2f}'.format(float(x)).zfill(4))  # convert the PF data's alpha column to the right format
combined_df = merger_pf_df.merge(big_df, how = 'inner', on = ['alpha', 'mean', 'radius_threshold'])  # merge the dataframes

# Rectify flow rate error by using values at second radius threshold for first
for i in range(len(combined_df)-1):
    if (combined_df.at[i, 'radius_threshold'] == 0):
        combined_df.at[i, 'PF'] = combined_df.at[i+1, 'PF']
        combined_df.at[i, 'n_perfused_vessels_composite'] = combined_df.at[i+1, 'n_perfused_vessels_composite']
        combined_df.at[i, 'n_unperfused_vessels_composite'] = combined_df.at[i+1, 'n_unperfused_vessels_composite']
'''

# Save/read dataframe for quicker execution
# combined_df.to_csv(main_folder_path + solver_name + 'Haematocrit/combined_df.csv', index=False)  # index=False prevents saving the index as a separate column
combined_df = pd.read_csv(main_folder_path + solver_name + 'Haematocrit/combined_df.csv')

# =============================================================================
# PLOT 3x3 FIGURE
# =============================================================================

# Set the figure design
fig, axs = plt.subplots(3, 3, figsize=(21, 18), tight_layout = {'pad': 2})
fig.subplots_adjust(hspace = .75, wspace = 0.5)
linestyles = ['solid', 'dashed', 'dotted', 'dashdot', (0,(5,10)), (0, (3, 1, 1, 1, 1, 1))]
linecolours = ['#0072B2', '#009E73','#D55E00']
width = 5

# Plot the distribution stats for a solver
# axs = axs.ravel()
offset = len(mean_list)
row_index = 0

# Plot the number of vessels
for row in range(len(alpha_list)):
    for col in range(len(mean_list)): 
        graph_df = combined_df
        # graph_df = graph_df.loc[(graph_df["mean"] == mean_list[col]) & (graph_df["alpha"] == alpha_list[row])]        
        ax = axs[row, col]
        graph_df = graph_df.loc[(graph_df["mean"] == int(mean_list[col])) & (graph_df["alpha"] == float(alpha_list[row]))]        
        line1, = ax.plot(graph_df['radius_threshold']*2, graph_df['n_perfused_vessels_composite'], ls='dashed', label='perfused', lw=width, c=linecolours[1])
        line2, = ax.plot(graph_df['radius_threshold']*2, graph_df['n_unperfused_vessels_composite'], ls='dotted', label='hypoperfused', lw=width, c=linecolours[2])
        ax.set_xlim(0, max_radius_threshold*2)
        ax.set_ylim(0, 250)
        ax.set_xlabel('dosage (diameter threshold in μm)')    
        ax.grid()
        ax2 = ax.twinx()        
        relative_pf = (graph_df['PF'] - graph_df['PF'].iloc[0])*100/graph_df['PF'].iloc[0]
        line3, = ax2.plot(graph_df['radius_threshold']*2, relative_pf, label='$\Delta_{\%} \mathcal{P}$', lw=width, c=linecolours[0])
        ax2.set_ylim(-100,250)  
        # ax2.set_ylim(-100,300)  
        # fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)
        if row==0:
            ax.title.set_text(r'$\overline{d}$ = ' + graph_mean_list[col])  
            if col==0:
                lines = [line1, line2, line3]
                labels = [line.get_label() for line in lines]
                ax.legend(lines, labels, loc='upper right')
        if col == 0:
            ax.text(-0.5, 0.5, '${α}$ = ' + graph_alpha_list[row], va='center', ha='right', transform=axs[row, col].transAxes)
            ax.set_ylabel('composition (number of vessels)') 
        if col == 2:
            ax2.set_ylabel(r'$\Delta_{\%} \mathcal{P}$')  #, c='C0
                  
# =============================================================================
# SAVE IMAGE
# =============================================================================

# Show plots
plt.show()

# Save image
file_path = Path('~/Desktop/Final Figures/' + description + '.tif').expanduser()
fig.savefig(file_path, dpi=300, bbox_inches = 'tight', pil_kwargs={"compression": "tiff_lzw"})
# file_path = Path('~/Desktop/Final Figures/' + description + '.svg').expanduser()
# fig.savefig(file_path, dpi=250, bbox_inches = 'tight')

# Prints execution time
print("\n--- Execution Time: %s seconds ---" % (time.time() - start_time))