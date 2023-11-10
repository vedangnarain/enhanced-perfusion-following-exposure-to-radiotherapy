#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Jun  1 18:23:28 2021

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)

Generate truncated normal and log normal distributions for hexagonal vessel radii.

The minimum is set based on the biological tumours. The max. is set by the length of 
the inlets

Tested in Python 3.7.4.
"""

# Initialise libraries
import errno
import matplotlib.pyplot as plt
import numpy as np
import os, os.path
import scipy.stats as stats

# Set LaTex-style font
from pathlib import Path
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 22})

#mu, sigma = 3., 1. # mean and standard deviation
#s = np.random.lognormal(mu, sigma, 1000)
#Display the histogram of the samples, along with the probability density function:
#
#import matplotlib.pyplot as plt
#count, bins, ignored = plt.hist(s, 100, density=True, align='mid')
#x = np.linspace(min(bins), max(bins), 10000)
#pdf = (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))
#       / (x * sigma * np.sqrt(2 * np.pi)))
#plt.plot(x, pdf, linewidth=2, color='r')
#plt.axis('tight')
#plt.show()

# Set the default mean, lower bound, upper bound, and number of vessels (in diameters)
#mean = 13
lower_bound, upper_bound = 7.39, 75

# Pick SD
sd_diameter_list = [8.68, 13.23, 17.49]
sigma = sd_diameter_list[2]

# Specify the list of means
mean_diameter_list = np.array([22.76, 28.5, 33.64])

#n_vessels = 407
#n_vessels = 386
#n_inlets = 0

#'''
n_inlets = 11
#n_outlets = n_inlets-1
n_vessels = 386
#n_vessels = 407-n_outlets
#'''

# Specify the list of heterogeneities (standard deviations)
#sigma_list = [8, 10, 13, 15]

# Convert diameters to radii
lower_bound, upper_bound = lower_bound/2, upper_bound/2
sigma = sigma/2
mean_list = mean_diameter_list/2

# Specify network name for 
file_name = 'hexagonal_diameter_log_normal_distribution_sigma_'

# =============================================================================
# FUNCTIONS
# =============================================================================

# Define a function to make a directory if needed
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

# Define a function to safely open a directory to write a file
def safe_open_w(path):
    mkdir_p(os.path.dirname(path))
    return open(path, 'w')
    
# Define a function to add vessel IDs to the radii list    
def add_vessel_ids(radii_list, samples=n_vessels):
    
    # Generate vessel IDs
    vessel_ids = [x for x in range(n_inlets, samples+n_inlets)]
    
    # Add IDs to radii list
    radii_id_list = np.matrix.transpose(np.vstack([vessel_ids, radii_list]))
    
    return(radii_id_list)
    
# Define a function to generate a log normal distribution
def generate_log_normal_distribution(mu, lower=lower_bound, upper=upper_bound, sigma=sigma, samples=n_vessels, add_ids=False):
        
    # Convert the desired mean to log normal mean
    log_mu = (2*np.log(mu)) - (0.5*np.log((sigma**2)+(mu**2)))
#        print(log_mu)
    
    # Convert the desired SD to log normal SD
    log_sigma = np.sqrt((-2*np.log(mu)) + np.log((sigma**2)+(mu**2)))
#        print(log_sigma)

    # Draw samples for all the non-inlet/outlet vessel radii
    radii_list = np.random.lognormal(log_mu, log_sigma, samples)

    # Replace values that exceed the upper and lower bounds
    for i in range(len(radii_list)):
        while radii_list[i]<lower or radii_list[i]>upper:    
                radii_list[i] = np.random.lognormal(log_mu, log_sigma, 1)

    # Check if bounds are exceeded in list
    if min(radii_list)<lower and max(radii_list)>upper:
            print('ERROR: BOUNDS EXCEEDED')

    # Round off values to five decimal places
    radii_list = np.around(radii_list, decimals=5)

    # Add vessel IDs
    if add_ids==True:
        radii_list = add_vessel_ids(radii_list)
    
    # Return the list of radii
    return(radii_list)
    
#'''
# Set the figure layout
fig, axs = plt.subplots(nrows=1, ncols=3, sharex=True,  figsize=(15, 5), tight_layout = {'pad': 2.5})
fig.subplots_adjust(hspace = 0.75, wspace=.25)
#plt.suptitle('Hexagonal Network Vessel Radius Distribution (μ = ' + str(mean) + ' μm, n_vessels = ' + str(n_vessels) + ', min = ' + str(lower_bound) + ' μm, max = ' + str(upper_bound) + ' μm)')

# Plot the stats
axs = axs.ravel()
for i in range(len(mean_list)):
    distribution = generate_log_normal_distribution(mean_list[i], add_ids=False)
    print(np.mean(distribution))
    # write to file here (beta_i/selection_y). Actually make plot function inside and generate outside.
    n, bins, patches = axs[i].hist(x=distribution*2, bins='auto', alpha=0.7)
    if i==0:
        axs[i].set_ylabel('number of vessels') 
#    if i==2 or i==3:
    axs[i].set_xlabel('vessel diameter (μm)')    
#    print(np.mean(distribution))
#    print(np.min(distribution), np.max(distribution))
#    axs[i].legend()
#    axs[i].grid()
    axs[i].title.set_text(r'$\overline{d}$ =' + str(mean_diameter_list[i]) + ' μm')
    axs[i].tick_params(labelbottom=True)
    axs[i].axvline(mean_diameter_list[i], c='black', ls='--', label='mean')
    axs[i].set_ylim(0, 100)
    if i==0:
        axs[i].legend()
plt.show()
#'''

#'''
# Write 100 radii lists for each mean value
form = "%.5f \n"
for mu_number in mean_list:
    for list_number in range(1,101):
        id_name = '/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/100VesselLength/' + file_name + str(sigma*2) + '/mu_' + str(mu_number*2) + '/id_list_' + str(list_number) + '.txt'
        radii_name = '/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/100VesselLength/' + file_name + str(sigma*2) + '/mu_' + str(mu_number*2) + '/radii_list_' + str(list_number) + '.txt'
        composite_list = generate_log_normal_distribution(mu_number, add_ids=True)  
        sorted_composite_list = composite_list[np.argsort(composite_list[:, 1])]
        sorted_id_list = sorted_composite_list[:,0]
        sorted_radii_list = sorted_composite_list[:,1]
        with safe_open_w(id_name) as i_f:
#            for i in range(len(radii_list)):
#                vector = 
            np.savetxt(id_name, sorted_id_list, fmt='%i')
        with safe_open_w(radii_name) as r_f:
#            for i in range(len(radii_list)):
#                vector = 
            np.savetxt(radii_name, sorted_radii_list, fmt='%1.5f')
            
# Save image
file_path = Path('/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/100VesselLength/' + file_name + str(sigma*2) + '/hexagonal_diameter_log_normal_distribution.svg').expanduser()
fig.savefig(file_path, dpi=500)
file_path = Path('/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/100VesselLength/' + file_name + str(sigma*2) + '/hexagonal_diameter_log_normal_distribution.png').expanduser()
fig.savefig(file_path, dpi=500)
#'''
                    