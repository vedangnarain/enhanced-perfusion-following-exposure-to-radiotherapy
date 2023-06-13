#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Fri Feb 18 01:54:42 2022

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)

Four Generation Network Analysis

This network plots the haematocrit distribution from a Chaste simulation to
compare it with the microchannel experiment of Hyakutake et al., Microvascular 
Research 140 (2022) 104281.

Tested in Python 3.7.4.

"""

# =============================================================================
# LIBRARIES & INITIALISATION
# =============================================================================

# Initialise libraries
import matplotlib.pyplot as plt
import pandas as pd
import time

# Starts stopwatch to clock execution time
start_time = time.time()

# =============================================================================
# DATA
# =============================================================================

# Import the data (pruned vessels numbered by quadrant)
data = pd.read_csv('/home/narain/Desktop/four_generation_network.csv')
haematocrit_series = data['Vessel Haematocrit']
haematocrit = haematocrit_series.values[[60,0,8,9,24,25,20,21,56,57,52,53,48,49,44,45]]  # arrange the vessel values in line with the experimental paper

# =============================================================================
# PLOTS
# =============================================================================

# Plot unpruned comparative histogram
plt.figure()
plt.style.use('fivethirtyeight')
plt.title('Tube Haematocrits')
tubes = ['P', 'A1', 'B1', 'B2', 'C1', 'C2', 'C3', 'C4', 'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8']
plt.bar(tubes,haematocrit)
plt.ylabel('H$_T$')
plt.ylim(0.21)
plt.legend()
plt.show()

# prints execution time
print("\n--- Execution Time: %s seconds ---" % (time.time() - start_time))
