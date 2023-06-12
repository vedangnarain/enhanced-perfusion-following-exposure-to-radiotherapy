#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 12:35:22 2022

@author: narain

Tested in Python 3.7.4.

Plot the PQs for a network with a varying threshold. Varying the thresholds helps pick a threshold that ensures that the initial PQ of a network varies enough.

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

# Define a function to calculate the perfusion quotient for a network with different thresholds
#def get_pq(vtk_path):    
    
#vtk_path = '/tmp/narain/testoutput/TestHexagonalNetwork/ConstantHaematocrit/Mu10/Selection1/Kills0/FinalHaematocrit.vtp'
#vtk_path = '/tmp/narain/testoutput/TestHexagonalNetwork/ConstantHaematocrit/Mu10/Selection1/Kills1/FinalHaematocrit.vtp'
#vtk_path = '/tmp/narain/testoutput/TestHexagonalNetwork/ConstantHaematocrit/Mu10/Selection1/Kills2/FinalHaematocrit.vtp'
vtk_path = '/tmp/narain/testoutput/TestHexagonalNetwork/ConstantHaematocrit/Mu19/Selection13/Kills0/FinalHaematocrit.vtp'

pq_thresholds = [1.e-14, 1.e-15]

point_coordinates_data_array, point_data_array, cell_data_array, polydata = get_vtk_data(vtk_path)
n_vessels = len(cell_data_array['Absolute Vessel Flow Rate m^3/s'])
for threshold in pq_thresholds:
    n_perfused_vessels = len(cell_data_array['Absolute Vessel Flow Rate m^3/s'][cell_data_array['Absolute Vessel Flow Rate m^3/s'] >= threshold])
    pqs = n_perfused_vessels/n_vessels
    print(pqs)
