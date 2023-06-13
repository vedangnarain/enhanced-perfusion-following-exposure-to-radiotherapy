#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 14:43:53 2022

@author: Bernadette 

Modified by Vedang.

This script takes in an adjacency matrix for an unpruned network and returns
the Betti curves.
"""

# =============================================================================
# LIBRARIES & INITIALISATION
# =============================================================================

# Initialise libraries
import numpy as np

# Import tools for TDA
from BettiCurves import PH_betti2


# =============================================================================
# BETTI CURVES
# =============================================================================

# Load the adjacency matrix
file_path = '/Users/user/Desktop/Vasculature Project/Project with Jakub/'
#adjacency_matrix = sio.loadmat(file_path+'AdjacencyMatrixDiameter.mat', squeeze_me=True)
adjacency_matrix = np.loadtxt(file_path+'TestMatrix.txt', delimiter=" ", dtype="float")

# Convert the diameters into radii
A = adjacency_matrix/2 

# Generate the pruning steps
max_radius= np.max(A)
min_radius = np.min(A)
thresholds= np.arange(min_radius, max_radius, 1)

# Get the beta0_radius, beta1_radius, biggest0, biggest1
#beta0_radius, beta1_radius, biggest0, biggest1 = PH_betti2(A, thresholds)
n_connected_components, n_loops, largest_connected_component, n_loops_in_largest_component = PH_betti2(A, thresholds)
