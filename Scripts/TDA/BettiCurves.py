#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 14:44:32 2022

@author: Bernadette

Function to get Betti curves based on Moo Chung's Matlab code which we previously modified

Input:  Weighted_adjacency_matrix - Adjacency matrix weighted by either lengths or diameters
        thresholds - increasing sequence of pruning thresholds

Output: beta0 - Betti 0 curve
        beta1 - Betti 1 curve
        size_largest_connected_component - number of nodes in largest connected component
        biggest1 - Betti1 curve for largest connected component (note that the largest connnected component can consist of completely different sets of nodes as we change the threshold)
"""

import numpy as np
import networkx as nx

def PH_betti2(Weighted_adjacency_matrix, thresholds):
    
    beta0=np.array([])
    beta1=np.array([])
    size_largest_connected_component = np.array([])
    biggest1 = np.array([])
    
    for threshold in thresholds:
        
        print('Threshold = ',threshold)
        adj = np.array(Weighted_adjacency_matrix)
        adj[adj <= threshold] = 0
        
        # we don't want any nodes showing up as separate connected components
        adj = adj[~np.all(adj == 0, axis=1), :] #rows
        adj = adj[:,~np.all(adj == 0, axis=0)] #cols - logically axis should be =1 here, but code only works if I set it to 0

        n_nodes = adj.shape[0]
        
        G = nx.Graph(adj)
        number_of_connected_components = nx.number_connected_components(G)
        
        beta0=np.append(beta0, number_of_connected_components)
        
        H = max(nx.connected_components(G), key=len)
        H_subgraph = G.subgraph(H).copy()
        size_largest_cc = H_subgraph.number_of_nodes()
        size_largest_connected_component = np.append(size_largest_connected_component, size_largest_cc)

        #computes Betti-1

        n_edges = G.number_of_edges()
        n_cycles = number_of_connected_components - n_nodes + n_edges
        beta1=np.append(beta1, n_cycles)
        
        
        #computes Beti-1 of largest connected component

        n_cycles_largest_connected_component = 1 - size_largest_cc + H_subgraph.number_of_edges()
        biggest1 = np.append(biggest1,n_cycles_largest_connected_component)
    
    
    return beta0, beta1, size_largest_connected_component, biggest1
