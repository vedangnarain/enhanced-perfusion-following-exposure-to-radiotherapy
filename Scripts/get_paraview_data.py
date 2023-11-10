#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Mar 31 16:44:36 2022

@author: narain

Tested in Python 3.7.4.

Read VTK and VTI files to process them in Python.

"""

# =============================================================================
# LIBRARIES & INITIALISATION
# =============================================================================

# Initialise libraries
from collections import deque
#from vtk import *
from vtk.util.numpy_support import vtk_to_numpy
import math
import networkx as nx
import numpy as np
import pandas as pd
import scipy.interpolate
import vtk

# Set up plotting
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 22})

# =============================================================================
# FUNCTIONS
# =============================================================================

# Define a function to read a .vtk file and return the point coordinates and data and the cell data
def get_vtk_data(vtk_path):
    print(vtk_path)

    # Set up the file reader
#    if 'vtp' in filename:
#        reader = vtk.vtkXMLPolyDataReader()
#    else:
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtk_path)
    reader.Update()
    polydata = reader.GetOutput()
#    polydata.ReleaseDataFlagOn()
    
    # Extract the point coordinates
    point_coordinates_data = polydata.GetPoints().GetData()
    point_coordinates_data_array = vtk_to_numpy(point_coordinates_data)
    
    # Extract the point data 
    point_data_array = {}
    point_data = polydata.GetPointData()
    for column_index in range(point_data.GetNumberOfArrays()):
       column_name =  point_data.GetArrayName(column_index)
       point_data_array[column_name] = vtk_to_numpy(point_data.GetArray(column_index))
    
    # Extract the cell data    
    cell_data_array = {}
    cell_data = polydata.GetCellData()
    for column_index in range(cell_data.GetNumberOfArrays()):
       column_name =  cell_data.GetArrayName(column_index)
       cell_data_array[column_name] = vtk_to_numpy(cell_data.GetArray(column_index))
    
    # Return a dictionary with the point coordinates and data and cell data
    return point_coordinates_data_array, point_data_array, cell_data_array, polydata
           
# Define a function to read a .vti file and return the data
def get_vti_data(vti_path):
    
    # Import the file
    extension = vti_path.split('.').pop()
    reader = None
    if extension == 'vtk':
        reader = vtk.vtkDataSetReader() 
    elif extension == 'vti':
        reader = vtk.vtkXMLImageDataReader() 
    else:
        raise RuntimeError('Unknown File Type: %s ' % vti_path)
    reader.SetFileName( "%s" % (vti_path) ) 
    reader.Update()
    image_data = reader.GetOutput()
    
    # Extract the dimensions, spacing, and origin
    spacing = image_data.GetSpacing()
    
    # Extract the point values 
    field_point_data = image_data.GetPointData() 
    field_values = vtk_to_numpy(field_point_data.GetArray(0)) 
    
    # Get the coordinates of each point
    position_list = deque()
    for index in range(len(field_values)):  # do something 
        position = image_data.GetPoint(index)
        position_list.append(position)
    position_array = np.array(position_list)
    
    # Return the field distribution
    distribution_array = np.column_stack((position_array,field_values))
    distribution = pd.DataFrame(data=distribution_array[:,[0,1,3]], columns=["x", "y", "oxygen"])
    return distribution, spacing#, dimensions, origin       

# Define a function to convert the .vti data to a plottable field with some basic stats
def get_plottable_field(vti_data):
    
    # Change the data type
    vti_data = vti_data.astype('float32')
    
    # Calculate concentration statistics
    O2_stats = vti_data['oxygen'].describe()
    
    # Downsample data if needed
    vti_data = vti_data[::int(1)]
    
    # Convert dataframe into NumPy matrix
    mat = vti_data.to_numpy()
    
    # Get the x and y axes
    x = np.unique(mat[:,0])
    y = np.unique(mat[:,1])
    
    # Create a mesh from the x and y axes 
    X,Y = np.meshgrid(x, y)
    
    # Interpolate the concentration values over the mesh
    Z = scipy.interpolate.griddata((mat[:,0], mat[:,1]), mat[:,2], (X,Y), method='nearest')

    # Return the oxygen stats
    return Z, O2_stats

# Define a function to extract the region of evaluation of the forking network
def get_forking_domain(field_dataset, middle_generation_number, generation_coordinates):
    x_start = generation_coordinates[middle_generation_number-1]
    x_end = generation_coordinates[-middle_generation_number]
    middle_field_dataset = field_dataset[(field_dataset['x'] >= int(x_start)) & (field_dataset['x'] <= int(x_end))]
    return middle_field_dataset

# Define a function to extract the the region of evaluation of the hexagonal network
def get_hex_domain(field_dataset, field_spacing):
    x_start = 10*field_spacing[0]
    x_end = 195*field_spacing[0]
    y_start = 17*field_spacing[1]
    y_end = 173*field_spacing[1]
    field_dataset = field_dataset[(field_dataset['x'] >= int(x_start)) & (field_dataset['x'] <= int(x_end)) & (field_dataset['y'] >= int(y_start)) & (field_dataset['y'] <= int(y_end))]
    return field_dataset
#    return vessel_network, oxygen_distribution, middle_x

# Define a function to extract the predictive metrics from the forking network
def get_forking_predictors(vtk_path, reference_rank_lengths, reference_node_coordinates, pq_threshold, diameter_threshold_list, flow_rate_threshold_list):
        
    # Get the .vtk data
    point_coordinates_data_array, point_data_array, cell_data_array, polydata = get_vtk_data(vtk_path)
    
    # Convert the radii and lengths in metres into the right units
    segment_diameters_um = cell_data_array['Vessel Radius m']*2*1000000  # convert radii (m) to diameters (um)
    reference_rank_lengths_um = reference_rank_lengths*1000000  # convert length (m to um)
    
    # Get the architectural metrics
    n_vessels = len(cell_data_array['Vessel Radius m'])
    n_perfused_vessels = len(cell_data_array['Absolute Vessel Flow Rate m^3/s'][cell_data_array['Absolute Vessel Flow Rate m^3/s'] >= pq_threshold])
    n_unperfused_vessels = n_vessels-n_perfused_vessels
    mean_diameter = sum(segment_diameters_um)/n_vessels
    segment_lengths = [reference_rank_lengths_um[int(vessel_rank)] for vessel_rank in cell_data_array['Vessel Owner Rank']]  # get vessel length based on rank    
    mean_geometric_resistance = sum(segment_lengths/(segment_diameters_um**4))/n_vessels
    diameter_binned = filter_by_thresholds(segment_diameters_um, diameter_threshold_list)
    flow_rate_binned = filter_by_thresholds(cell_data_array['Absolute Vessel Flow Rate m^3/s'], flow_rate_threshold_list)

    # Get a list of vessel segments with node IDs for the adjacency matrices
    cellIds = vtk.vtkIdList()  # cell ids store to
    numberOfCells = polydata.GetNumberOfCells()
    segment_nodes = np.array([])
    for cellIndex in range(numberOfCells):  # for every cell
    #    print('new cell')
        polydata.GetCellPoints(cellIndex, cellIds)  # get IDs of nodes of the given cell
        cell_nodes = np.array([])
        for i in range(0, cellIds.GetNumberOfIds()):  # for every node of the given cell
            coord = polydata.GetPoint(cellIds.GetId(i))  # get coordinates of the node, type: class 'tuple'
            x = np.around(coord[0], 2)  # get x-coordinate of the node, type: class 'float'
    #        print(x)
            y = np.around(coord[1], 2)  # get y-coordinate of the node, type: class 'float'
    #        print(y)
            node_id = np.where((reference_node_coordinates[:,0] == x) & (reference_node_coordinates[:,1] == y))[0]
    #        print(node_id)
            cell_nodes = np.hstack([cell_nodes, node_id])
        segment_nodes = np.vstack([segment_nodes, cell_nodes]) if segment_nodes.size else cell_nodes
    segment_nodes = segment_nodes.astype('int')
    
    # Get the number of nodes
    number_of_nodes = len(reference_node_coordinates)
    
    # Initialise two nxn matrices of zeros, where n is the number of nodes
    diameter_adjacency_matrix = np.zeros((number_of_nodes, number_of_nodes))
    length_adjacency_matrix = np.zeros((number_of_nodes, number_of_nodes))

    # Fill in the adjacency matrices
    for segment_index, segment in enumerate(segment_nodes):
        
        # Get the segment nodes
        row, col = segment
        
        # Fill in the adjacency matrices
        diameter_adjacency_matrix[row,col] = segment_diameters_um[segment_index]
        diameter_adjacency_matrix[col,row] = segment_diameters_um[segment_index]
        length_adjacency_matrix[row,col] = reference_rank_lengths_um[int(cell_data_array['Vessel Owner Rank'][segment_index])]
        length_adjacency_matrix[col,row] = reference_rank_lengths_um[int(cell_data_array['Vessel Owner Rank'][segment_index])]

    # Get the diameter TDA characteristics
    n_connected_components, n_cycles, size_largest_connected_component, n_cycles_largest_connected_component = PH_betti2(diameter_adjacency_matrix)

    # Return the architectural features and the adjacency matrices
    return n_vessels, n_perfused_vessels, n_unperfused_vessels, mean_diameter, mean_geometric_resistance, diameter_binned, flow_rate_binned, n_connected_components, n_cycles, size_largest_connected_component, n_cycles_largest_connected_component

# Define a function to extract the predictive metrics from the hexagonal network
def get_hex_predictors(vtk_path, vessel_length_m, reference_node_coordinates, inlet_radius_m, pq_threshold, diameter_threshold_list, flow_rate_threshold_list):
    
    # Get the .vtk data
    point_coordinates_data_array, point_data_array, cell_data_array, polydata = get_vtk_data(vtk_path)
    
    # Convert the radii and lengths in metres into the right units
    segment_diameters_um = cell_data_array['Vessel Radius m']*2*1000000  # convert radii (m) to diameters (um)
#    non_inlet_segment_radii_m = cell_data_array['Vessel Radius m'][cell_data_array['Vessel Radius m'] != inlet_radius_m]
#    non_inlet_segment_diameters_um = cell_data_array['Vessel Radius m'][cell_data_array['Vessel Radius m'] != inlet_radius_m]*2*1000000  # convert radii (m) to diameters (um) and disregard inlet/outlet vessels
    segment_diameters_um = cell_data_array['Vessel Radius m']*2*1000000  # convert radii (m) to diameters (um)
    segment_length_um = vessel_length_m*1000000  # convert length (m to um)
       
    # Get the architectural metrics    
    n_vessels = len(segment_diameters_um)
    n_perfused_vessels = len(cell_data_array['Absolute Vessel Flow Rate m^3/s'][cell_data_array['Absolute Vessel Flow Rate m^3/s'] >= pq_threshold])
    n_unperfused_vessels = n_vessels-n_perfused_vessels
    mean_diameter = sum(segment_diameters_um)/n_vessels
    mean_geometric_resistance = sum(segment_length_um/(segment_diameters_um**4))/n_vessels
    diameter_binned = filter_by_thresholds(segment_diameters_um, diameter_threshold_list)
    flow_rate_binned = filter_by_thresholds(cell_data_array['Absolute Vessel Flow Rate m^3/s'], flow_rate_threshold_list)
    
    # Get a list of vessel segments with node IDs for the adjacency matrices
    cellIds = vtk.vtkIdList()  # cell ids store to
    numberOfCells = polydata.GetNumberOfCells()
    segment_nodes = np.array([])
    for cellIndex in range(numberOfCells):  # for every cell
    #    print('new cell')
        polydata.GetCellPoints(cellIndex, cellIds)  # get IDs of nodes of the given cell
        cell_nodes = np.array([])
        for i in range(0, cellIds.GetNumberOfIds()):  # for every node of the given cell
            coord = polydata.GetPoint(cellIds.GetId(i))  # get coordinates of the node, type: class 'tuple'
            x = np.around(coord[0], 2)  # get x-coordinate of the node, type: class 'float'
    #        print(x)
            y = np.around(coord[1], 2)  # get y-coordinate of the node, type: class 'float'
    #        print(y)
            node_id = np.where((reference_node_coordinates[:,0] == x) & (reference_node_coordinates[:,1] == y))[0]
    #        print(node_id)
            cell_nodes = np.hstack([cell_nodes, node_id])
        segment_nodes = np.vstack([segment_nodes, cell_nodes]) if segment_nodes.size else cell_nodes
    segment_nodes = segment_nodes.astype('int')
    
    # Get the number of nodes
    number_of_nodes = len(reference_node_coordinates)
    
    # Initialise two nxn matrices of zeros, where n is the number of nodes
    diameter_adjacency_matrix = np.zeros((number_of_nodes, number_of_nodes))
    length_adjacency_matrix = np.zeros((number_of_nodes, number_of_nodes))
    
    # Fill in the adjacency matrices
    for segment_index, segment in enumerate(segment_nodes):
        
        # Get the segment nodes
        row, col = segment
        
        # Fill in the adjacency matrices
        diameter_adjacency_matrix[row,col] = segment_diameters_um[segment_index]
        diameter_adjacency_matrix[col,row] = segment_diameters_um[segment_index]
        length_adjacency_matrix[row,col] = segment_length_um
        length_adjacency_matrix[col,row] = segment_length_um
    
    # Get the diameter TDA characteristics
    n_connected_components, n_cycles, size_largest_connected_component, n_cycles_largest_connected_component = PH_betti2(diameter_adjacency_matrix)
 
    return n_vessels, n_perfused_vessels, n_unperfused_vessels, mean_diameter, mean_geometric_resistance, diameter_binned, flow_rate_binned, n_connected_components, n_cycles, size_largest_connected_component, n_cycles_largest_connected_component
    
# Define a function to extract the predictive metrics from the Voronoi network (INLET/OUTLETS STILL EXCLUDED)
def get_voronoi_predictors(vtk_path, reference_node_coordinates, inlet_radius_m, pq_threshold, diameter_threshold_list, flow_rate_threshold_list):
    
    # Get the .vtk data
    point_coordinates_data_array, point_data_array, cell_data_array, polydata = get_vtk_data(vtk_path)
    
    # Convert the radii and lengths in metres into the right units
    segment_diameters_m = cell_data_array['Vessel Radius m']*2  # convert radii (m) to diameters (um)
    segment_diameters_um = segment_diameters_m*1000000  # convert radii (m) to diameters (um)
    segment_viscosities = cell_data_array['Vessel Viscosity Pa.s'] 
    segment_impedances = cell_data_array['Vessel Impedance kg/m^4/s']
    segment_lengths_um = ((segment_impedances*math.pi*((segment_diameters_um/2)**4))/(8.0*segment_viscosities))*1000000  # convert length (m to um)

#    non_inlet_segment_radii_m = cell_data_array['Vessel Radius m'][cell_data_array['Vessel Radius m'] != inlet_radius_m]
    non_inlet_segment_viscosities = cell_data_array['Vessel Viscosity Pa.s'][cell_data_array['Vessel Radius m'] != inlet_radius_m]  # disregard inlet/outlet vessels
    non_inlet_segment_diameters_m = cell_data_array['Vessel Radius m'][cell_data_array['Vessel Radius m'] != inlet_radius_m]*2  # disregard inlet/outlet vessels
    non_inlet_segment_diameters_um = non_inlet_segment_diameters_m*1000000  # convert radii (m) to diameters (um)
    non_inlet_segment_impedances = cell_data_array['Vessel Impedance kg/m^4/s'][cell_data_array['Vessel Radius m'] != inlet_radius_m]
    non_inlet_segment_lengths_um = ((non_inlet_segment_impedances*math.pi*((non_inlet_segment_diameters_m/2)**4))/(8.0*non_inlet_segment_viscosities))*1000000  # convert length (m to um)
       
    # Get the architectural metrics    
    n_vessels = len(segment_diameters_um)
    n_perfused_vessels = len(cell_data_array['Absolute Vessel Flow Rate m^3/s'][cell_data_array['Absolute Vessel Flow Rate m^3/s'] >= pq_threshold])
    n_unperfused_vessels = n_vessels-n_perfused_vessels
    mean_diameter = sum(non_inlet_segment_diameters_um)/n_vessels
    mean_geometric_resistance = sum(non_inlet_segment_lengths_um/(non_inlet_segment_diameters_um**4))/n_vessels
    diameter_binned = filter_by_thresholds(segment_diameters_um, diameter_threshold_list)
    flow_rate_binned = filter_by_thresholds(cell_data_array['Absolute Vessel Flow Rate m^3/s'], flow_rate_threshold_list)
    
    # Get a list of vessel segments with node IDs for the adjacency matrices
    cellIds = vtk.vtkIdList()  # cell ids store to
    numberOfCells = polydata.GetNumberOfCells()
    segment_nodes = np.array([])
    for cellIndex in range(numberOfCells):  # for every cell
    #    print('new cell')
        polydata.GetCellPoints(cellIndex, cellIds)  # get IDs of nodes of the given cell
        cell_nodes = np.array([])
        for i in range(0, cellIds.GetNumberOfIds()):  # for every node of the given cell
            coord = polydata.GetPoint(cellIds.GetId(i))  # get coordinates of the node, type: class 'tuple'
            x = np.around(coord[0], 2)  # get x-coordinate of the node, type: class 'float'
    #        print(x)
            y = np.around(coord[1], 2)  # get y-coordinate of the node, type: class 'float'
    #        print(y)
            node_id = np.where((reference_node_coordinates[:,0] == x) & (reference_node_coordinates[:,1] == y))[0]
    #        print(node_id)
            cell_nodes = np.hstack([cell_nodes, node_id])
        segment_nodes = np.vstack([segment_nodes, cell_nodes]) if segment_nodes.size else cell_nodes
    segment_nodes = segment_nodes.astype('int')
    
    # Get the number of nodes
    number_of_nodes = len(reference_node_coordinates)
    
    # Initialise two nxn matrices of zeros, where n is the number of nodes
    diameter_adjacency_matrix = np.zeros((number_of_nodes, number_of_nodes))
    length_adjacency_matrix = np.zeros((number_of_nodes, number_of_nodes))

    # Fill in the adjacency matrices
    for segment_index, segment in enumerate(segment_nodes):
        
        # Get the segment nodes
        row, col = segment
        
        # Fill in the adjacency matrices
        diameter_adjacency_matrix[row,col] = segment_diameters_um[segment_index]
        diameter_adjacency_matrix[col,row] = segment_diameters_um[segment_index]
        length_adjacency_matrix[row,col] = segment_lengths_um[segment_index]
        length_adjacency_matrix[col,row] = segment_lengths_um[segment_index]
    
    # Get the diameter TDA characteristics
    n_connected_components, n_cycles, size_largest_connected_component, n_cycles_largest_connected_component = PH_betti2(diameter_adjacency_matrix)
 
    return n_vessels, n_perfused_vessels, n_unperfused_vessels, mean_diameter, mean_geometric_resistance, diameter_binned, flow_rate_binned, n_connected_components, n_cycles, size_largest_connected_component, n_cycles_largest_connected_component
    
# Import the generation coordinates, rank lengths, and node coordinates for the forking network 
def get_reference_forking_network(generations, inlet_diameter):
    alpha_value = '1.00'  # only the diameter varies with alpha 
#    reference_network_path = '/home/narain/Desktop/Scripts/reference_networks/forking_inlet_diameter_' + str(inlet_diameter) + '_um/Alpha' + alpha_value + '/FinalHaematocrit.vtp'
#    reference_network_path = '/home/narain/Desktop/Scripts/reference_networks/forking_inlet_diameter_50_um_five_generations/FinalHaematocrit.vtp'
    reference_network_path = '/home/narain/Desktop/Scripts/reference_networks/forking_networks/inlet_diameter_' + str(inlet_diameter) + '_um_generations_' + str(generations) + '/Alpha' + alpha_value + '/FinalHaematocrit.vtp'
    reference_point_coordinates_data_array, point_data_array, reference_cell_data_array, polydata = get_vtk_data(reference_network_path)
    generation_coordinates = np.unique(reference_point_coordinates_data_array[:,0])
    rank_diameters = -np.sort(-np.unique(reference_cell_data_array['Vessel Radius m'])*2)
    reference_rank_lengths = rank_diameters*4
    
    # Round the coordinates to two decimal places
    reference_point_coordinates_data_array = reference_point_coordinates_data_array.astype('float64')  # convert to the right data type for later comparison
    reference_node_coordinates = np.around(reference_point_coordinates_data_array[:,0:2], 2)
 
    return generation_coordinates, reference_rank_lengths, reference_node_coordinates

# Import the node coordinates for the hexagonal network 
def get_reference_hexagonal_network(reference_network_path):
#def get_reference_hexagonal_network(vessel_length_m=100*(10**-6)):
#    sigma = 1
#    selection = 1
#    radius_threshold = 0
#    reference_network_path = '/home/narain/Desktop/Scripts/reference_networks/Sigma' + str(sigma) + '/Selection' + str(selection) +'/RadiusThreshold' + str(radius_threshold) + '/FinalHaematocrit.vtp'
    reference_point_coordinates_data_array, _, _, _ = get_vtk_data(reference_network_path)
        
    # Round the coordinates to two decimal places
    reference_point_coordinates_data_array = reference_point_coordinates_data_array.astype('float64')  # convert to the right data type for later comparison
    reference_node_coordinates = np.around(reference_point_coordinates_data_array[:,0:2], 2)
 
    return reference_node_coordinates

'''
# Define a function to filter the simulation stats based on alpha values and return individual matrices (WORK NEEDED)
def filter_by_alpha(alpha_list, solver_stats):
    hypoxic_fraction_composite = np.array([])
    field_composite_collection = np.array([])
    for alpha_value in alpha_list:
        alpha_array = solver_stats[(solver_stats[:,0]==float(alpha_value))]
        hypoxic_fraction_data = alpha_array[:,2:3]  # extract data between identifiers and basic stats (i.e., hypoxic fractions)
        hypoxic_fraction_composite = np.vstack([hypoxic_fraction_composite, hypoxic_fraction_data]) if hypoxic_fraction_composite.size else hypoxic_fraction_data
        for column_index in range(3,15):
            field_data = alpha_array[:,column_index]
            field_composite = np.vstack([field_composite, field_data]) if field_composite.size else field_data
            field_composite_collection = np.hstack([field_composite, field_composite]) if field_composite_collection.size else field_composite
        a
    return alpha_array, hypoxic_fraction_composite, mean_composite, min_composite, half_composite, max_composite, sd_composite, n_vessels_composite, n_perfused_vessels_composite, n_unperfused_vessels_composite, mean_diameter_composite, mean_geometric_resistance_composite, diameter_binned_composite, flow_rate_binned_composite
  '''
        
# Define a function to filter the simulation stats based on alpha values and return individual matrices
def filter_by_alpha(alpha_list, solver_stats):
    hypoxic_fraction_composite = np.array([])
    mean_composite = np.array([])
    min_composite = np.array([])
    half_composite = np.array([])
    max_composite = np.array([])
    sd_composite = np.array([])  
    n_vessels_composite = np.array([])
    n_perfused_vessels_composite = np.array([])
    n_unperfused_vessels_composite = np.array([])
    mean_diameter_composite = np.array([])
    mean_geometric_resistance_composite = np.array([])
    diameter_binned_composite = np.array([])  
    flow_rate_binned_composite = np.array([])
    n_connected_components_composite = np.array([]) 
    n_cycles_composite = np.array([])
    size_largest_connected_component_composite = np.array([])
    n_cycles_largest_connected_component_composite = np.array([])
    for alpha_value in alpha_list:
        alpha_array = solver_stats[(solver_stats[:,0]==float(alpha_value))]
        hypoxic_fraction_data = alpha_array[:,2:3]  # extract data between identifiers and basic stats (i.e., hypoxic fractions)
        mean_data = alpha_array[:,3]
        min_data = alpha_array[:,4]
        half_data = alpha_array[:,5]
        max_data = alpha_array[:,6]
        sd_data = alpha_array[:,7]
        n_vessels_data = alpha_array[:,8]
        n_perfused_vessels_data = alpha_array[:,9]
        n_unperfused_vessels_data = alpha_array[:,10]
        mean_diameter_data = alpha_array[:,11]
        mean_geometric_resistance_data = alpha_array[:,12]  
        diameter_binned_data = alpha_array[:,13]  
        flow_rate_binned_data = alpha_array[:,14]  
        n_connected_components_data = alpha_array[:,15]  
        n_cycles_data = alpha_array[:,16]  
        size_largest_connected_component_data = alpha_array[:,17]  
        n_cycles_largest_connected_component_data = alpha_array[:,18]  
        hypoxic_fraction_composite = np.vstack([hypoxic_fraction_composite, hypoxic_fraction_data]) if hypoxic_fraction_composite.size else hypoxic_fraction_data
        mean_composite = np.vstack([mean_composite, mean_data]) if mean_composite.size else mean_data
        min_composite = np.vstack([min_composite, min_data]) if min_composite.size else min_data
        half_composite = np.vstack([half_composite, half_data]) if half_composite.size else half_data
        max_composite = np.vstack([max_composite, max_data]) if max_composite.size else max_data
        sd_composite = np.vstack([sd_composite, sd_data]) if sd_composite.size else sd_data
        n_vessels_composite = np.vstack([n_vessels_composite, n_vessels_data]) if n_vessels_composite.size else n_vessels_data
        n_perfused_vessels_composite = np.vstack([n_perfused_vessels_composite, n_perfused_vessels_data]) if n_perfused_vessels_composite.size else n_perfused_vessels_data
        n_unperfused_vessels_composite = np.vstack([n_unperfused_vessels_composite, n_unperfused_vessels_data]) if n_unperfused_vessels_composite.size else n_unperfused_vessels_data
        mean_diameter_composite = np.vstack([mean_diameter_composite, mean_diameter_data]) if mean_diameter_composite.size else mean_diameter_data
        mean_geometric_resistance_composite = np.vstack([mean_geometric_resistance_composite, mean_geometric_resistance_data]) if mean_geometric_resistance_composite.size else mean_geometric_resistance_data
        diameter_binned_composite = np.vstack([diameter_binned_composite, diameter_binned_data]) if diameter_binned_composite.size else diameter_binned_data
        flow_rate_binned_composite = np.vstack([flow_rate_binned_composite, flow_rate_binned_data]) if flow_rate_binned_composite.size else flow_rate_binned_data
        n_connected_components_composite = np.vstack([n_connected_components_composite, n_connected_components_data]) if n_connected_components_composite.size else n_connected_components_data
        n_cycles_composite = np.vstack([n_cycles_composite, n_cycles_data]) if n_cycles_composite.size else n_cycles_data
        size_largest_connected_component_composite = np.vstack([size_largest_connected_component_composite, size_largest_connected_component_data]) if size_largest_connected_component_composite.size else size_largest_connected_component_data
        n_cycles_largest_connected_component_composite = np.vstack([n_cycles_largest_connected_component_composite, n_cycles_largest_connected_component_data]) if n_cycles_largest_connected_component_composite.size else n_cycles_largest_connected_component_data
    return alpha_array, hypoxic_fraction_composite, mean_composite, min_composite, half_composite, max_composite, sd_composite, n_vessels_composite, n_perfused_vessels_composite, n_unperfused_vessels_composite, mean_diameter_composite, mean_geometric_resistance_composite, diameter_binned_composite, flow_rate_binned_composite, n_connected_components_composite, n_cycles_composite, size_largest_connected_component_composite, n_cycles_largest_connected_component_composite 

# Define a function to filter the simulation stats based on alpha values and return individual matrices
def filter_by_alpha_fork_stoch(alpha_list, solver_stats):
    hypoxic_fraction_composite = np.array([])
    anoxic_fraction_composite = np.array([])
    mean_composite = np.array([])
    min_composite = np.array([])
    half_composite = np.array([])
    max_composite = np.array([])
    sd_composite = np.array([])  
    n_vessels_composite = np.array([])
    n_perfused_vessels_composite = np.array([])
    n_unperfused_vessels_composite = np.array([])
    mean_diameter_composite = np.array([])
    mean_geometric_resistance_composite = np.array([])
    diameter_binned_composite = np.array([])  
    flow_rate_binned_composite = np.array([])
    n_connected_components_composite = np.array([]) 
    n_cycles_composite = np.array([])
    size_largest_connected_component_composite = np.array([])
    n_cycles_largest_connected_component_composite = np.array([])
    for alpha_value in alpha_list:
        alpha_array = solver_stats[(solver_stats[:,0]==float(alpha_value))]
        hypoxic_fraction_data = alpha_array[:,2]  # extract data between identifiers and basic stats (i.e., hypoxic fractions)
        anoxic_fraction_data = alpha_array[:,3]  # extract data between identifiers and basic stats (i.e., hypoxic fractions)
        mean_data = alpha_array[:,4]
        min_data = alpha_array[:,5]
        half_data = alpha_array[:,6]
        max_data = alpha_array[:,7]
        sd_data = alpha_array[:,8]
        n_vessels_data = alpha_array[:,9]
        n_perfused_vessels_data = alpha_array[:,10]
        n_unperfused_vessels_data = alpha_array[:,11]
        mean_diameter_data = alpha_array[:,12]
        mean_geometric_resistance_data = alpha_array[:,13]  
        diameter_binned_data = alpha_array[:,14]  
        flow_rate_binned_data = alpha_array[:,15]  
        n_connected_components_data = alpha_array[:,16]  
        n_cycles_data = alpha_array[:,17]  
        size_largest_connected_component_data = alpha_array[:,18]  
        n_cycles_largest_connected_component_data = alpha_array[:,19]  
        hypoxic_fraction_composite = np.vstack([hypoxic_fraction_composite, hypoxic_fraction_data]) if hypoxic_fraction_composite.size else hypoxic_fraction_data
        anoxic_fraction_composite = np.vstack([anoxic_fraction_composite, anoxic_fraction_data]) if anoxic_fraction_composite.size else anoxic_fraction_data
        mean_composite = np.vstack([mean_composite, mean_data]) if mean_composite.size else mean_data
        min_composite = np.vstack([min_composite, min_data]) if min_composite.size else min_data
        half_composite = np.vstack([half_composite, half_data]) if half_composite.size else half_data
        max_composite = np.vstack([max_composite, max_data]) if max_composite.size else max_data
        sd_composite = np.vstack([sd_composite, sd_data]) if sd_composite.size else sd_data
        n_vessels_composite = np.vstack([n_vessels_composite, n_vessels_data]) if n_vessels_composite.size else n_vessels_data
        n_perfused_vessels_composite = np.vstack([n_perfused_vessels_composite, n_perfused_vessels_data]) if n_perfused_vessels_composite.size else n_perfused_vessels_data
        n_unperfused_vessels_composite = np.vstack([n_unperfused_vessels_composite, n_unperfused_vessels_data]) if n_unperfused_vessels_composite.size else n_unperfused_vessels_data
        mean_diameter_composite = np.vstack([mean_diameter_composite, mean_diameter_data]) if mean_diameter_composite.size else mean_diameter_data
        mean_geometric_resistance_composite = np.vstack([mean_geometric_resistance_composite, mean_geometric_resistance_data]) if mean_geometric_resistance_composite.size else mean_geometric_resistance_data
        diameter_binned_composite = np.vstack([diameter_binned_composite, diameter_binned_data]) if diameter_binned_composite.size else diameter_binned_data
        flow_rate_binned_composite = np.vstack([flow_rate_binned_composite, flow_rate_binned_data]) if flow_rate_binned_composite.size else flow_rate_binned_data
        n_connected_components_composite = np.vstack([n_connected_components_composite, n_connected_components_data]) if n_connected_components_composite.size else n_connected_components_data
        n_cycles_composite = np.vstack([n_cycles_composite, n_cycles_data]) if n_cycles_composite.size else n_cycles_data
        size_largest_connected_component_composite = np.vstack([size_largest_connected_component_composite, size_largest_connected_component_data]) if size_largest_connected_component_composite.size else size_largest_connected_component_data
        n_cycles_largest_connected_component_composite = np.vstack([n_cycles_largest_connected_component_composite, n_cycles_largest_connected_component_data]) if n_cycles_largest_connected_component_composite.size else n_cycles_largest_connected_component_data
    return alpha_array, hypoxic_fraction_composite, anoxic_fraction_composite, mean_composite, min_composite, half_composite, max_composite, sd_composite, n_vessels_composite, n_perfused_vessels_composite, n_unperfused_vessels_composite, mean_diameter_composite, mean_geometric_resistance_composite, diameter_binned_composite, flow_rate_binned_composite, n_connected_components_composite, n_cycles_composite, size_largest_connected_component_composite, n_cycles_largest_connected_component_composite 

# Define a function to filter the simulation stats based on alpha values and return individual matrices
def filter_by_mean_hex(alpha_list, solver_stats):
    hypoxic_fraction_composite = np.array([])
    mean_composite = np.array([])
    min_composite = np.array([])
    half_composite = np.array([])
    max_composite = np.array([])
    sd_composite = np.array([])  
    n_vessels_composite = np.array([])
    n_perfused_vessels_composite = np.array([])
    n_unperfused_vessels_composite = np.array([])
    mean_diameter_composite = np.array([])
    mean_geometric_resistance_composite = np.array([])
    diameter_binned_composite = np.array([])  
    flow_rate_binned_composite = np.array([])
    n_connected_components_composite = np.array([]) 
    n_cycles_composite = np.array([])
    size_largest_connected_component_composite = np.array([])
    n_cycles_largest_connected_component_composite = np.array([])
    for alpha_value in alpha_list:
        alpha_array = solver_stats[(solver_stats[:,0]==float(alpha_value))]
        hypoxic_fraction_data = alpha_array[:,2:3]  # extract data between identifiers and basic stats (i.e., hypoxic fractions)
        mean_data = alpha_array[:,3]
        min_data = alpha_array[:,4]
        half_data = alpha_array[:,5]
        max_data = alpha_array[:,6]
        sd_data = alpha_array[:,7]
        n_vessels_data = alpha_array[:,8]
        n_perfused_vessels_data = alpha_array[:,9]
        n_unperfused_vessels_data = alpha_array[:,10]
        mean_diameter_data = alpha_array[:,11]
        mean_geometric_resistance_data = alpha_array[:,12]  
        diameter_binned_data = alpha_array[:,13:17]  
        flow_rate_binned_data = alpha_array[:,17:21]  
        n_connected_components_data = alpha_array[:,21]  
        n_cycles_data = alpha_array[:,22]  
        size_largest_connected_component_data = alpha_array[:,23]  
        n_cycles_largest_connected_component_data = alpha_array[:,24]  
        hypoxic_fraction_composite = np.vstack([hypoxic_fraction_composite, hypoxic_fraction_data]) if hypoxic_fraction_composite.size else hypoxic_fraction_data
        mean_composite = np.vstack([mean_composite, mean_data]) if mean_composite.size else mean_data
        min_composite = np.vstack([min_composite, min_data]) if min_composite.size else min_data
        half_composite = np.vstack([half_composite, half_data]) if half_composite.size else half_data
        max_composite = np.vstack([max_composite, max_data]) if max_composite.size else max_data
        sd_composite = np.vstack([sd_composite, sd_data]) if sd_composite.size else sd_data
        n_vessels_composite = np.vstack([n_vessels_composite, n_vessels_data]) if n_vessels_composite.size else n_vessels_data
        n_perfused_vessels_composite = np.vstack([n_perfused_vessels_composite, n_perfused_vessels_data]) if n_perfused_vessels_composite.size else n_perfused_vessels_data
        n_unperfused_vessels_composite = np.vstack([n_unperfused_vessels_composite, n_unperfused_vessels_data]) if n_unperfused_vessels_composite.size else n_unperfused_vessels_data
        mean_diameter_composite = np.vstack([mean_diameter_composite, mean_diameter_data]) if mean_diameter_composite.size else mean_diameter_data
        mean_geometric_resistance_composite = np.vstack([mean_geometric_resistance_composite, mean_geometric_resistance_data]) if mean_geometric_resistance_composite.size else mean_geometric_resistance_data
        diameter_binned_composite = np.vstack([diameter_binned_composite, diameter_binned_data]) if diameter_binned_composite.size else diameter_binned_data
        flow_rate_binned_composite = np.vstack([flow_rate_binned_composite, flow_rate_binned_data]) if flow_rate_binned_composite.size else flow_rate_binned_data
        n_connected_components_composite = np.vstack([n_connected_components_composite, n_connected_components_data]) if n_connected_components_composite.size else n_connected_components_data
        n_cycles_composite = np.vstack([n_cycles_composite, n_cycles_data]) if n_cycles_composite.size else n_cycles_data
        size_largest_connected_component_composite = np.vstack([size_largest_connected_component_composite, size_largest_connected_component_data]) if size_largest_connected_component_composite.size else size_largest_connected_component_data
        n_cycles_largest_connected_component_composite = np.vstack([n_cycles_largest_connected_component_composite, n_cycles_largest_connected_component_data]) if n_cycles_largest_connected_component_composite.size else n_cycles_largest_connected_component_data
    return alpha_array, hypoxic_fraction_composite, mean_composite, min_composite, half_composite, max_composite, sd_composite, n_vessels_composite, n_perfused_vessels_composite, n_unperfused_vessels_composite, mean_diameter_composite, mean_geometric_resistance_composite, diameter_binned_composite, flow_rate_binned_composite, n_connected_components_composite, n_cycles_composite, size_largest_connected_component_composite, n_cycles_largest_connected_component_composite 

# Define a function to plot the diameter distribution of a network and return the diameter data
def get_diameter_distribution(reference_network_path):
    
    # Get the .vtk data
    point_coordinates_data_array, point_data_array, cell_data_array, polydata = get_vtk_data(reference_network_path)
    
    # Convert the radii in metres into the right units
    segment_diameters_um = cell_data_array['Vessel Radius m']*2*1000000  # convert radii (m) to diameters (um)
    mean_diameter = np.mean(segment_diameters_um)
    std_dev = np.std(segment_diameters_um)
    
    # Set the figure layout
    fig, axs = plt.subplots(nrows=1, ncols=1, sharex=True)
    fig.subplots_adjust(hspace = 0.75, wspace=.25)

    # Plot the distribution
    n, bins, patches = axs.hist(x=segment_diameters_um, bins='auto', alpha=0.7)
    axs.set_ylabel('number of vessels') 
    axs.set_xlabel('vessel diameter (μm)')    
    axs.title.set_text('${σ}$ = ' + str(round(std_dev,2)) + ' μm')
    axs.tick_params(labelbottom=True)
    axs.axvline(mean_diameter, c='black', ls='--', label='${µ}$ = ' + str(round(mean_diameter,2)) + ' μm')
    axs.set_ylim(0, 75)
    axs.legend()
    
    # Return the diameter stats
    return segment_diameters_um, mean_diameter, std_dev

# Define a function to bin an array according to three thresholds
def filter_by_thresholds(array, threshold_list):

    # Extract the thresholds
    threshold_1 = threshold_list[0]
    threshold_2 = threshold_list[1]
    threshold_3 = threshold_list[2]
    
    # Bin the array
    subset_1 = sum(x < threshold_1 for x in array)
    subset_2 = sum(threshold_1 <= x < threshold_2 for x in array)
    subset_3 = sum(threshold_2 <= x < threshold_3 for x in array)
    subset_4 = sum(threshold_3 <= x for x in array)

    # Return the bins
    return subset_1, subset_2, subset_3, subset_4

# Define a function to return the Betti curves for a network (based on Bernadette Stolz's code)
def PH_betti2(Weighted_adjacency_matrix):
    
#    beta0=np.array([])
#    beta1=np.array([])
#    size_largest_connected_component = np.array([])
#    biggest1 = np.array([])
    
#    for threshold in thresholds:
        
#    print('Threshold = ',step)
    adj = np.array(Weighted_adjacency_matrix)
#    adj[adj <= step] = 0
    
    # we don't want any nodes showing up as separate connected components
    adj = adj[~np.all(adj == 0, axis=1), :] #rows
    adj = adj[:,~np.all(adj == 0, axis=0)] #cols - logically axis should be =1 here, but code only works if I set it to 0

    n_nodes = adj.shape[0]
    
    G = nx.Graph(adj)
    number_of_connected_components = nx.number_connected_components(G)
    
#    beta0=np.append(beta0, number_of_connected_components)
    
    H = max(nx.connected_components(G), key=len)
    H_subgraph = G.subgraph(H).copy()
    size_largest_connected_component = H_subgraph.number_of_nodes()
#        size_largest_connected_component = np.append(size_largest_connected_component, size_largest_cc)

    #computes Betti-1
    n_edges = G.number_of_edges()
    n_cycles = number_of_connected_components - n_nodes + n_edges
#        beta1=np.append(beta1, n_cycles)
    
    
    #computes Beti-1 of largest connected component
    n_cycles_largest_connected_component = 1 - size_largest_connected_component + H_subgraph.number_of_edges()
#        biggest1 = np.append(biggest1,n_cycles_largest_connected_component)
    
    return number_of_connected_components, n_cycles, size_largest_connected_component, n_cycles_largest_connected_component
