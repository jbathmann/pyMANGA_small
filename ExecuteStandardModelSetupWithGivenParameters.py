#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 20 14:35:40 2018

@author: bathmann
"""
import sys
sys.path.append('./pybettina/')
import Land
import Flora
import SALT

class Run:
    def __init__(self, working_directory, setup_name, land_name, flora_name,
                 output_midstring, postfix_vtu_files, land_length_x, 
                 land_length_y, land_length_z, land_origin_x, land_origin_y, 
                 land_layers_x, land_layers_y, land_layers_z, 
                 pressure_variable_name, concentration_variable_name, 
                 pressure_initial_name, concentration_initial_name, 
                 darcy_velocity_initial_name, node_id_name, ini_darcy_function,
                 ini_pressure_function, ini_concentration_function,
                 tree_species, initial_plants, flora_plant_function, 
                 bettina_delta_t, number_of_bettina_timesteps, 
                 ogs_time_delta_ts, ogs_timerepeats, ogs_outputdeltaN,
                 ogs_outputrepeats):
        land = Land.Land(setup_name + land_name, working_directory)
        land.create3DRectangularLand(setup_name + land_name, land_origin_x, 
                                     land_origin_y, -land_length_z,
                                     land_length_x, land_length_y, 
                                     land_length_z, land_layers_x + 1, 
                                     land_layers_y + 1, land_layers_z + 1)
        
        land.setCIniPIniAndNodeIds(c_name = concentration_initial_name, 
                          p_name = pressure_initial_name, 
                          q_name = darcy_velocity_initial_name,
                          node_id_name = node_id_name, 
                          darcy_velocity_function = ini_darcy_function, 
                          pressure_function = ini_pressure_function, 
                          concentration_function= ini_concentration_function)    
        
        land.setSurfacePointLocations()
        land.outputLand()
        
        flora = Flora.Flora(setup_name + flora_name, land, working_directory)
        flora_plant_function(flora, land)
        
        
        model = SALT.SaltSetup(setup_name, working_directory, land, flora, 
                               True, output_midstring)
        model.setVariableNames(pressure_variable_name, 
                               concentration_variable_name)
        model.setInitialConditionName(pressure_initial_name, 
                                      concentration_initial_name, 
                                      darcy_velocity_initial_name)
        
        model.createBoundarySurface("left")
        model.createBoundarySurface("right")
        model.updateBoundaryConditions()
        model.createMeshCollection(output_midstring, postfix_vtu_files)
        model.createTreeCollection( tree_species, postfix_vtu_files)
        model.createFloraCollection(flora_name, postfix_vtu_files)
        
        for i in range(number_of_bettina_timesteps):
        
            print("Simulating Month ", (i+1))
            model.updateBoundaryConditions()
            t_ini = i * bettina_delta_t
            t_end = (i + 1) * bettina_delta_t
            model.setAndRunOgs(t_ini, t_end, ogs_timerepeats, 
                       ogs_time_delta_ts, ogs_outputrepeats, ogs_outputdeltaN)
            model.readAndPassTMeshFileNames(t_ini, t_end)
            
            model.updateModel()