#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example file. This script generates a population of trees on a 3 dimensional
land domain and evolves it in time.
@date: 2018-Today
@author: jasper.bathmann@ufz.de
"""

import ExecuteStandardModelSetupWithGivenParameters
# ExecuteStandardModelSetupWithGivenParameters.Run(args) starts the salt model

# -- Parameters on local working directories and setup naming ###

# working_directory: Directory, where simulation results are saved. Please make
# sure the location exists on your local machine
working_directory = "./testruns/"
# setup_name: This string is contained in all output files generated
setup_name = "testruns"
# land_name: specific name, which is contained on all output files for the
# meshes representing the land domain
land_name = "_exampleland"
# flora_name: specific name, which is contained on all output files for the
# meshes associated with the flora
flora_name = "_exampleflora"
# output_midstring: specific name, which is contained on all output files for
# the meshes associated with the ogs output
output_midstring = "_ogsOutput"
# postfix_vtu_files: specific string, which is at the end of all vtu files
# saved on within the working directory
postfix_vtu_files = ".vtu"


# -- Parameters on land domain properties ###

# - Geometry ##
# land_length_x,y,z: the dimensions (in m) of the simulated land domain
land_length_x, land_length_y, land_length_z = 5, 5, 3
# land_origin_x,y: the coordinates of the upper-left-bottom corner of the land
# mesh. The z_coordinate of the upper_left_bottom corner is always set to
# z = -land_length_z
land_origin_x, land_origin_y = 0, 0
# land_layers_x,y,z: the number of layers in each dimension.
land_layers_x, land_layers_y, land_layers_z = 25, 25, 1

# - Names for primary variables in subsurface processes ##
# pressure_variable_name, concentration_variable_name: Names of the primary
# variables themselves
pressure_variable_name, concentration_variable_name = "pressure",\
    "concentration"
# pressure_initial_name, concentration_initial_name: names for the initial con-
# ditions for the primary variables
pressure_initial_name, concentration_initial_name = "p_ini", "c_ini"
# darcy_velocity_initial_name: name of the initial darcy velocity distribution.
# Necessary to define 2nd type boundary conditions
darcy_velocity_initial_name = "q_ini"
# nade_id_name: name of the property vector containing the node ids. Necessary
# for ogs to run properly.
node_id_name = "bulk_node_ids"
# ini_darcy/pressure/concentration_function(point): functions to create intial
# conditions. Note: all functions must depend on point, which is a (3,1) tuple
# in order to allow for spatially depending initial distributions
dp_dx = 5e-3


def ini_darcy_function(point, dp_dx=dp_dx):
    return -1000*9.81*dp_dx*1.239*1e-5


def ini_pressure_function(point, dp_dx=dp_dx):
    return -1000*9.81*(point[2]+dp_dx*point[0])


def ini_concentration_function(point):
    return 0.035


# -- Parameters on flora properties ###

# tree species: list containing the species which considered in the initial
# tree distribution. Given in the form: [species1, species2, ...], with
# species_i being a string with the species name
tree_species = ["Avicennia"]
# initial_plants: number of individums planted for each species for the initial
# plant distributions. Given in the form [n_species1, n_species2, ...], with
# n_speciesi being the number of individuums of each species
initial_plants = [150]
# flora_plant_function: function defining how the initial plant distribution
# has to look like


def flora_plant_function(flora, land):
    flora.randomlyPlantTreesInRectangularDomain(
            initial_plants, tree_species, land.bounding_box)


# -- Parameters on time loop ###
# bettina_timesteps: length of one timestep in bettina in [s]. Note: half a
# year corresponds to 15778800.0 seconds
bettina_delta_t = 15778800.0/(6.)
# number_of_bettina_timesteps: total number of iterations of the bettina model
number_of_bettina_timesteps = 50 * 12
# ogs_timerepeats
# ogs_time_delta_ts: list containing different timestep lengths in s for the
# ogs simulation
ogs_time_delta_ts = [1e-1, 1e0, 60, 300, 900, 1800,
                     3600, 3600*2, 3600*4, 3600*8]
# ogs_timerepeats: list of the same length as ogs_time_delta_ts containing the
# number of iterations with the different step sizes given in ogs_time_delta_ts
ogs_timerepeats = [1,   1,  10, 10,  10,  50,   100,   100,     100,  1000]
# ogs_outputdeltaN: list containing number of steps for after which ogs is
# writing an output file
ogs_outputdeltaN = [10000]
# ogs_outputrepeats: list of same shape as ogs_outputdeltaN containing the
# number of iterations with the different output intervals given in
# ogs_outputdeltaN
ogs_outputrepeats = [1]


# -- Salt model execution file import and run ###
# with the configuration given above
ExecuteStandardModelSetupWithGivenParameters.Run(
        working_directory, setup_name,
        land_name, flora_name, output_midstring, postfix_vtu_files,
        land_length_x, land_length_y, land_length_z, land_origin_x,
        land_origin_y, land_layers_x, land_layers_y, land_layers_z,
        pressure_variable_name, concentration_variable_name,
        pressure_initial_name, concentration_initial_name,
        darcy_velocity_initial_name, node_id_name, ini_darcy_function,
        ini_pressure_function, ini_concentration_function,
        tree_species, initial_plants, flora_plant_function,
        bettina_delta_t, number_of_bettina_timesteps,
        ogs_time_delta_ts, ogs_timerepeats, ogs_outputdeltaN,
        ogs_outputrepeats)
