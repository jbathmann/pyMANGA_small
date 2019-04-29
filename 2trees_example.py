#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example script. This file passes parameters to the salt model, such that 2
trees with one of them downstream to the other are located on a 3 dimensional
land domain.
@date: 2018-Today
@author: jasper.bathmann@ufz.de
"""
import ExecuteStandardModelSetupWithGivenParameters
# ExecuteStandardModelSetupWithGivenParameters.Run(args) starts the salt model
# with the configuration given below

# -- Parameters on local working directories and setup naming ###

# working_directory: Directory, where simulation results are saved. Please make
# sure the location exists on your local machine -- type:string
working_directory = "./testcases/2trees/"
# setup_name: This string is contained in all output files generated
#  -- type:string
setup_name = "2trees"
# land_name: specific name, which is contained on all output files for the
# meshes representing the land domain -- type:string
land_name = "_exampleland"
# flora_name: specific name, which is contained on all output files for the
# meshes associated with the flora -- type:string
flora_name = "_exampleflora"
# output_midstring: specific name, which is contained on all output files for
# the meshes associated with the ogs output -- type:string
output_midstring = "_ogsOutput"
# postfix_vtu_files: specific string, which is at the end of all vtu files
# saved on within the working directory -- type:string
postfix_vtu_files = ".vtu"


# -- Parameters on land domain properties ###

# - Geometry ##
# land_length_x,y,z: the dimensions (in m) of the simulated land domain
# -- type:float
land_length_x, land_length_y, land_length_z = 12, 6, 5
# land_origin_x,y: the coordinates of the upper-left-bottom corner of the land
# mesh. The z_coordinate of the upper_left_bottom corner is always set to
# z = -land_length_z -- type:float
land_origin_x, land_origin_y = 0, 0
# land_layers_x,y,z: the number of layers in each dimension. -- type:int
land_layers_x, land_layers_y, land_layers_z = 48, 24, 10

# - Names for primary variables in subsurface processes ##
# pressure_variable_name, concentration_variable_name: Names of the primary
# variables themselves -- type:string
pressure_variable_name, concentration_variable_name = \
                                                "pressure", "concentration"
# pressure_initial_name, concentration_initial_name: names for the initial con-
# ditions for the primary variables -- type:string
pressure_initial_name, concentration_initial_name = "p_ini", "c_ini"
# darcy_velocity_initial_name: name of the initial darcy velocity distribution.
# Necessary to define 2nd type boundary conditions -- type:string
darcy_velocity_initial_name = "q_ini"
# node_id_name: name of the property vector containing the node ids. Necessary
# for ogs to run properly. -- type:string
node_id_name = "bulk_node_ids"
# ini_darcy/pressure/concentration_function(point): functions to create intial
# conditions. Note: all functions must depend on point, which is a (3,1) tuple
# in order to allow for spatially depending initial distributions
# -- type:python function declarations
dp_dx = 1e-3


def ini_darcy_function(point, dp_dx=dp_dx):
    return -1000*(1+.701*.035)*9.81*dp_dx*1.239*1e-5


def ini_pressure_function(point, dp_dx=dp_dx):
    return -1000*(1+.701*.035)*9.81*(point[2]+dp_dx*point[0])


def ini_concentration_function(point):
    return 0.035


# ## Parameters on flora properties ###

# tree species: list containing the species which considered in the initial
# tree distribution. Given in the form: [species1, species2, ...], with
# species_i being a string with the species name -- type:list of strings
tree_species = ["Avicennia"]
# initial_plants: number of individums planted for each species for the initial
# plant distributions. Given in the form [n_species1, n_species2, ...], with
# n_speciesi being the number of individuums of each species
# -- type:list of ints
initial_plants = [150]
# flora_plant_function: function defining how the initial plant distribution
# has to look like -- type:python function declarations


def flora_plant_function(flora, land):
    from pybettina import Tree
    new_tree = Tree.Tree(3., 3., flora.land, "Avicennia", 0,
                         flora.flora_name)
    new_tree.plantTree(flora.working_directory)
    flora.trees.append(new_tree)

    second_tree = Tree.Tree(6, 3, flora.land, "Avicennia", 1,
                            flora.flora_name)
    second_tree.plantTree(flora.working_directory)
    flora.trees.append(second_tree)


# ## Parameters on time loop ###
# bettina_timesteps: length of one timestep in bettina in [s]. Note: half a
# year corresponds to 15778800.0 seconds -- type:double
bettina_delta_t = 15778800.0/(6.)
# number_of_bettina_timesteps: total number of iterations of the bettina model
# -- type:int
number_of_bettina_timesteps = 50 * 12
# ogs_timerepeats
# fixed_output_times: list containing different fixed output times in s for the
# ogs simulation -- type:list of doubles
ogs_fixed_output_times = [60 * 60 * 24 * 5, 60 * 60 * 24 * 10,
                          60 * 60 * 24 * 15, 60 * 60 * 24 * 20,
                          60 * 60 * 24 * 25, 60 * 60 * 24 * 30]
# ogs_outputdeltaN: list containing number of steps for after which ogs is
# writing an output file -- type:list of ints
ogs_outputdeltaN = [250]
# ogs_outputrepeats: list of same shape as ogs_outputdeltaN containing the
# number of iterations with the different output intervals given in
# ogs_outputdeltaN -- type:list of ints
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
        ogs_fixed_output_times, ogs_outputdeltaN,
        ogs_outputrepeats)
