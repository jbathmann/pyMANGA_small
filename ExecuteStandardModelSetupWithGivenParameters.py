# !/usr/bin/env python3
#  -*- coding: utf-8 -*-
"""
Class to run model with given parameters:
# #  Parameters on local working directories and setup naming # # #
working_directory: Directory, where simulation results are saved. Please make
sure the location exists on your local machine
setup_name: This string is contained in all output files generated
land_name: specific name, which is contained on all output files for the
meshes representing the land domain
flora_name: specific name, which is contained on all output files for the
meshes associated with the flora
output_midstring: specific name, which is contained on all output files for the
meshes associated with the ogs output
postfix_vtu_files: specific string, which is at the end of all vtu files
saved on within the working directory

# #  Parameters on land domain properties # # #
#  Geometry # #
land_length_x,y,z: the dimensions (in m) of the simulated land domain
land_origin_x,y: the coordinates of the upper-left-bottom corner of the land
mesh. The z_coordinate of the upper_left_bottom corner is always set to
z = -land_length_z
land_layers_x,y,z: the number of layers in each dimension.
#  Names for primary variables in subsurface processes # #
pressure_variable_name, concentration_variable_name: Names of the primary
variables themselves
pressure_initial_name, concentration_initial_name: names for the initial con-
ditions for the primary variables
darcy_velocity_initial_name: name of the initial darcy velocity distribution.
Necessary to define 2nd type boundary conditions
nade_id_name: name of the property vector containing the node ids. Necessary
for ogs to run properly.
ini_darcy/pressure/concentration_function(point): functions to create intial
conditions. Note: all functions must depend on point, which is a (3,1) tuple
in order to allow for spatially depending initial distributions

# #  Parameters on flora properties # # #
tree species: list containing the species which considered in the initial
tree distribution. Given in the form: [species1, species2, ...], with speciesi
being a string with the species name
initial_plants: number of individums planted for each species for the initial
plant distributions. Given in the form [n_species1, n_species2, ...], with
n_speciesi being the number of individuums of each species
flora_plant_function: function defining how the initial plant distribution has
to look like

# #  Parameters on time loop # # #
bettina_timesteps: length of one timestep in bettina in [s]. Note: half a year
corresponds to 15778800.0 seconds
number_of_bettina_timesteps: total number of iterations of the bettina model
ogs_timerepeats
ogs_time_delta_ts: list containing different timestep lengths in s for the ogs
simulation
ogs_timerepeats: list of the same length as ogs_time_delta_ts containing the
number of iterations with the different step sizes given in ogs_time_delta_ts
ogs_outputdeltaN: list containing number of steps for after which ogs is
writing an output file
ogs_outputrepeats: list of same shape as ogs_outputdeltaN containing the
number of iterations with the different output intervals given in
ogs_outputdeltaN

@date: 2018-Today
@author: jasper.bathmann@ufz.de
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
        # land domain creation
        land = Land.Land(setup_name + land_name, working_directory)
        land.create3DRectangularLand(setup_name + land_name, land_origin_x,
                                     land_origin_y, -land_length_z,
                                     land_length_x, land_length_y,
                                     land_length_z, land_layers_x + 1,
                                     land_layers_y + 1, land_layers_z + 1)
        # initial conditions for land domain
        land.setCIniPIniAndNodeIds(
                c_name=concentration_initial_name,
                p_name=pressure_initial_name,
                q_name=darcy_velocity_initial_name,
                node_id_name=node_id_name,
                darcy_velocity_function=ini_darcy_function,
                pressure_function=ini_pressure_function,
                concentration_function=ini_concentration_function)

        land.setSurfacePointLocations()
        land.outputLand()
        # flora creation
        flora = Flora.Flora(setup_name + flora_name, land, working_directory)
        # initial plant distribution
        flora_plant_function(flora, land)

        # salt model setup
        model = SALT.SaltSetup(setup_name, working_directory, land, flora,
                               True, output_midstring)
        # variable name definition
        model.setVariableNames(pressure_variable_name,
                               concentration_variable_name)
        # initial condition array name definition
        model.setInitialConditionName(pressure_initial_name,
                                      concentration_initial_name,
                                      darcy_velocity_initial_name)
        # creation of boundary surfaces at given faces ("left", "right", "top",
        # "bottom", "front", "back") of land domain box
        model.createBoundarySurface("left")
        model.createBoundarySurface("right")
        # storing of boundary conditions in the model
        model.updateBoundaryConditions()

        # pvd_file-output definitions
        model.createMeshCollection(output_midstring, postfix_vtu_files)
        model.createTreeCollection(tree_species, postfix_vtu_files)
        model.createFloraCollection(flora_name + "_grid", postfix_vtu_files)

        # timeloop for model evolution
        for i in range(number_of_bettina_timesteps):

            print("Simulating Month ", (i+1))
            model.updateBoundaryConditions()
            t_ini = i * bettina_delta_t
            t_end = (i + 1) * bettina_delta_t
            # ogs execution
            model.setAndRunOgs(
                    t_ini, t_end, ogs_timerepeats,
                    ogs_time_delta_ts, ogs_outputrepeats, ogs_outputdeltaN)
            # analysing output and storing it for the bettina execution
            model.readAndPassTMeshFileNames(t_ini, t_end)
            # bettina timestep and model update
            model.updateModel()
