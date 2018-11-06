#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 15:35:00 2018

@author: bathmann

Test file creator for simple one tree scenario:
This file script creates a testcase for the time dependant boundary conditions in OGS. 
The subsurface mesh is created in 3D with dimensions 10m,2m,1m (x,y,z) with a grid spacing
of 2cm. The subsurface mesh is stored as "SmallBox1Tree.vtu".
Parameters for the boundary mesh are given according to the current definition of avicinnia
parameters. The boundary mesh is stored as "0Avicennia000000000.vtu".
Additionally it will run simulations with feedback between the current bettina model and
ogs, if the bool useOgs is set True. Hereby subsurface data is stored after each 1e5 seconds,
tree data is updated every 1e6 seconds and the simulationtime is set to 1e6 seconds.
"""
import sys
sys.path.insert(0, '../')
sys.path.insert(0, '../mesh_tools')
sys.path.insert(0, '../ogs_interface')
import Land
import MeshEditor
import numpy as np


projectdir = "./TestCase_SteadyState/"#directory, where output is stored
landName = "HC_Testcase-SteadyState"#name of subsurface mesh
prefix = "Output_"+landName #prefix for ogs output storage
iterations = 1 #number of iterations
delta_t = 1e3#data storage frequency
t_ini = 0 #initial time
t_end = 1e7#first flora update time
l=10 #defines reference length in meters
parts=51 #defines grid spacing
geometricalRatio=[1,.1,.1]#defines ratio between (xlen,ylen,zlen)
horizontal_gradient = -1e-3 # defines horizontal landscape slope
vertical_gradient = 0 # useless
g= 9.81 # gravitational acceleration
n_avi = 1 #number of avicennias planted
n_rhi = 0 #number of rhi planted
random = False # bool for type of concentration distribution
meanc = .035 #mean component concentration
useOgs = True #bool for interaction with OGS


def calculateAnalyticalSolution(darcyvel, Dm, beta_t, L, c_ini, mesh):
    D_h = Dm*0.001+beta_t*abs(darcyvel)
    print(darcyvel)
    Kstar= -darcyvel/D_h
    prefac = c_ini * L /(L-1/Kstar*(1-np.exp(-Kstar*L)))
    N = mesh.GetPoints().GetNumberOfPoints()
    values = []
    for i in range(N):
        x = mesh.GetPoints().GetPoint(i)[0]
        wx = prefac * (1-np.exp(-Kstar*x))
        values.append(wx)
    return values


#Setup of land mesh and project
land = Land.Land(landName,l, parts, geometricalRatio, horizontal_gradient, g, vertical_gradient, useOgs)
land.setProjectdir(projectdir)
land.setMeanC(meanc, random)
land.landMesh.linear=True
land.initializeMesh()

landMeshEdit = MeshEditor.MeshEditor(land.landMesh.pd)


Dm = 2e-3
beta_t = 0.5
land.project.parameter_values[1]=str(Dm)
land.project.parameter_values[4]=str(beta_t)
land.project.parameter_values[5]=".5"
land.project.parameter_values[-1]="1.e-11 0 0 0 1.e-11 0 0 0 1.e-11"
land.project.convergence_criterion_reltols = "5e-8 5e-8"
land.project.solver_linear_solver_lis = "-i bicgstab -p ilut -tol 1e-12 -maxiter 200000"
land.project.solver_linear_solver_petsc_parameters = "-hc_ksp_type bcgs -hc_pc_type bjacobi -hc_ksp_rtol 1e-10 -hc_ksp_max_it 20000"
land.project.solver_linear_solver_eigen_error_tol = "1e-11"
land.initializeHorizontalFlow()
darcyvel=land.getHotizontalDarcyVelBoundary()
omega = calculateAnalyticalSolution(darcyvel, Dm, beta_t, l, meanc, land.landMesh.pd)
landMeshEdit.AddPropertyVector(omega,"Analytical Steady State Solution","double")

#Initialization of flora
#First ogs time stepping
land.configureTimeSteppingProject(str(delta_t), str(t_end), str(t_ini), prefix)
land.project.output_each_steps = "100"

#Initialization of 1 Tree
#newTree= Tree.Tree(coordx, coordy, land, "Avicennia", 0)
#newTree.plantTree()
#flora.trees.append(newTree)
#flora.createTreeBCs()
land.landMesh.OutputMesh(projectdir)
#        self.parameter_values = [rho_v, Dm_v, retardation_v, decay_v, beta_l_v, beta_t_v, c_ini_v, 
#                                 p_ini_v, constant_porosity_parameter_v, kappa1_v]

land.project.writeProjectFile()    
land.getSalinitiesAndPressures()    

    
