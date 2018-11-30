#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 12:21:17 2018

@author: bathmann
"""
import sys
sys.path.append('./pybettina/')
import Land
import Flora
import SALT
import numpy as np
working_directory = "./testruns/"
setup_name ="testruns"
prefix = setup_name + "_Output_pcs"
postfix = ".vtu"
land = Land.Land("testtestrunsmesh", working_directory)
land.create3DRectangularLand("testtestrunsmesh", 0, 0, -1, 10, 10, 1, 51, 51, 2)
n = land.initial_mesh.grid.GetNumberOfPoints()
c,p ,q, iD = [], [], [], []
dp_dx = 1e-3
for i in range(n):
    point = land.initial_mesh.grid.GetPoints().GetPoint(i)
    
    q.append(-(1000+0*35*.7)*9.81*dp_dx*1.239*1e-5)
    p.append(-(1000+0*35*.7)*9.81*(point[2]+dp_dx*point[0]))
    c.append(.035)
    iD.append(i)
land.setCIniPIniAndNodeIds(["c_ini", "p_ini"], "bulk_node_ids", [np.array(c), np.array(p)], np.array(iD))
land.outputLand()

flora = Flora.Flora("testrunsflora", "testrunsconstants", land, working_directory)
flora.randomlyPlantTreesInRectangularDomain([10],["Avicennia"],land.bounding_box)


model = SALT.SaltSetup(setup_name, working_directory, land, flora)
model.setVariableNames("pressure","concentration")
model.setInitialConditionName("p_ini", "c_ini", "q_ini")
model.setQiniArray(np.array(q))
model.createBoundarySurface("left")
model.createBoundarySurface("right")
model.updateBoundaryConditions()
model.createMeshCollection(prefix, postfix)
model.createTreeCollection( "Avicennia", ".vtu")
model.createFloraCollection( "prefix" + "flora_grid", ".vtu")
#1/2 Jahr = 15778800.0 Sekunde34n
dt = 15778800.0/(6.)
for i in range(50*12):

    print("Simulating Month ", (i+1))
    model.updateBoundaryConditions()
    t_ini = i * dt
    t_end = (i + 1) * dt
    #1 Tag = [100, 50, 59, 23 mit [1e-1, 1e0, 60, 3600, ]
    timerepeats = [1,   1,  10, 10,  10,  50,   100,   100,     100,  1000]#[100, 50, 59, 23, 29*24, 30*5*24]#
    timedeltaTs = [1e-1, 1e0, 60, 300, 900, 1800, 3600, 3600*2, 3600*4,3600*8]#[1e-1, 1e0, 60, 3600, 3600, 3600]#    outputrepeats =[1]# [1, 29, 30*5]
    outputrepeats =[1]# [1, 29, 30*5]
    outputdeltaN = [10000]#[232, 24, 24]
    model.setAndRunOgs(t_ini, t_end, timerepeats, timedeltaTs, outputrepeats,
                       outputdeltaN)
    model.readAndPassTMeshFileNames(t_ini, t_end)
    
    model.updateModel()
    