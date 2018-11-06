#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 15:35:00 2018

@author: bathmann

"""
import sys
sys.path.insert(0, './pyogsproject')
sys.path.insert(0, './pybettina')

sys.path.insert(0, './pymeshinteraction')
import MeshInteractor
import numpy as np



        
timerepeats = [50,  595,  1140,3600]#,440,8500,1]
timedeltaTs = [1e-2,1e-1,1e-0,1e0]#,1e-1,1e-1,1]

outputrepeats = [1, 1,6]#,500,1]
outputdeltaN = [1185,600,600]#,5000,1]

projectdir = "./TestCase_VDBP_OutFlowConcentrationRegulated/"#directory, where output is stored
landName = "goswami_input"#name of subsurface mesh
prefix = "Goswami_Component_Transport" #prefix for ogs output storage
iterations = 1 #number of iterations
delta_t = 1e2#data storage frequency
t_ini = 0 #initial time
t_end = 4800#first flora update time
l=.53 #defines reference length in meters
parts=27#53 #defines grid spacing

factor = 2
geometricalRatio=[1,.05,0.49056603773584906]#defines ratio between (xlen,ylen,zlen)
horizontal_gradient = -1e-1 # defines horizontal landscape slope
vertical_gradient = 0 # useless
g= 9.81 # gravitational acceleration
n_avi = 1 #number of avicennias planted
n_rhi = 0 #number of rhi planted
random = False # bool for type of concentration distribution
linear = False
meanc = .035 #mean component concentration
target_c = .05
useOgs = True #bool for interaction with OGS
conc_difference_ratio = 0.701
class Land:
    def __init__(self,land_name):
        self.name = land_name
    def create3DRectangularLand(self, name, ox, oy, oz, lx, ly, lz, partsx, partsy, partsz):
        self.landMesh = MeshInteractor.MeshInteractor(name)
        self.landMesh.create3DRectangularGrid(partsx, ox, lx, partsy, oy, ly, partsz, oz, lz)
        
        
    def addCIniPIniAndNodeIds(self, c_name, p_name, node_id_name, c_ini, p_ini, node_id):
        self.landMesh.AddPropertyVector(p_ini,p_name,"double")
        self.landMesh.AddPropertyVector(c_ini,c_name,"double")
        self.landMesh.AddPropertyVector(node_id,node_id_name,"unsigned_long")

    def outputLand(self, location):
        self.landMesh.OutputMesh(location)
land = Land("testland")
land.create3DRectangularLand("test", 0, 0, 0, 10, 10, 10, 11, 11, 11)
c,p,iD = [],[],[]
n = land.landMesh.pd.GetNumberOfPoints()
for i in range(n):
    point = land.landMesh.pd.GetPoints().GetPoint(i)
    c.append(point[0])
    p.append(point[1])
    iD.append(33)
land.addCIniPIniAndNodeIds("c_ini", "p_ini", "bulk_node_ids", np.array(c), np.array(p), np.array(iD))
land.outputLand("./")