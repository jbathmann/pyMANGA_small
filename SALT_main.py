#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 15:35:00 2018

@author: bathmann

"""
import sys
sys.path.append('./pymeshinteraction/')
import MeshPointFinder as MPF
import MeshInteractor
import numpy as np


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 14:27:35 2018

@author: bathmann
"""
sys.path.append('./pybettina/')
import Land
import Flora
import Bettina
import SortFileNames as SFN
sys.path.append('./pyogsproject/')

import OGSProject
#import SortFilenames as SFn
#TODO: Check, where the script file is gone. Maybe it is saved on the ssd...
import numpy as np
##Dummy files for Jasper to test interaction with OGS:
class SaltSetup:
    def __init__(self, name, working_directory):
        self.setup_name = name
        self.working_directory = working_directory
        self.boundary_surfaces = []
        
    def setLand(self, land):
        self.land = land
    
    def setFlora(self, flora):
        self.flora = flora
    
    def initializeBettina(self, name):
        self.bettina = Bettina.Bettina(name, self.land, self.flora)
        
    def createOGSProject(self):
        self.ogsPrj = OGSProject.OGSProject(self.working_directory, self.setup_name + "_OGSproject")
        self.ogsPrj.setLandName(self.land.name)
        self.ogsPrj.initializeProject()
        
    def setInitialConditionNames(self, p_ini, c_ini):
        self.c_ini_name = c_ini
        self.p_ini_name = p_ini
        self.ogsPrj.setInitialConditionName(p_ini, c_ini)
        
    def setVariableNames(self, p_var, c_var):
        self.c_var_name = c_var
        self.p_var_name = p_var
        self.ogsPrj.setVariableNames(p_var, c_var)
        self.land.setCurrentPropertyNames([c_var,p_var])
        
    def setInitialConditionName(self, p_ini, c_ini):
        self.c_ini_name = c_ini
        self.p_ini_name = p_ini
        self.ogsPrj.setInitialConditionName(p_ini, c_ini)
        
    def addBoundaryConditionForSurfaces(self, variable, *args):
        self.ogsPrj.createBoundaryConditionsFromList(self.boundary_surfaces, variable, *args)
        
    def writeOgsProject(self):
        self.ogsPrj.writeProject()

    def createBoundarySurface(self, location):
        land_mesh = self.land.initial_mesh
        land_grid = land_mesh.grid
        origin = [land_mesh.ox, land_mesh.oy, land_mesh.oz]
        deltas = [land_mesh.dx, land_mesh.dy, land_mesh.dz]
        steps = [land_mesh.extension_points_X, land_mesh.extension_points_Y, land_mesh.extension_points_Z]
        shift = 0
        if location == "left" or location == "right":
            steps[0] = 1
            shift = 1
        if location == "right":
            origin[0] = origin[0] + land_mesh.lengthX
        if location == "front" or location == "back":
            steps[1] = 1
            shift = 2
        if location == "back":
            origin[1] = origin[1] + land_mesh.lengthY
        if location == "bottom" or location == "top": steps[2] = 1
        if location == "top":
            origin[2] = origin[2] + land_mesh.lengthZ
        boundary_mesh_name = land_mesh.meshName + "_" + location + "Boundary"
        boundary_creator = MPF.MeshPointFinder(land_grid)
        points, ids = boundary_creator.findPointsOnPlane(steps[0], steps[1], steps[2], deltas[0], deltas[1], deltas[2], origin[0], origin[1], origin[2]) 
        temppoints = []
        for i in range(len(points)):
            point = points[i]
            point = point[(0+shift)%3], point[(1+shift)%3], point[(2+shift)%3]
            temppoints.append(point)
        temp_boundary = MeshInteractor.MeshInteractor(boundary_mesh_name)
        temp_boundary.CreateMeshFromPoints([temppoints, ids])
        temp_boundary.CreateMultipleTriangles()
        cells = temp_boundary.grid.GetCells()
        boundary = MeshInteractor.MeshInteractor(boundary_mesh_name)
        boundary.CreateMeshFromPoints([points, ids])
        boundary.grid.SetCells(5,cells)
        full_c_ini = land_grid.GetPointData().GetArray(self.c_ini_name)
        full_p_ini = land_grid.GetPointData().GetArray(self.c_ini_name)
        boundary_c_ini, boundary_p_ini = np.zeros(len(ids)), np.zeros(len(ids)) 
        for i in range(len(ids)):
            iD = ids[i]
            boundary_c_ini[i]=full_c_ini .GetTuple(iD)[0]
            boundary_p_ini[i]=(full_p_ini .GetTuple(iD)[0])
        boundary.AddPropertyVector(np.array(boundary_c_ini),self.c_ini_name, "double")
        boundary.AddPropertyVector(np.array(boundary_p_ini),self.p_ini_name, "double")
        boundary.OutputMesh(self.working_directory)
        self.boundary_surfaces.append(boundary_mesh_name)
        
    def setTimeSteppingAndOutputLoopsForOgs(self,timerepeats, timedeltaTs, outputrepeats, outputdeltaN):
        self.ogsPrj.setTimeSteppingAndOutputLoops(timerepeats, timedeltaTs, outputrepeats, outputdeltaN)

    
working_directory = "./testruns/"
prefix = "Output_HC_Testcase-SteadyState_pcs"
postfix = ".vtu"
land = Land.Land("testlandmesh", working_directory)
land.create3DRectangularLand("testlandmesh", 0, 0, -1, 10, 10, 1, 11, 11, 3)
n = land.initial_mesh.grid.GetNumberOfPoints()
c,p , iD = [], [], []
for i in range(n):
    point = land.initial_mesh.grid.GetPoints().GetPoint(i)
    c.append(35)
    p.append(-1000*9.81*(point[2]-1e-2*point[0]))
    iD.append(i)
land.setCIniPIniAndNodeIds(["c_ini", "p_ini"], "bulk_node_ids", [np.array(c), np.array(p)], np.array(iD))
land.outputLand()
  
model = SaltSetup("testrun", working_directory)
model.setLand(land)
model.createOGSProject()
model.setVariableNames("pressure","concentration")
model.setInitialConditionName("p_ini", "c_ini")
model.createBoundarySurface("left")
model.createBoundarySurface("right")
model.addBoundaryConditionForSurfaces("pressure", "NonuniformDirichlet")
model.addBoundaryConditionForSurfaces("concentration", "NonuniformDirichlet")

timerepeats = [50,  595,  1140,3600]#,440,8500,1]
timedeltaTs = [1e-2,1e-1,1e-0,1e0]#,1e-1,1e-1,1]

outputrepeats = [1, 1,6]#,500,1]
outputdeltaN = [1185,600,600]#,5000,1]

model.setTimeSteppingAndOutputLoopsForOgs(timerepeats, timedeltaTs, outputrepeats, outputdeltaN)




model.writeOgsProject()

"""
n = land.initial_mesh.grid.GetNumberOfPoints()
c,p , iD = [], [], []
for i in range(n):
    point = land.initial_mesh.grid.GetPoints().GetPoint(i)
    c.append(point[0])
    p.append(point[1])
    iD.append(i)
land.setCIniPIniAndNodeIds(["c_ini", "p_ini"], "bulk_node_ids", [np.array(c), np.array(p)], np.array(iD))
land.setCurrentPropertyNames(["concentration","pressure"])
land.outputLand()

flora = Flora.Flora("testflora", "testconstants", land, "./")
flora.randomlyPlantTreesInRectangularDomain([1],["Avicennia"],land.bounding_box)

bettina = Bettina.Bettina("Testbettina", land, flora)

bettina.setConstantSubsurfaceProperties([0,1000,2000],[c,p])
#bettina.setVariableSubsurfaceProperties(files_t, working_directory)
bettina.evolveSystem()
"""
"""
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
        self.landMesh.fillRectangularGridWithVoxels()
        
        
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
land.outputLand("./")"""