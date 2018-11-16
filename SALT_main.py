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
    def __init__(self, name, working_directory, land, flora):
        self.setup_name = name
        self.working_directory = working_directory
        self.boundary_surfaces = []
        self.land = land
        self.flora = flora
        self.createOGSProject()
        self.initializeBettina()
        
    def setLand(self, land):
        self.land = land
    
    def setFlora(self, flora):
        self.flora = flora
    
    def initializeBettina(self):
        self.bettina = Bettina.Bettina(self.setup_name+"Bettina", self.land, self.flora)
        
    def createOGSProject(self):
        self.ogsPrj = OGSProject.OGSProject(self.working_directory, self.setup_name + "_OGSproject")
        self.ogsPrj.setLandName(self.land.initial_mesh_name)
        self.ogsPrj.initializeProject()
        self.ogsPrj.project.output_prefix = self.setup_name + "_Output"
        
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
        full_p_ini = land_grid.GetPointData().GetArray(self.p_ini_name)
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

    def updateFloraBoundaryConditions(self):
        self.root_surfaces = flora.getAllRootNames()
        args = "NonuniformVariableDependentNeumann", "constant", "coeff1", "coeff2", "coeff3"
        self.ogsPrj.createTreeBoundaryConditionsFromList(self.root_surfaces, "pressure", *args)
    
    def updateBoundaryConditions(self):
        self.ogsPrj.resetBoundaryConditions()
        self.addBoundaryConditionForSurfaces("pressure", "NonuniformDirichlet")
        self.addBoundaryConditionForSurfaces("concentration", "NonuniformDirichlet")
        self.updateFloraBoundaryConditions()
        
    def runOGS(self):
        self.ogsPrj.runOgs()
        
    def setOgsTiniTend(self, t_ini, t_end):
        self.ogsPrj.setTiniTend(t_ini, t_end)
        
    def progressBettina(self, files_t):
        self.bettina.setVariableSubsurfaceProperties(files_t, self.working_directory)
        self.bettina.evolveSystem()
        
    def updateModel(self):
        self.ogsPrj.setLandName(self.land.initial_mesh_name)
        self.ogsPrj.initializeProject()
        self.ogsPrj.project.output_prefix = self.setup_name + "_Output"
working_directory = "./testruns/"
setup_name ="testrun"
prefix = setup_name + "_Output_pcs"
postfix = ".vtu"
land = Land.Land("testlandmesh", working_directory)
land.create3DRectangularLand("testlandmesh", 0, 0, -1, 100, 100, 1, 101, 101, 3)
n = land.initial_mesh.grid.GetNumberOfPoints()
c,p , iD = [], [], []
for i in range(n):
    point = land.initial_mesh.grid.GetPoints().GetPoint(i)
    c.append(.035)
    p.append(-(1000+35*.7)*9.81*(point[2]-1e-2*point[0]))
    iD.append(i)
land.setCIniPIniAndNodeIds(["c_ini", "p_ini"], "bulk_node_ids", [np.array(c), np.array(p)], np.array(iD))
land.outputLand()

flora = Flora.Flora("testflora", "testconstants", land, working_directory)
flora.randomlyPlantTreesInRectangularDomain([10],["Avicennia"],land.bounding_box)

model = SaltSetup(setup_name, working_directory, land, flora)
model.setVariableNames("pressure","concentration")
model.setInitialConditionName("p_ini", "c_ini")
model.createBoundarySurface("left")
model.createBoundarySurface("right")
model.updateBoundaryConditions()


file_reader = SFN.ReadAndSortFileNames(working_directory, prefix, postfix)
#1/2 Jahr = 15778800.0 Sekunden
file_reader.createPvDFile( "land_meshes.pvd")
dt = 15778800.0
for i in range(3):
    model.updateBoundaryConditions()
    t_ini = i * dt
    t_end = (i + 1) * dt
    #1 Tag = [100, 50, 59, 23 mit [1e-1, 1e0, 60, 3600, ]
    timerepeats = [100, 50, 59, 23, 29*24, 30*24*5]
    timedeltaTs = [1e-1, 1e0, 60, 3600, 3600, 3600]
    outputrepeats = [1, 29, 30*5]
    outputdeltaN = [232, 24, 24]
    model.setOgsTiniTend(t_ini, t_end)
    model.setTimeSteppingAndOutputLoopsForOgs(timerepeats, timedeltaTs, outputrepeats, outputdeltaN)
    
    model.writeOgsProject()
    model.runOGS()
    t_files = file_reader.getFilesInTimeIntervall(t_ini, t_end)[1:]
    print(t_files)
    file_reader.addLandMeshesToPvdFile(t_files)
    model.progressBettina(t_files)
    model.updateModel()
file_reader.finishPvDFile()
    
