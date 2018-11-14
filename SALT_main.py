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


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 14:27:35 2018

@author: bathmann
"""
#sys.path.append('../pySALT/FileOperations')
import Land
import Flora
import Bettina
import pyOgsProject
import SortFileNames as SFN
#import SortFilenames as SFn
#TODO: Check, where the script file is gone. Maybe it is saved on the ssd...
import numpy as np
##Dummy files for Jasper to test interaction with OGS:
class OGSproject:
    def __init__(self,project_dir, project_name):
        self.project_dir = project_dir
        self.project_name = project_name
        
    def setProjectdir(self, project_dir):
        self.project_dir = project_dir
    
    def setProjectName(self, project_name):
        self.project_name = project_name
        
    def setLandName(self, land_name):
        self.land_name = land_name
        
    def iniParameters(self):
        self.parameter_names = []
        self.parameter_types = []
        self.parameter_values = []
        
    def addParameter(self, name, tYpe, value):
        self.parameter_names.append(name)
        self.parameter_types.append(tYpe)
        self.parameter_values.append(value)

    def initializeProject(self):
        self.project = pyOgsProject.GenerateProject(self.project_dir+self.project_name + ".prj")
        self.project.setMesh(self.land_name+".vtu")
        self.project.setStandardProcessInformation()
        self.project.setStandartTimeLoop()
        self.project.setStandartParameters()
        self.project.setStandardDensityModel()
        self.project.setStandardNonlinearSolvers()
        self.project.resetInitialConditions("p_ini","c_ini")
        self.project.processspeci_bo_force = "0 0 -"+str(self.g)
        self.project.resetBoundaryConditions()
working_directory = "/home/bathmann/Dokumente/UFZ/code/pySALT/testcases/TestCase_SteadyState/"
prefix = "Output_HC_Testcase-SteadyState_pcs"
postfix = ".vtu"
sFN = SFN.ReadAndSortFileNames(working_directory, prefix, postfix)
all_files_sorted = sFN.getSortedFiles()
files_t = sFN.getFilesInTimeIntervall(100000, 800000, all_files_sorted)
print(files_t)

"""
working_directory = "/home/bathmann/Dokumente/UFZ/code/pySALT/testcases/TestCase_SteadyState/"
files_t = ["Output_HC_Testcase-SteadyState_pcs_0_ts_8900_t_8900000.000000.vtu",
           "Output_HC_Testcase-SteadyState_pcs_0_ts_9000_t_9000000.000000.vtu",
           "Output_HC_Testcase-SteadyState_pcs_0_ts_9100_t_9100000.000000.vtu",
           "Output_HC_Testcase-SteadyState_pcs_0_ts_9200_t_9200000.000000.vtu"]


working_directory = "./"
land = Land.Land("testland", working_directory)
land.create3DRectangularLand("testlandmesh", 0, 0, -1, 10, 10, 1, 11, 11, 3)
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