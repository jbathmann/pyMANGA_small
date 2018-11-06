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
import MeshInteractor
import vtk
import numpy as np

def calculatePressure(g,conc_difference_ratio,mesh):
    
    array = vtk.vtkDoubleArray()
    array.SetName("c_ini")
    array2 = vtk.vtkDoubleArray()
    array2.SetName("p_ini")
    for i in range(mesh.GetNumberOfPoints()):
        point = mesh.GetPoints().GetPoint(i)
        cvalue = 0
        pvalue = 0
        #[8]: (26.7-25.5)/52#=0.023076923076923064
        #[8]: (26.2-25.5)/52#=0.013461538461538448
        #(26.55-25.5)/53 Out[3]: 0.019811320754716994

        hoz_gradient = 0.023076923076923064
        if(point[0] == 0):
            cvalue = 0.0371
            
        patmo = -1000*(1+conc_difference_ratio*cvalue)*g*(point[1]-0.49056603773584906*.53)
        px = g*1000*(1+conc_difference_ratio*.0371*(.53-point[0])/.53)*hoz_gradient*point[0]

        #if(point[0] == 0):#f(patmo>(-1000*1*g*(point[1]-headleft)+g*1000*(1)*hoz_gradient*(.53))):
        #    cvalue=0

        pvalue = patmo+px
        #if(point[2]>-.1):
        #    cvalue=0
        array.InsertNextTuple1(cvalue)
        array2.InsertNextTuple1(pvalue)
    return array, array2

        
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

main_origin_x = 0.04
main_origin_y = 0.16
partsx = 50
partsz = 11
parts_per_m = 100 
lx0, lx1, lx2, lx3 = 0.00, 0.00, 0.015, (0.015 )#+ 0.115)
#lz0, lz1, lz2, lz3 = 0.01, 0.01, 0.00, 0.01
#ux0, ux1, ux2, ux3 = 0.01, 0.01, 0.01, 0.01 
uz0, uz1, uz2, uz3 = 0.000, 0.00, 0.085, (0.015 )#+ 0.115) 
lx =[lx0,lx1,lx2,lx3] 
#lz =[lz0,lz1,lz2,lz3] 
#ux =[ux0,ux1,ux2,ux3] 
uz =[uz0,uz1,uz2,uz3] 
factor = np.array([16,8,4,2]) * parts_per_m

#Setup of land mesh and project
land = Land.Land(landName,l, parts, geometricalRatio, horizontal_gradient, g, vertical_gradient, useOgs)
land.setProjectdir(projectdir)
land.setMeanC(meanc, random)
land.landMesh.linear=False
#land.lz=.26
origin_x = 0.00
origin_z = 0.00
partsx = int((land.lx-sum(lx)) * parts_per_m + 1)
partsz = int((land.lz-sum(uz)) * parts_per_m + 1)
land.landMesh.Create2DMeshFromGeometry(partsx, sum(lx), land.lx-sum(lx), partsz, sum(uz), land.lz-sum(uz),
                               land.horizontal_gradient, land.g,  
                               boundary=False)
unterer_bereich = MeshInteractor.MeshInteractor("unterer_rechts")
unterer_bereich.Create2DMeshFromGeometry(38, .53-.37, .37 , int(sum(uz)*parts_per_m+1), 0, sum(uz), 
                                  land.horizontal_gradient , land.g,boundary = False)
land.landMesh.CombineMultipleGrids([unterer_bereich.pd])
unterer_bereich.Create2DMeshFromGeometry(7, .53-.40, .03 , int(sum(uz)*parts_per_m*2+1), 0, sum(uz), 
                                  land.horizontal_gradient , land.g,boundary = False)
land.landMesh.CombineMultipleGrids([unterer_bereich.pd])
for i in range(4):
    j = 2**(i+1)
    lxi = lx[i]
    #lzi = lz[i]
    #uxi = ux[i]
    uzi = uz[i]

    factori = factor[i]
    partsxL = int(lxi * factori + 1)
    partszL = int((land.lz - origin_z) * factori + 1)
    partsxU = int((land.lx-.4- origin_x - lxi) * factori +1)
    partszU = int(uzi * factori + 1)
    if(lxi!=0):
        linker_bereich = MeshInteractor.MeshInteractor("linker_bereich")
        linker_bereich.Create2DMeshFromGeometry(partsxL, origin_x, lxi, partszL, origin_z, land.lz-origin_z, 
                                          land.horizontal_gradient , land.g,boundary = False)
        land.landMesh.CombineMultipleGrids([linker_bereich.pd])
    if(uzi!=0):
        unterer_bereich = MeshInteractor.MeshInteractor("unterer_bereich")
        unterer_bereich.Create2DMeshFromGeometry(partsxU, origin_x + lxi, land.lx-.4 - lxi - origin_x, partszU, origin_z, uzi, 
                                          land.horizontal_gradient , land.g,boundary = False)
        land.landMesh.CombineMultipleGrids([unterer_bereich.pd])
    origin_x += lxi
    origin_z += uzi





land.initializeProject()
landmesh = land.landMesh.pd
cland, pland = calculatePressure(g,conc_difference_ratio,landmesh)
landmesh.GetPointData().RemoveArray("c_ini")
landmesh.GetPointData().RemoveArray("p_ini")
landmesh.GetPointData().AddArray(cland)
landmesh.GetPointData().AddArray(pland)
land.landMesh.pd = landmesh

land.project.setTimeSteppingAndOutputLoops(timerepeats, timedeltaTs, outputrepeats, outputdeltaN)

Dm =0#1e-9 
beta_t = 5e-4
land.project.parameter_values[1]=str(Dm)
land.project.parameter_values[5]=str(beta_t)
land.project.parameter_values[4]="5e-3" #beta_l
land.project.parameter_values[-1]="1.2388e-9 0 0 1.2388e-9"#kappa # Goswami paper says: 1.2388e-9
#self.parameter_values = [rho_v, Dm_v, retardation_v, decay_v, beta_l_v, beta_t_v, c_ini_v, 
#                                 p_ini_v, constant_porosity_parameter_v, kappa1_v]
land.project.parameter_values[-2]="0.385"
land.project.convergence_criterion_reltols = "5e-6 5e-6"
land.project.processspeci_bo_force = "0 -"+str(g)
land.project.densityModel = "ConcentrationDependent"
land.project.solver_linear_solver_lis = "-i bicgstab -p ilut -tol 1e-12 -maxiter 200000"
land.project.solver_linear_solver_petsc_parameters = "-hc_ksp_type bcgs -hc_pc_type bjacobi -hc_ksp_rtol 1e-10 -hc_ksp_max_it 20000"
land.project.solver_linear_solver_eigen_error_tol = "1e-11"
#land.initializeHorizontalFlow()
land.project.output_prefix=prefix
darcyvel=land.getHotizontalDarcyVelBoundary()
#omega = calculateAnalyticalSolution(darcyvel, Dm, beta_t, l, meanc, land.landMesh.pd)
#landMeshEdit.AddPropertyVector(omega,"Analytical Steady State Solution","double")

#Initializes horizontal flow by means of boundary conditions
links = MeshInteractor.MeshInteractor(land.landName+"_leftBoundary")
links.cvalue = land.landMesh.cvalue
links.random = land.landMesh.random
links.linear = land.landMesh.linear
links.pd = land.landMesh.GetBoundaryOf2DMesh("left")

links.AddPropertyVector(1000*darcyvel,"darcy_ini","double")

clinks, plinks = calculatePressure(g,conc_difference_ratio,links.pd)
links.pd.GetPointData().RemoveArray("c_ini")
links.pd.GetPointData().RemoveArray("p_ini")
links.pd.GetPointData().AddArray(clinks)
links.pd.GetPointData().AddArray(plinks)


links.OutputMesh(land.projectdir)
rechts = MeshInteractor.MeshInteractor(land.landName+"_rightBoundary")
rechts.cvalue = land.landMesh.cvalue
rechts.random = land.landMesh.random
rechts.linear = land.landMesh.linear
rechts.pd = land.landMesh.GetBoundaryOf2DMesh("right")

#rechts.Create2DMeshFromGeometry(partsx, land.lx, partsz, land.lz, 
#                                  land.horizontal_gradient , land.g,boundary = "right")

rechts.AddPropertyVector(-1000*darcyvel,"darcy_ini","double")
rechts.AddPropertyVector(-1000*darcyvel*1,"constant","double")
rechts.AddPropertyVector(0,"prefac1","double")
rechts.AddPropertyVector(1000*darcyvel/target_c,"prefac2","double")
crechts, prechts = calculatePressure(g,conc_difference_ratio,rechts.pd)
rechts.pd.GetPointData().RemoveArray("c_ini")
rechts.pd.GetPointData().RemoveArray("p_ini")
rechts.pd.GetPointData().AddArray(crechts)
rechts.pd.GetPointData().AddArray(prechts)


rechts.OutputMesh(land.projectdir)


"""
oben = MeshInteractor.MeshInteractor(land.landName+"_topBoundary")
oben.cvalue = land.landMesh.cvalue
oben.random = land.landMesh.random
oben.linear = land.landMesh.linear
oben.Create2DMeshFromGeometry(land.partsx, land.lx, land.partsz, land.lz, 
                                  land.horizontal_gradient , land.g,boundary = "top")

oben.AddPropertyVector(-1000*darcyvel,"darcy_ini","double")
oben.AddPropertyVector(-1000*darcyvel*1,"constant","double")
oben.AddPropertyVector(0,"prefac1","double")
oben.AddPropertyVector(1000*darcyvel/target_c,"prefac2","double")
coben, poben = calculatePressure(g,conc_difference_ratio,oben.pd)
oben.pd.GetPointData().RemoveArray("c_ini")
oben.pd.GetPointData().RemoveArray("p_ini")
oben.pd.GetPointData().AddArray(coben)
oben.pd.GetPointData().AddArray(poben)


oben.OutputMesh(land.projectdir)"""
land.project.createBoundaryCondition("p","NonuniformDirichlet", "geometry","left",land.landName+"_leftBoundary","p_ini")
land.project.createBoundaryCondition("p","NonuniformDirichlet", "geometry","right",land.landName+"_rightBoundary","p_ini")
land.project.createBoundaryCondition("c","NonuniformDirichlet", "geometry","right",land.landName+"_leftBoundary","c_ini")
land.project.createBoundaryCondition("c","NonuniformDirichlet", "geometry","right",land.landName+"_rightBoundary","c_ini")
#land.project.createBoundaryCondition("c","NonuniformDirichlet", "geometry","right",land.landName+"_topBoundary","c_ini")


#land.project.createBoundaryCondition("p","NonuniformVariableDependantNeumann", land.landName+"_rightBoundary","constant", "prefac1", "prefac2" )
#Initialization of 1 Tree
#newTree= Tree.Tree(coordx, coordy, land, "Avicennia", 0)
#newTree.plantTree()
#flora.trees.append(newTree)
#flora.createTreeBCs()
land.project.time_stepping_t_end = str(t_end)
land.landMesh.OutputMesh(projectdir)
#        self.parameter_values = [rho_v, Dm_v, retardation_v, decay_v, beta_l_v, beta_t_v, c_ini_v, 
#                                 p_ini_v, constant_porosity_parameter_v, kappa1_v]

land.project.writeProjectFile()    
land.getSalinitiesAndPressures()    