#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Class managing the whola salt seltup. This class is containing 3 important ob-
jects: land, flora and the ogs-project information. The built-in functions are
designed to manage information transfer between bettina and ogs model
simulations.
@date: 2018 - Today
@author: jasper.bathmann@ufz.de
"""
import numpy as np
import sys
sys.path.append('./pymeshinteraction/')
import MeshPointFinder as MPF
import MeshInteractor

sys.path.append('./pybettina/')
import Bettina
import SortFileNames as SFN

sys.path.append('./pyogsproject/')
import OGSProject


class SaltSetup:
    def __init__(self, name, working_directory, land, flora, constant_density,
                 output_midstring):
        self.setup_name = name
        self.working_directory = working_directory
        self.boundary_surfaces = []
        self.land = land
        self.flora = flora
        self.output_midstring = output_midstring
        self.constant_density = constant_density
        self.createOGSProject(constant_density)
        self.initializeBettina()

    def setLand(self, land):
        self.land = land

    def setFlora(self, flora):
        self.flora = flora

    def initializeBettina(self):
        self.bettina = Bettina.Bettina(self.setup_name+"Bettina", self.land,
                                       self.flora)

    def createOGSProject(self, constant_density):
        self.ogsPrj = OGSProject.OGSProject(self.working_directory,
                                            self.setup_name + "_OGSproject")
        self.ogsPrj.setLandName(self.land.initial_mesh_name)
        self.ogsPrj.initializeProject(constant_density)
        self.ogsPrj.project.output_prefix = (self.setup_name +
                                             self.output_midstring)

    def setInitialConditionNames(self, p_ini, c_ini, q_ini):
        self.c_ini_name = c_ini
        self.p_ini_name = p_ini
        self.q_ini_name = q_ini
        self.ogsPrj.setInitialConditionName(p_ini, c_ini)

    def setVariableNames(self, p_var, c_var):
        self.c_var_name = c_var
        self.p_var_name = p_var
        self.ogsPrj.setVariableNames(p_var, c_var)
        self.land.setCurrentPropertyNames([c_var, p_var])

    def setInitialConditionName(self, p_ini, c_ini, q_ini):
        self.c_ini_name = c_ini
        self.p_ini_name = p_ini
        self.q_ini_name = q_ini
        self.ogsPrj.setInitialConditionName(p_ini, c_ini)

    def setQiniArray(self, q_ini):
        self.q_ini = q_ini

    def addBoundaryConditionForSurfaces(self, variable, *args):
        self.ogsPrj.createBoundaryConditionsFromList(self.boundary_surfaces,
                                                     variable, *args)

    def writeOgsProject(self):
        self.ogsPrj.writeProject()

    def createBoundarySurface(self, location):
        # This functions only works for land domains with the shape of regular
        # 3D boxes.
        land_mesh = self.land.initial_mesh
        land_grid = land_mesh.grid
        origin = [land_mesh.ox, land_mesh.oy, land_mesh.oz]
        deltas = [land_mesh.dx, land_mesh.dy, land_mesh.dz]
        steps = [land_mesh.extension_points_X, land_mesh.extension_points_Y,
                 land_mesh.extension_points_Z]
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
        if location == "bottom" or location == "top":
            steps[2] = 1
        if location == "top":
            origin[2] = origin[2] + land_mesh.lengthZ
        boundary_mesh_name = land_mesh.meshName + "_" + location + "Boundary"
        boundary_creator = MPF.MeshPointFinder(land_grid)
        points, ids = boundary_creator.findPointsOnPlane(steps[0], steps[1],
                                                         steps[2], deltas[0],
                                                         deltas[1], deltas[2],
                                                         origin[0], origin[1],
                                                         origin[2])
        temppoints = []
        for i in range(len(points)):
            point = points[i]
            x = (0+shift) % 3
            y = (1+shift) % 3
            z = (2+shift) % 3
            point = point[x], point[y], point[z]
            temppoints.append(point)
        boundary = MeshInteractor.MeshInteractor(boundary_mesh_name)
        boundary.CreateMeshFromPoints([points, ids])
        temp_boundary = MeshInteractor.MeshInteractor(boundary_mesh_name)
        temp_boundary.CreateMeshFromPoints([temppoints, ids])
        temp_boundary.CreateMultipleTriangles()
        cells = temp_boundary.grid.GetCells()
        if(land_grid.GetBounds()[-1] - land_grid.GetBounds()[-2] == 0):
            boundary.grid.SetCells(3, cells)
        else:
            boundary.grid.SetCells(5, cells)
        full_c_ini = land_grid.GetPointData().GetArray(self.c_ini_name)
        full_p_ini = land_grid.GetPointData().GetArray(self.p_ini_name)
        full_q_ini = land_grid.GetPointData().GetArray(self.q_ini_name)
        boundary_c_ini = np.zeros(len(ids))
        boundary_p_ini = np.zeros(len(ids))
        boundary_q_ini = np.zeros(len(ids))
        for i in range(len(ids)):
            iD = ids[i]
            boundary_c_ini[i] = full_c_ini.GetTuple(iD)[0]
            boundary_p_ini[i] = full_p_ini.GetTuple(iD)[0]
            boundary_q_ini[i] = full_q_ini.GetTuple(iD)[0]
        boundary.AddPropertyVector(np.array(boundary_c_ini), self.c_ini_name,
                                   "double")
        boundary.AddPropertyVector(np.array(boundary_p_ini), self.p_ini_name,
                                   "double")
        boundary.AddPropertyVector(np.array(boundary_q_ini), self.q_ini_name,
                                   "double")
        boundary.AddPropertyVector(0, "zero", "double")
        boundary.OutputMesh(self.working_directory)
        self.boundary_surfaces.append(boundary_mesh_name)

    def setTimeSteppingAndOutputLoopsForOgs(self, timerepeats, timedeltaTs,
                                            outputrepeats, outputdeltaN):
        self.ogsPrj.setTimeSteppingAndOutputLoops(timerepeats, timedeltaTs,
                                                  outputrepeats, outputdeltaN)

    def updateFloraBoundaryConditions(self):
        self.root_surfaces = self.flora.getAllRootNames()
        c1, c2, c3 = "coeff1", "coeff2", "coeff3"
        args = "NonuniformVariableDependentNeumann", "constant", c1, c2, c3
        self.ogsPrj.createTreeBoundaryConditionsFromList(self.root_surfaces,
                                                         "pressure", *args)

    def updateBoundaryConditions(self):
        self.ogsPrj.resetBoundaryConditions()
        self.addBoundaryConditionForSurfaces("pressure")
        self.addBoundaryConditionForSurfaces("concentration")
        self.updateFloraBoundaryConditions()

    def runOGS(self):
        self.ogsPrj.runOgs()

    def setOgsTiniTend(self, t_ini, t_end):
        self.ogsPrj.setTiniTend(t_ini, t_end)

    def progressBettina(self, files_t):
        self.bettina.setVariableSubsurfaceProperties(files_t,
                                                     self.working_directory)
        self.bettina.evolveSystem()

    def updateModel(self):
        self.ogsPrj.setLandName(self.land.initial_mesh_name)
        self.ogsPrj.initializeProject(self.constant_density)
        self.ogsPrj.project.output_prefix = self.setup_name + \
            self.output_midstring

    def createMeshCollection(self, prefix, postfix):
        self.file_reader_meshes = SFN.ReadAndSortFileNames(
                self.working_directory, prefix, postfix)

    def createFloraCollection(self, prefix, postfix):
        self.file_reader_flora = SFN.ReadAndSortFileNames(
                self.working_directory, prefix, postfix)

    def updateFloraCollection(self):
        self.file_reader_flora.createPvDFile(
                self.setup_name + "_flora_meshes" + ".pvd")
        files = self.file_reader_flora.getSortedFiles()
        self.file_reader_flora.addMeshesToPvdFile(files)

    def updateMeshCollection(self):
        self.file_reader_meshes.createPvDFile(
                self.setup_name + "_land_meshes" + ".pvd")
        files = self.file_reader_meshes.getSortedFiles()
        rm_files = []
        for file in files:
            if("_pcs_0_ts_0" in file):
                rm_files.append(file)
        for rm in rm_files:
            files.remove(rm)
        self.file_reader_meshes.addMeshesToPvdFile(files)

    def createTreeCollection(self, species_list, postfix):
        self.file_reader_trees_list = []
        for species in species_list:
            self.file_reader_trees_list.append(
                SFN.ReadAndSortFileNames(self.working_directory,
                                         species, postfix))

    def updateTreeCollection(self):
        i = 1
        for file_reader_trees in self.file_reader_trees_list:
            file_reader_trees.createPvDFile(
                self.setup_name + "_tree_meshes_" + file_reader_trees.prefix
                + ".pvd")
            i += 1
            files = file_reader_trees.getSortedFiles()
            file_reader_trees.addMeshesToPvdFile(files)

    def readAndPassTMeshFileNames(self, t_ini, t_end):
        t_files = self.file_reader_meshes.getFilesInTimeIntervall(t_ini, t_end)
        self.progressBettina(t_files)
        self.updateMeshCollection()
        self.updateTreeCollection()
        self.updateFloraCollection()

    def setAndRunOgs(self, t_ini, t_end, timerepeats, timedeltaTs,
                     outputrepeats, outputdeltaN):
        self.setOgsTiniTend(t_ini, t_end)
        self.setTimeSteppingAndOutputLoopsForOgs(timerepeats, timedeltaTs,
                                                 outputrepeats, outputdeltaN)

        self.writeOgsProject()
        self.runOGS()
