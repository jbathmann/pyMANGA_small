#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from pymeshinteraction import *
from pybettina import *
import OGSProject
import os

class SaltSetup:

    ## Class managing the whole salt seltup. This class is containing 3
    #  important objects: land, flora and the ogs-project information. The
    #  built-in functions are designed to manage information transfer between
    #  bettina and ogs model simulations.
    #  @date: 2018 - Today
    #  @author: jasper.bathmann@ufz.de

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
        steps = [land_mesh.extension_points_X, land_mesh.extension_points_Y,
                 land_mesh.extension_points_Z]
        shift = 0
        if location == "left" or location == "right":
            steps[0] = 1
            shift = 1
        if location == "right":
            origin[0] = origin[0] + land_mesh.length_x
        if location == "front" or location == "back":
            steps[1] = 1
            shift = 2
        if location == "back":
            origin[1] = origin[1] + land_mesh.length_y
        if location == "bottom" or location == "top":
            steps[2] = 1
        if location == "top":
            origin[2] = origin[2] + land_mesh.length_z
        boundary_mesh_name = land_mesh.meshName + "_" + location + "Boundary"
        boundary_creator = MeshPointFinder.MeshPointFinder(land_grid, land_mesh.z)
        points, ids = boundary_creator.findPointsOnPlane(steps[0], steps[1],
                                                         steps[2],
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
        if location == "right":
            os.system(
                "ExtractSurface -x -1 -y 0 -z 0 -a 0 -o " + self.working_directory + boundary_mesh_name + ".vtu -i " + self.working_directory + self.land.initial_mesh_name +
                ".vtu")

            cell_data = boundary.readMesh(self.working_directory, boundary_mesh_name)
            boundary.setTempMeshAsMainMesh()
            boundary.readMesh(self.working_directory, self.land.initial_mesh_name)
            boundary.resampleDataset()
            boundary.grid.GetCellData().AddArray(cell_data.GetArray("bulk_element_ids"))
            boundary.grid.GetCellData().AddArray(cell_data.GetArray("bulk_face_ids"))
            boundary.outputMesh(self.working_directory)
        elif location == "top":
            os.system(
                "ExtractSurface -x -1 -y 0 -z 0 -a 30 -o " + self.working_directory + boundary_mesh_name + ".vtu -i " + self.working_directory + self.land.initial_mesh_name +
                ".vtu")

            cell_data = boundary.readMesh(self.working_directory, boundary_mesh_name)
            boundary.setTempMeshAsMainMesh()
            boundary.readMesh(self.working_directory, self.land.initial_mesh_name)
            boundary.resampleDataset()
            boundary.grid.GetCellData().AddArray(cell_data.GetArray("bulk_element_ids"))
            boundary.grid.GetCellData().AddArray(cell_data.GetArray("bulk_face_ids"))
            boundary.outputMesh(self.working_directory)
        else:
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
            boundary.addPropertyVector(np.array(boundary_c_ini), self.c_ini_name,
                                       "double")
            boundary.addPropertyVector(np.array(boundary_p_ini), self.p_ini_name,
                                       "double")
            boundary.addPropertyVector(np.array(boundary_q_ini), self.q_ini_name,
                                       "double")
        boundary.addPropertyVector(1, "one", "double")

        boundary.outputMesh(self.working_directory)
        self.boundary_surfaces.append(boundary_mesh_name)

    def updateFloraBoundaryConditions(self):
        self.root_surfaces = self.flora.getAllRootNames()
        constant, c1, c2, c3 = "constant", "coeff1", "coeff2", "coeff3"
        args = "VariableDependentNeumann", constant, c1, c2, c3
        self.ogsPrj.createTreeBoundaryConditionsFromList(self.root_surfaces,
                                                         "pressure", *args)

    def updateBoundaryConditions(self):
        self.ogsPrj.resetBoundaryConditions()
        self.boundary_surfaces = []
        self.createBoundarySurface("left")
        self.createBoundarySurface("right")
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
        kappa_t_old = self.ogsPrj.project.parameter_values[-1]
        self.ogsPrj.initializeProject(self.constant_density)
        self.ogsPrj.project.parameter_values[-1] = kappa_t_old
        self.ogsPrj.project.output_prefix = self.setup_name + \
            self.output_midstring

    def createMeshCollection(self, prefix, postfix):
        self.file_reader_meshes = SortFileNames.ReadAndSortFileNames(
                self.working_directory, prefix, postfix)

    def createFloraCollection(self, prefix, postfix):
        self.bettina.createFloraCollection(prefix, postfix)

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
        self.bettina.createTreeCollection(species_list, postfix)

    def readAndPassTMeshFileNames(self, t_ini, t_end):
        t_files = self.file_reader_meshes.getFilesInTimeIntervall(t_ini, t_end)
        self.progressBettina(t_files)
        self.updateMeshCollection()

    def setAndRunOgs(self, t_ini, t_end, outputrepeats, outputdeltaN,
                     fixed_output_times):
        self.setOgsTiniTend(t_ini, t_end)
        self.ogsPrj.setOutputLoops(outputrepeats, outputdeltaN)
        for i in range(len(fixed_output_times)):
            fixed_output_times[i] += t_ini
        print(fixed_output_times)
        self.ogsPrj.setFixedOutputTimes(fixed_output_times)
        self.writeOgsProject()
        self.runOGS()
