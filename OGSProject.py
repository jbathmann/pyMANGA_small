#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pyogsproject import pyOgsProject
import os


class OGSProject:
    ## Tool to organise ogs prj files.
    #  @VAR project_dir: path to folder, where project files are saved \n
    #  @VAR project_name: project file name
    #  @author: bathmann

    def __init__(self, project_dir, project_name):
        self.project_dir = project_dir
        self.project_name = project_name
        self.g = 9.81
        self.c_ini_name = "c_ini"
        self.p_ini_name = "p_ini"
        self.q_ini_name = "q_ini"

    def setProjectdir(self, project_dir):
        self.project_dir = project_dir

    def setProjectName(self, project_name):
        self.project_name = project_name

    def setLandName(self, land_name):
        self.land_name = land_name

    def setGravitation(self, value):
        self.g = value

    def setInitialConditionName(self, p_ini, c_ini):
        self.c_ini_name = c_ini
        self.p_ini_name = p_ini
        self.project.resetInitialConditions(self.p_ini_name, self.c_ini_name)

    def getInitialConditionName(self):
        return self.c_ini_name, self.p_ini_name

    def setVariableNames(self, p_var, c_var):
        self.c_var_name = c_var
        self.p_var_name = p_var

    def getVariableNames(self):
        return self.c_var_name, self.p_var_name

    def iniParameters(self):
        self.parameter_names = []
        self.parameter_types = []
        self.parameter_values = []

    def addParameter(self, name, tYpe, value):
        self.parameter_names.append(name)
        self.parameter_types.append(tYpe)
        self.parameter_values.append(value)

    def initializeProject(self, constant_density):
        self.project = pyOgsProject.GenerateProject(self.project_dir +
                                                    self.project_name + ".prj")
        self.project.setMesh(self.land_name+".vtu")
        self.project.setStandardDensityModel()
        if constant_density:
            self.project.densityModel = "Constant"
        self.project.setStandardProcessInformation()
        self.project.setStandartAdaptiveTimeLoop()
        self.project.setStandartParameters()

        self.project.setStandardNonlinearSolvers()
        self.project.resetInitialConditions(self.p_ini_name, self.c_ini_name)
        self.project.processspeci_bo_force = "0 0 -"+str(self.g)
        self.project.resetBoundaryConditions()

    def writeProject(self):
        self.project.writeProjectFile()

    def createBoundaryConditionsFromList(self, mesh_list, variable):
        for mesh in mesh_list:
            if "right" in mesh:
                if variable == self.p_var_name:
                    self.project.createBoundaryCondition(
                            "p", "Dirichlet",
                            mesh, self.p_ini_name)
                if variable == self.c_var_name:
                    self.project.createBoundaryCondition(
                            "c", "Dirichlet",
                            mesh, self.c_ini_name)
            elif "left" in mesh:
                if variable == self.c_var_name:
                    self.project.createBoundaryCondition(
                            "c", "Dirichlet",
                            mesh, self.c_ini_name)
                if variable == self.p_var_name:
                    self.project.createBoundaryCondition(
                            "p", "Neumann",
                            mesh, self.q_ini_name)

    def createTreeBoundaryConditionsFromList(
            self, mesh_list, variable, tYpe, constant, coeff1, coeff2, coeff3):
        for mesh in mesh_list:
            if variable == self.p_var_name:
                self.project.createBoundaryCondition(
                        "p", tYpe, mesh, constant, coeff1, coeff2, coeff3)

    def setTimeSteppingOutputLoops(self, timerepeats, timedeltaTs):
        self.project.setFixedTimeStepping(timerepeats, timedeltaTs)

    def setOutputLoops(self, outputrepeats, outputdeltaN):
        self.project.setOutputLoops(outputrepeats, outputdeltaN)

    def setFixedOutputTimes(self, fixed_output_times):
        self.project.setFixedOutputTimes(fixed_output_times)

    def resetBoundaryConditions(self):
        self.project.resetBoundaryConditions()

    def setTiniTend(self, t_ini, t_end):
        self.project.time_stepping_t_ini = str(t_ini)
        self.project.time_stepping_t_end = str(t_end)

    def runOgs(self):
        os.system(
                "ogs " + self.project_dir + self.project_name +
                ".prj -o " + self.project_dir)
