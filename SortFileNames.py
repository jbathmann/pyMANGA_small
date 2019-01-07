#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tool to identify specific files in a folder.
working_directory: path to folder, where the files are stored
prefix: characteristic string contained in all files of interest
postfix: file ending, e.g. ".vtu", ".pvd"...
@date: 2018-Today
@author: jasper.bathmann@ufz.de
"""
import os
import re


class ReadAndSortFileNames:
    def __init__(self, working_directory, prefix, postfix):
        self.working_directory = working_directory
        self.prefix = prefix
        self.postfix = postfix
        self.numbers = re.compile(r'(\d+)')

    def numericalSort(self, value):
        # function to sort given list of strings numerically, e.g. "file1",
        # "file2",...
        parts = self.numbers.split(value)
        parts[1::2] = map(int, parts[1::2])
        return parts

    def getSortedFiles(self):
        # returns all filenames containing the prefix and the postfix numerical
        # sorted
        relevant_files = []
        files = (os.listdir(self.working_directory))
        for file in files:
            if self.postfix in file and self.prefix in file:
                relevant_files.append(file)
        sorted_relevant_files = sorted(relevant_files, key=self.numericalSort)
        return sorted_relevant_files

    def getFilesInTimeIntervall(self, t_begin, t_end):
        # returns all filenames containing the prefix and the postfix and a
        # number large than t_begin and smaller than t_end in numerical
        # sorted order
        relevant_files = []
        files = self.getSortedFiles()
        for file in files:
            time = float(file.strip(self.postfix).split("_t_")[-1])
            if time == t_begin and "_ts_0_t_" in file:
                relevant_files.append(file)
            if t_begin < time <= t_end:
                relevant_files.append(file)
        return relevant_files

    def createPvDFile(self, pvd_file_name):
        # creates an empty pvd file with given pvd_file_name
        self.pvd_file_name = self.working_directory + pvd_file_name
        file = open(self.pvd_file_name, "w")
        file.write("""<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">
  <Collection>
                   \n""")
        file.close()

    def addMeshesToPvdFile(self, meshes):
        # adds a list of mesh files to the empty pvd file created in
        # self.createPvDFile
        tempfile = open(self.pvd_file_name, "r")
        lines = tempfile.readlines()
        if(len(lines) > 4):
            lines = lines[:-2]
        tempfile.close()
        file = open(self.pvd_file_name, "w")

        for line in lines:
            file.write(line)
        for mesh in meshes:
            t = float(mesh.strip(".vtu").split("_t_")[-1])
            file.write('        <DataSet timestep="' + str(t) +
                       '" group="" part="0" file="' + mesh + '"/>\n')
        file.write("""  </Collection>
</VTKFile>""")
        file.close()

    def finishPvDFile(self):
        # adds the closing lines to a given pvd file.
        file = open(self.pvd_file_name, "a")
        file.write("""  </Collection>
</VTKFile>""")
        file.close()
