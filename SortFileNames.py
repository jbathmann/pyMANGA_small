#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 08:38:51 2018

@author: bathmann
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
        parts = self.numbers.split(value)
        parts[1::2] = map(int, parts[1::2])
        return parts
    
    def getSortedFiles(self):
        relevant_files = []
        files = (os.listdir(self.working_directory))
        for file in files:
            if self.postfix in file and self.prefix in file:
                relevant_files.append(file)
        sorted_relevant_files = sorted(relevant_files, key=self.numericalSort)
        return sorted_relevant_files
    
    def getFilesInTimeIntervall(self, t_begin, t_end):
        relevant_files = []
        files = self.getSortedFiles()
        for file in files:
            time = float(file.strip(self.postfix).split("_t_")[-1])
            if t_begin <= time <= t_end:
                relevant_files.append(file)
        return relevant_files

    def createPvDFile(self, pvd_file_name):
        self.pvd_file_name = self.working_directory + pvd_file_name
        file = open(self.pvd_file_name, "w")
        file.write("""<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">
  <Collection>
                   \n""")
        file.close()
        
    def addLandMeshesToPvdFile(self, meshes):
        tempfile = open(self.pvd_file_name, "r")
        lines = tempfile.readlines()
        if(len(lines)>4):
            lines = lines[:-2]
        tempfile.close()
        file = open(self.pvd_file_name, "w")

        for line in lines:
            file.write(line)
        for mesh in meshes[:-1]:
            t = float(mesh.strip(".vtu").split("_t_")[-1])
            file.write('        <DataSet timestep="' + str(t) + '" group="" part="0" file="'+ mesh + '"/>\n')
        file.write("""  </Collection>
</VTKFile>""")
        file.close()
        
    def finishPvDFile(self):
        file = open(self.pvd_file_name, "a")
        file.write("""  </Collection>
</VTKFile>""")
        file.close()