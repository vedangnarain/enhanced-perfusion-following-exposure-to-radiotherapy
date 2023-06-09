#!/usr/bin/env python


"""Copyright (c) 2005-2017, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""


"""Script to manage wiki page generation
"""
import fnmatch
import os
import sys
import ntpath
import subprocess

if __name__ == '__main__':

    # Find all the tutorial files.
    search_paths = ["test/tutorials/", "test/python/tutorials/"]
    
    tutorial_files = []
    for eachSearchPath in search_paths:
        for root, dirs, files in os.walk(eachSearchPath):
            for file in files:
                if fnmatch.fnmatch(file, 'Test*LiteratePaper*') or fnmatch.fnmatch(file, 'Test*Tutorial*'):
                    if not fnmatch.fnmatch(file, '*.pyc'):
                        tutorial_files.append([eachSearchPath, file])
                        
    output_format = "jupyter"
    
    if output_format == "markdown":           
        # Generate the markdown for each
        for eachFile in tutorial_files:
            outfile = " doc/tutorials/" + os.path.splitext(ntpath.basename(eachFile[1]))[0] +"_md.md"
            inputfile = eachFile[0] + "/" + eachFile[1]
            if fnmatch.fnmatch(inputfile, '*.hpp'):
                launch_string = "infra/CreateMarkdownTutorial.py " + inputfile + outfile 
                os.system(launch_string)
    
    elif output_format == "jupyter": 
        # Generate the jupyter notebooks for each
        for eachFile in tutorial_files:
            outfile = " doc/tutorials/" + os.path.splitext(ntpath.basename(eachFile[1]))[0] +".ipynb"
            inputfile = eachFile[0] + "/" + eachFile[1]
            if not fnmatch.fnmatch(inputfile, '*.hpp'):
                print "Processing: ", inputfile
                launch_string = "infra/CreateJupyterNotebookTutorial.py " + inputfile + outfile 
                os.system(launch_string)
                
                subprocess.call("jupyter nbconvert " + outfile, shell=True)

