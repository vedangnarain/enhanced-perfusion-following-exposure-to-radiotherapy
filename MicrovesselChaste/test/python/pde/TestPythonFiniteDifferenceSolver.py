
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

import unittest
import math
import os
import vtk
import numpy as np
import chaste
import chaste.core
chaste.init()
import microvessel_chaste
import microvessel_chaste.geometry
import microvessel_chaste.mesh 
import microvessel_chaste.pde 
from microvessel_chaste.utility import * # bring in all units for convenience

class TestFiniteDifferenceSolver(microvessel_chaste.utility.AbstractMicrovesselTestSuite):
          
    def test_fixed_outer_boundary(self):
        
        file_handler = chaste.core.OutputFileHandler("Python/TestFiniteDifferenceSolver/test_fixed_outer_boundary")
        
        domain = microvessel_chaste.geometry.Part3()
        length_scale = 1.e-6*metres
        domain.AddCuboid(100.e-6*metres, 100.e-6*metres, 100.e-6*metres)
        
        grid = microvessel_chaste.mesh.RegularGrid3()
        grid.GenerateFromPart(domain, 10.e-6*metres)
         
        pde = microvessel_chaste.pde.DiscreteContinuumLinearEllipticPde3_3()
        pde.SetIsotropicDiffusionConstant(0.003*metre_squared_per_second)
        pde.SetContinuumLinearInUTerm(-1.0*per_second)
         
        bc = microvessel_chaste.pde.DiscreteContinuumBoundaryCondition3()
        bc.SetValue(30.0*mole_per_metre_cubed)
         
        solver = microvessel_chaste.pde.SimpleLinearEllipticFiniteDifferenceSolver3()
        solver.SetGrid(grid)
        solver.SetPde(pde)
        solver.AddBoundaryCondition(bc)
        solver.SetFileHandler(file_handler)
        solver.SetWriteSolution(True)
        solver.Solve()

if __name__ == '__main__':
    unittest.main()
        
        