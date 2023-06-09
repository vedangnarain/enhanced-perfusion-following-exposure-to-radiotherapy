/*

Copyright (c) 2005-2017, University of Oxford.
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

*/

#ifndef TESTSOLVINGPDESLITERATEPAPER_HPP_
#define TESTSOLVINGPDESLITERATEPAPER_HPP_

/* # Solving PDEs in the Microvessel Project
 * This tutorial demonstrates methods for solving PDEs in the Microvessel Project. It is noted
 * that the way to set up PDEs differs from that of Cell Based Chaste, although the same solver
 * can be used behind the scenes.
 *
 * The following is covered:
 * * Solving a linear reaction-diffusion PDE with finite differences.
 * * Solving a linear reaction-diffusion PDE with finite finite elements and discrete sinks and sources.
 * * Solving a non-linear reaction-diffusion PDE with finite differences.
 * * Solving a linear reaction-diffusion PDE using Green's functions.
 * * Interacting with regular grids and finite element meshes.
 *
 * # The Test
 * Start by introducing the necessary header files. The first contain functionality for setting up unit tests,
 * smart pointer tools and output management.
 */
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "SmartPointers.hpp"
#include "OutputFileHandler.hpp"
#include "FileFinder.hpp"
/*
 * Dimensional analysis.
 */
#include "Vertex.hpp"
#include "UnitCollection.hpp"
#include "Owen11Parameters.hpp"
#include "GenericParameters.hpp"
#include "ParameterCollection.hpp"
#include "BaseUnits.hpp"
/*
 * Geometry tools.
 */
#include "MappableGridGenerator.hpp"
#include "Part.hpp"
/*
 * Vessel networks.
 */
#include "VesselNetwork.hpp"
#include "VesselNetworkGenerator.hpp"
/*
 * Grids and PDEs.
 */
#include "DiscreteContinuumMesh.hpp"
#include "VtkMeshWriter.hpp"
#include "SimpleLinearEllipticFiniteDifferenceSolver.hpp"
#include "DiscreteSource.hpp"
#include "VesselBasedDiscreteSource.hpp"
#include "DiscreteContinuumBoundaryCondition.hpp"
#include "DiscreteContinuumLinearEllipticPde.hpp"
#include "MichaelisMentenSteadyStateDiffusionReactionPde.hpp"
/*
 * This should appear last.
 */
#include "PetscSetupAndFinalize.hpp"
class TestSolvingPdesLiteratePaper : public AbstractCellBasedWithTimingsTestSuite
{
public:
    /*
     * ## Test 1 - Linear Reaction Diffusion With Finite Differences =
     * In the first example we will solve a steady-state linear reaction diffusion
     * PDE with finite differences.
     */
    void TestLinearReactionDiffusionPdeWithFiniteDifferences()
    {
        auto p_handler =
                std::make_shared<OutputFileHandler>("TestSolvingPdesLiteratePaper/TestLinearReactionDiffusionPdeWithFiniteDifferences");
        /*
         * We will work in microns
         */
        QLength reference_length(1_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);
        /*
         * Set up a simulation domain, which will be a cuboid.
         */
        QLength domain_width(100_um);
        QLength domain_height(100_um);
        QLength domain_depth(20_um);
        auto p_domain = Part<3>::Create();
        p_domain->AddCuboid(domain_width, domain_height, domain_depth);
        /*
         * Make a regular grid on the domain
         */
        auto p_grid = RegularGrid<3>::Create();
        p_grid->GenerateFromPart(p_domain, 10.0*reference_length);
        /*
         * Set up a PDE, we will model oxygen diffusion.
         */
        auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<3>::Create();
        QDiffusivity oxygen_diffusivity(1.e-6*unit::metre_squared_per_second);
        p_oxygen_pde->SetIsotropicDiffusionConstant(oxygen_diffusivity);
        /*
         * Add continuum sink term for cells
         */
        QRate oxygen_consumption_rate(1.e-6*unit::per_second);
        p_oxygen_pde->SetContinuumLinearInUTerm(-oxygen_consumption_rate);
        /*
         * Add a Dirichlet boundary condition on the left face of the domain.
         */
        p_domain->AddAttributeToPolygonIfFound(Vertex<3>(0.0_um, domain_height/2.0, domain_depth/2.0), "boundary_1", 1.0);

        auto p_left_face_boundary = DiscreteContinuumBoundaryCondition<3>::Create();
        p_left_face_boundary->SetType(BoundaryConditionType::POLYGON);
        p_left_face_boundary->SetDomain(p_domain);
        p_left_face_boundary->SetValue(10e-3_M);
        p_left_face_boundary->SetLabel("boundary_1");
        /*
         * Set up the PDE solvers for the oxygen problem
         */
        auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<3>::Create();
        p_oxygen_solver->SetPde(p_oxygen_pde);
        p_oxygen_solver->SetGrid(p_grid);
        p_oxygen_solver->AddBoundaryCondition(p_left_face_boundary);
        p_oxygen_solver->SetLabel("oxygen");
        p_oxygen_solver->SetFileHandler(p_handler);
        p_oxygen_solver->SetFileName("fd_solution.vti");
        p_oxygen_solver->SetWriteSolution(true);
        p_oxygen_solver->Setup();
        p_oxygen_solver->Solve();
    }
};

#endif /*TESTSOLVINGPDESLITERATEPAPER_HPP_*/
