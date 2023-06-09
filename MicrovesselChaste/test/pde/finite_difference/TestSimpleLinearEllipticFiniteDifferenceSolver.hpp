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

#ifndef TESTSIMPLELINEARELLIPTICFINITEDIFFERENCESOLVER_HPP_
#define TESTSIMPLELINEARELLIPTICFINITEDIFFERENCESOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#include <string>
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include "DiscreteContinuumLinearEllipticPde.hpp"
#include "SmartPointers.hpp"
#include "Part.hpp"
#include "SimpleLinearEllipticFiniteDifferenceSolver.hpp"
#include "VesselNetwork.hpp"
#include "VesselNetworkGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "RegularGrid.hpp"
#include "BaseUnits.hpp"

#include "PetscAndVtkSetupAndFinalize.hpp"

class TestSimpleLinearEllipticFiniteDifferenceSolver : public CxxTest::TestSuite
{

public:

    void TestRectangleDomain()
    {
        // Set up the grid
        BaseUnits::Instance()->SetReferenceLengthScale(1.0*unit::metres);
        BaseUnits::Instance()->SetReferenceConcentrationScale(1.0*unit::mole_per_metre_cubed);
        BaseUnits::Instance()->SetReferenceTimeScale(1.0*unit::seconds);

        std::shared_ptr<Part<2> > p_domain = Part<2>::Create();
        p_domain->AddRectangle(5.0*unit::metres,
                               5.0*unit::metres,
                               Vertex<2>(0.0, 0.0, 0.0));
        std::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
        p_grid->GenerateFromPart(p_domain, 1.0*unit::metres);

        // Choose the PDE
        std::shared_ptr<DiscreteContinuumLinearEllipticPde<2> > p_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
        QDiffusivity diffusivity(1.0* unit::metre_squared_per_second);
        QConcentrationFlowRate consumption_rate(-0.05 * unit::mole_per_metre_cubed_per_second);
        p_pde->SetIsotropicDiffusionConstant(diffusivity);
        p_pde->SetContinuumConstantInUTerm(consumption_rate);

        // Prescribe a value on the domain's left boundary
        std::shared_ptr<DiscreteContinuumBoundaryCondition<2> > p_boundary_condition = DiscreteContinuumBoundaryCondition<2>::Create();
        QConcentration boundary_concentration(1.0* unit::mole_per_metre_cubed);
        p_boundary_condition->SetValue(boundary_concentration);
        p_boundary_condition->SetType(BoundaryConditionType::POINT);
        vtkSmartPointer<vtkPoints> p_boundary_points = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkPoints> p_points = p_grid->GetPoints();
        for(unsigned idx=0; idx<p_points->GetNumberOfPoints(); idx++)
        {
            if(p_points->GetPoint(idx)[0]==0.0)
            {
                p_boundary_points->InsertNextPoint(p_points->GetPoint(idx));
            }
        }
        p_boundary_condition->SetPoints(p_boundary_points);

        // Set up and run the solver
        SimpleLinearEllipticFiniteDifferenceSolver<2> solver;
        solver.SetGrid(p_grid);
        solver.SetPde(p_pde);
        solver.AddBoundaryCondition(p_boundary_condition);

        auto p_output_file_handler =
        		std::make_shared<OutputFileHandler>("TestSimpleLinearEllipticFiniteDifferenceSolver/RectangleDomain");

        solver.SetFileHandler(p_output_file_handler);
        solver.SetWriteSolution(true);
        solver.Solve();

        std::vector<QConcentration > solution = solver.GetConcentrations();

        // Analytical c = k*x*x/(2*D) - k*x*w/D+c_0
        for(unsigned idx=0; idx<6; idx++)
        {
            QLength x = double(idx)*1.0*unit::metres;
            QLength w = 5.0*unit::metres;
            QConcentration c = -consumption_rate*x*x/(2.0*diffusivity)-
                    x*-consumption_rate*w/diffusivity + boundary_concentration;
            double norm_analytical = c/(1.0* unit::mole_per_metre_cubed);
            double norm_numerical = solution[idx]/(1.0* unit::mole_per_metre_cubed);
            TS_ASSERT_DELTA(norm_analytical, norm_numerical, 1.e-6)
        }
    }

    void TestCuboidalDomain()
    {
        // Set up the grid
        std::shared_ptr<Part<3> > p_domain = Part<3>::Create();
        p_domain->AddCuboid(200.0*1_um,
                            10.0*1_um,
                            10.0*1_um,
                            Vertex<3>(0.0, 0.0, 0.0));
        std::shared_ptr<RegularGrid<3> > p_grid = RegularGrid<3>::Create();
        p_grid->GenerateFromPart(p_domain, 10.0*1_um);

        // Choose the PDE
        std::shared_ptr<DiscreteContinuumLinearEllipticPde<3> > p_pde = DiscreteContinuumLinearEllipticPde<3>::Create();
        QDiffusivity diffusivity(1.e-3 * unit::metre_squared_per_second);
        QConcentrationFlowRate consumption_rate(-20.0 * unit::mole_per_metre_cubed_per_second);
        p_pde->SetIsotropicDiffusionConstant(diffusivity);
        p_pde->SetContinuumConstantInUTerm(consumption_rate);

        // Prescribe a value on the domain boundaries
        std::shared_ptr<DiscreteContinuumBoundaryCondition<3> > p_boundary_condition = DiscreteContinuumBoundaryCondition<3>::Create();
        QConcentration boundary_concentration(1.e-3 * unit::mole_per_metre_cubed);
        p_boundary_condition->SetValue(boundary_concentration);
        p_boundary_condition->SetType(BoundaryConditionType::POINT);
        vtkSmartPointer<vtkPoints> p_boundary_points = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkPoints> p_points = p_grid->GetPoints();
        for(unsigned idx=0; idx<p_points->GetNumberOfPoints(); idx++)
        {
            if(p_points->GetPoint(idx)[0]==0.0)
            {
                p_boundary_points->InsertNextPoint(p_points->GetPoint(idx));
            }
        }
        p_boundary_condition->SetPoints(p_boundary_points);
        // Set up and run the solver

        SimpleLinearEllipticFiniteDifferenceSolver<3> solver;
        solver.SetGrid(p_grid);
        solver.SetPde(p_pde);
        solver.AddBoundaryCondition(p_boundary_condition);

        auto p_output_file_handler =
        		std::make_shared<OutputFileHandler>("TestSimpleLinearEllipticFiniteDifferenceSolver/CuboidalDomain");
        solver.SetFileHandler(p_output_file_handler);
        solver.SetWriteSolution(true);
        solver.Solve();

        std::vector<QConcentration > solution = solver.GetConcentrations();

        // Analytical c = k*x*x/(2*D) - k*x*w/D+c_0
        for(unsigned idx=0; idx<6; idx++)
        {
            QLength x = double(idx)*10.0e-6*unit::metres;
            QLength w = 200.0e-6*unit::metres;
            QConcentration c = -consumption_rate*x*x/(2.0*diffusivity)-
                    x*-consumption_rate*w/diffusivity + boundary_concentration;
            double norm_analytical = c/(1.0* unit::mole_per_metre_cubed);
            double norm_numerical = solution[idx]/(1.0* unit::mole_per_metre_cubed);
            TS_ASSERT_DELTA(norm_analytical, norm_numerical, 1.e-6)
        }
    }

    void TestWithVesselBoundaryConditions()
    {
        // Set up the vessel network
        QLength vessel_length = 100.0 * 1_um;
        VesselNetworkGenerator<3> generator;
        std::shared_ptr<VesselNetwork<3> > p_network = generator.GenerateSingleVessel(vessel_length,
                                                                                        Vertex<3>(0.0, 0.0, 0.0));

        // Set up the grid
        std::shared_ptr<Part<3> > p_domain = Part<3>::Create();
        p_domain->AddCuboid(vessel_length, vessel_length, vessel_length, Vertex<3>(0.0, 0.0, 0.0));
        std::shared_ptr<RegularGrid<3> > p_grid = RegularGrid<3>::Create();
        p_grid->GenerateFromPart(p_domain, 10.0*1_um);

        // Choose the PDE
        std::shared_ptr<DiscreteContinuumLinearEllipticPde<3> > p_pde = DiscreteContinuumLinearEllipticPde<3>::Create();
        QDiffusivity diffusivity(0.0033 * unit::metre_squared_per_second);
        QRate consumption_rate(-2.e6 * unit::per_second);
        p_pde->SetIsotropicDiffusionConstant(diffusivity);
        p_pde->SetContinuumLinearInUTerm(consumption_rate);

        // Set up the boundary condition
        std::shared_ptr<DiscreteContinuumBoundaryCondition<3> > p_vessel_boundary_condition = DiscreteContinuumBoundaryCondition<3>::Create();
        QConcentration boundary_concentration(1.e-3 * unit::mole_per_metre_cubed);
        p_vessel_boundary_condition->SetValue(boundary_concentration);
        p_vessel_boundary_condition->SetType(BoundaryConditionType::VESSEL_LINE);
        p_vessel_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);

        // Set up and run the solver
        SimpleLinearEllipticFiniteDifferenceSolver<3> solver;
        solver.SetGrid(p_grid);
        solver.SetPde(p_pde);
        solver.AddBoundaryCondition(p_vessel_boundary_condition);
        solver.SetVesselNetwork(p_network);

        auto p_output_file_handler =
        		std::make_shared<OutputFileHandler>("TestSimpleLinearEllipticFiniteDifferenceSolver/WithVessels", true);
        solver.SetFileHandler(p_output_file_handler);
        solver.SetWriteSolution(true);
        solver.Solve();
    }
};

#endif /*TESTSIMPLELINEARELLIPTICFINITEDIFFERENCESOLVER_HPP_*/
