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

#ifndef TESTDISCRETESOURCE_HPP_
#define TESTDISCRETESOURCE_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#include <string>
#include "DiscreteContinuumLinearEllipticPde.hpp"
#include "SmartPointers.hpp"
#include "Part.hpp"
#include "MichaelisMentenSteadyStateDiffusionReactionPde.hpp"
#include "SimpleLinearEllipticFiniteDifferenceSolver.hpp"
#include "SimpleNonLinearEllipticFiniteDifferenceSolver.hpp"
#include "VesselNetwork.hpp"
#include "VesselNetworkGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "RegularGrid.hpp"
#include "DiscreteSource.hpp"
#include "Vertex.hpp"
#include "DiscreteContinuumMeshGenerator.hpp"
#include "DiscreteContinuumMesh.hpp"
#include "FunctionMap.hpp"
#include "SimpleLinearEllipticFiniteElementSolver.hpp"
#include "VtkMeshWriter.hpp"
#include "GridCalculator.hpp"

#include "PetscAndVtkSetupAndFinalize.hpp"

class TestDiscreteSource : public CxxTest::TestSuite
{

public:

    void TestGridFunction()
    {
        QLength length(100.0_um);

        // Set up the grid
        std::shared_ptr<Part<2> > p_domain = Part<2>::Create();
        p_domain->AddRectangle(length, length);

        std::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
        QLength grid_spacing(5.0_um);
        p_grid->GenerateFromPart(p_domain, grid_spacing);

        // Set up a density map
        std::shared_ptr<DensityMap<2> > p_density_map = DensityMap<2>::Create();
        p_density_map->SetGrid(p_grid);

        // Set up the discrete source
        std::vector<Vertex<2> > linear_consumption_points;
        linear_consumption_points.push_back(Vertex<2>(50.0_um, 50.0_um));
        std::shared_ptr<DiscreteSource<2> > p_linear_point_source = DiscreteSource<2>::Create();

        p_linear_point_source->SetLinearInUValue(1.0 * unit::per_second);
        p_linear_point_source->SetPoints(linear_consumption_points);
        p_linear_point_source->SetDensityMap(p_density_map);

        std::shared_ptr<DiscreteSource<2> > p_const_point_source = DiscreteSource<2>::Create();
        QConcentrationFlowRate consumption_rate(2.0 * unit::mole_per_metre_cubed_per_second);
        p_const_point_source->SetConstantInUValue(consumption_rate);
        std::vector<Vertex<2> > constant_consumption_points;
        constant_consumption_points.push_back(Vertex<2>(25.0_um, 25.0_um, 0.0));
        constant_consumption_points.push_back(Vertex<2>(75.0_um, 25.0_um, 0.0));
        constant_consumption_points.push_back(Vertex<2>(75.0_um, 75.0_um, 0.0));
        constant_consumption_points.push_back(Vertex<2>(25.0_um, 75.0_um, 0.0));
        p_const_point_source->SetPoints(constant_consumption_points);
        p_const_point_source->SetDensityMap(p_density_map);

        // Set up a function map
        FunctionMap<2> solver;
        solver.SetGrid(p_grid);

        // Get the source values at each point on the grid
        std::vector<QRate > point_rates = p_linear_point_source->GetLinearInUValues();
        std::vector<QConcentrationFlowRate > point_conc_rates = p_const_point_source->GetConstantInUValues();
        std::vector<double> solution;
        for(unsigned idx=0; idx<p_density_map->GetGridCalculator()->GetGrid()->GetNumberOfPoints(); idx++)
        {
            solution.push_back(double(point_rates[idx].GetValue() + point_conc_rates[idx].GetValue()));
        }

        solver.UpdateSolution(solution);
        auto p_output_file_handler =
                std::make_shared<OutputFileHandler>("TestDiscreteSource/TestGridFunction");

        solver.SetFileHandler(p_output_file_handler);
        solver.Write();
    }

    void TestMeshFunction()
    {
        QLength length(100.0*unit::microns);

        // Set up the grid
        std::shared_ptr<Part<2> > p_domain = Part<2>::Create();
        p_domain->AddRectangle(length, length);

        auto p_mesh_generator = DiscreteContinuumMeshGenerator<2>::Create();
        p_mesh_generator->SetDomain(p_domain);
        p_mesh_generator->SetMaxElementArea(Qpow3(0.02*length));
        p_mesh_generator->Update();

        std::shared_ptr<DensityMap<2> > p_density_map = DensityMap<2>::Create();
        p_density_map->SetGrid(p_mesh_generator->GetMesh());

        // Set up the discrete source
        std::vector<Vertex<2> > linear_consumption_points;
        linear_consumption_points.push_back(Vertex<2>(50.0_um, 50.0_um));
        auto p_linear_point_source = DiscreteSource<2>::Create();
        p_linear_point_source->SetLinearInUValue(1.0 * unit::per_second);
        p_linear_point_source->SetPoints(linear_consumption_points);
        p_linear_point_source->SetDensityMap(p_density_map);

        auto p_const_point_source = DiscreteSource<2>::Create();
        QConcentrationFlowRate consumption_rate(2.0 * unit::mole_per_metre_cubed_per_second);
        p_const_point_source->SetConstantInUValue(consumption_rate);
        std::vector<Vertex<2> > constant_consumption_points;
        constant_consumption_points.push_back(Vertex<2>(25.0_um, 25.0_um, 25.0_um));
        constant_consumption_points.push_back(Vertex<2>(75.0_um, 25.0_um, 25.0_um));
        constant_consumption_points.push_back(Vertex<2>(75.0_um, 75.0_um, 25.0_um));
        constant_consumption_points.push_back(Vertex<2>(25.0_um, 75.0_um, 25.0_um));
        p_const_point_source->SetPoints(constant_consumption_points);
        p_const_point_source->SetDensityMap(p_density_map);

        // Set up a function map
        FunctionMap<2> solver;
        solver.SetGrid(p_mesh_generator->GetMesh());

        // Get the source values at each point on the grid
        std::vector<QRate > point_rates = p_linear_point_source->GetLinearInUValues();
        std::vector<QConcentrationFlowRate > point_conc_rates = p_const_point_source->GetConstantInUValues();
        std::vector<double> solution;
        for(unsigned idx=0; idx<point_conc_rates.size(); idx++)
        {
            solution.push_back(double(point_rates[idx].GetValue() + point_conc_rates[idx].GetValue()));
        }
        solver.UpdateElementSolution(solution);
        auto p_output_file_handler =
                std::make_shared<OutputFileHandler>("TestDiscreteSource/TestMeshFunction");
        solver.SetFileHandler(p_output_file_handler);
        solver.Write();
    }

    void TestLinearGridPde()
    {
        BaseUnits::Instance()->SetReferenceLengthScale(1.0*unit::metres);
        BaseUnits::Instance()->SetReferenceConcentrationScale(1.0*unit::mole_per_metre_cubed);

        QLength length(100.0*unit::metres);

        // Set up the grid
        std::shared_ptr<Part<2> > p_domain = Part<2>::Create();
        p_domain->AddRectangle(length, length, Vertex<2>());
        std::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
        QLength grid_spacing(5.0*unit::metres);
        p_grid->GenerateFromPart(p_domain, grid_spacing);

        // Choose the PDE
        std::shared_ptr<DiscreteContinuumLinearEllipticPde<2> > p_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
        QDiffusivity diffusivity(0.0033 * unit::metre_squared_per_second);
        p_pde->SetIsotropicDiffusionConstant(diffusivity);

        // Set up the discrete source
        std::vector<Vertex<2> > linear_consumption_points;
        linear_consumption_points.push_back(Vertex<2>(50.0_um, 50.0_um));
        std::shared_ptr<DiscreteSource<2> > p_linear_point_source = DiscreteSource<2>::Create();
        p_linear_point_source->SetLinearInUValue(-1.0 * unit::per_second);
        p_linear_point_source->SetPoints(linear_consumption_points);

        std::shared_ptr<DiscreteSource<2> > p_const_point_source = DiscreteSource<2>::Create();
        QConcentrationFlowRate consumption_rate(2.0 * unit::mole_per_metre_cubed_per_second);
        p_const_point_source->SetConstantInUValue(consumption_rate);
        std::vector<Vertex<2> > constant_consumption_points;
        constant_consumption_points.push_back(Vertex<2>(25.0_um, 25.0_um));
        constant_consumption_points.push_back(Vertex<2>(75.0_um, 25.0_um));
        constant_consumption_points.push_back(Vertex<2>(75.0_um, 75.0_um));
        constant_consumption_points.push_back(Vertex<2>(25.0_um, 75.0_um));
        p_const_point_source->SetPoints(constant_consumption_points);

        p_pde->AddDiscreteSource(p_const_point_source);
        p_pde->AddDiscreteSource(p_linear_point_source);

        auto p_boundary2 = DiscreteContinuumBoundaryCondition<2>::Create();
        p_boundary2->SetValue(3.0*unit::mole_per_metre_cubed);

        // Set up and run the solver
        SimpleLinearEllipticFiniteDifferenceSolver<2> solver;
        solver.SetGrid(p_grid);
        solver.SetPde(p_pde);
        solver.AddBoundaryCondition(p_boundary2);
        auto p_output_file_handler =
                std::make_shared<OutputFileHandler>("TestDiscreteSource/TestLinearGridPde");
        solver.SetFileHandler(p_output_file_handler);
        solver.SetWriteSolution(true);
        solver.Solve();
    }

    void TestNonLinearGridPde()
    {
        QLength length(100.0*unit::microns);

        // Set up the grid
        std::shared_ptr<Part<2> > p_domain = Part<2>::Create();
        p_domain->AddRectangle(length, length);
        std::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
        QLength grid_spacing(5.0*unit::microns);
        p_grid->GenerateFromPart(p_domain, grid_spacing);

        // Choose the PDE
        auto p_pde = MichaelisMentenSteadyStateDiffusionReactionPde<2>::Create();
        QDiffusivity diffusivity(0.0033 * unit::metre_squared_per_second);

        p_pde->SetIsotropicDiffusionConstant(diffusivity);
        p_pde->SetMichaelisMentenThreshold(2.0 * unit::mole_per_metre_cubed);

        // Set up the discrete source
        std::vector<Vertex<2> > linear_consumption_points;
        linear_consumption_points.push_back(Vertex<2>(50.0_um, 50.0_um, 50.0_um));
        auto p_linear_point_source = DiscreteSource<2>::Create();
        p_linear_point_source->SetLinearInUValue(-1.0 * unit::per_second);
        p_linear_point_source->SetPoints(linear_consumption_points);

        auto p_const_point_source = DiscreteSource<2>::Create();
        QConcentrationFlowRate consumption_rate(2.e-4 * unit::mole_per_metre_cubed_per_second);
        p_const_point_source->SetConstantInUValue(consumption_rate);
        std::vector<Vertex<2> > constant_consumption_points;
        constant_consumption_points.push_back(Vertex<2>(25.0_um, 25.0_um, 25.0_um));
        constant_consumption_points.push_back(Vertex<2>(75.0_um, 25.0_um, 25.0_um));
        constant_consumption_points.push_back(Vertex<2>(75.0_um, 75.0_um, 25.0_um));
        constant_consumption_points.push_back(Vertex<2>(25.0_um, 75.0_um, 25.0_um));
        p_const_point_source->SetPoints(constant_consumption_points);

        auto p_boundary2 = DiscreteContinuumBoundaryCondition<2>::Create();
        p_boundary2->SetValue(3.e-6*unit::mole_per_metre_cubed);
        p_pde->AddDiscreteSource(p_linear_point_source);
        p_pde->AddDiscreteSource(p_const_point_source);

        // Set up and run the solver
        SimpleNonLinearEllipticFiniteDifferenceSolver<2> solver;
        solver.SetGrid(p_grid);
        solver.SetPde(p_pde);
        solver.AddBoundaryCondition(p_boundary2);
        auto p_output_file_handler =
                std::make_shared<OutputFileHandler>("TestDiscreteSource/TestNonlinearGridPde");
        solver.SetFileHandler(p_output_file_handler);
        solver.SetWriteSolution(true);
        solver.Solve();
    }
};

#endif /*TESTDISCRETESOURCE_HPP_*/
