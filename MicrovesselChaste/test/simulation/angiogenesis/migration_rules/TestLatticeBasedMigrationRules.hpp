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

#ifndef TESTLATTICEBASEDMIGRATIONRULES_HPP
#define TESTLATTICEBASEDMIGRATIONRULES_HPP

#include <cxxtest/TestSuite.h>
#include "VesselImpedanceCalculator.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "VesselNetworkGenerator.hpp"
#include "VesselNode.hpp"
#include "VesselSegment.hpp"
#include "Vessel.hpp"
#include "VesselNetwork.hpp"
#include "Part.hpp"
#include "AngiogenesisSolver.hpp"
#include "VesselSegment.hpp"
#include "LatticeBasedMigrationRule.hpp"
#include "FunctionMap.hpp"
#include "Owen2011MigrationRule.hpp"
#include "FlowSolver.hpp"
#include "UnitCollection.hpp"
#include "GridCalculator.hpp"
#include "VesselNetworkPropertyManager.hpp"

#include "PetscAndVtkSetupAndFinalize.hpp"

class TestLatticeBasedMigrationRules : public AbstractCellBasedWithTimingsTestSuite
{

public:

    void TestSingleVessel()
    {
        // Set the grid to move on
        auto p_grid = RegularGrid<2>::Create();
        QLength spacing = 100_um;
        p_grid->SetSpacing(spacing);
        c_vector<double, 3> dimensions;
        dimensions[0] = 7; // num x
        dimensions[1] = 5; // num_y
        dimensions[2] = 1; // num_z
        p_grid->SetDimensions(dimensions);

        auto p_grid_calc = GridCalculator<2>::Create();
        p_grid_calc->SetGrid(p_grid);

        // Make a vessel
        auto p_node1 = VesselNode<2>::Create(0.0_m, 2.0*spacing);
        auto p_node2 = VesselNode<2>::Create(spacing, 2.0*spacing);
        p_node2->SetIsMigrating(true);
        auto p_vessel = Vessel<2>::Create(p_node1, p_node2);
        auto p_network = VesselNetwork<2>::Create();
        p_network->AddVessel(p_vessel);
        p_grid_calc->SetVesselNetwork(p_network);

        // Set up the migration rule
        auto p_migration_rule = LatticeBasedMigrationRule<2>::Create();
        p_migration_rule->SetGridCalculator(p_grid_calc);
        p_migration_rule->SetMovementProbability(0.1);
        p_migration_rule->SetNetwork(p_network);

        // Test that we move into the correct locations and that sometimes, but not always, we don't move
        RandomNumberGenerator::Instance()->Reseed(522525);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        unsigned not_moved = 0;
        for(unsigned idx=0; idx<100; idx++)
        {
            std::vector<int> indices = p_migration_rule->GetIndices(std::vector<VesselNodePtr<2> > (1, p_node2));
            if (indices[0] == -1)
            {
                not_moved++;
            }
            else
            {
                TS_ASSERT(indices[0] == 16 or indices[0] == 22 or indices[0] == 8)
            }
        }
        TS_ASSERT(not_moved>0)
        TS_ASSERT(not_moved<100)
    }

    void TestOwen2011SingleVessel()
    {
        // Set the grid to move on
        auto p_grid = RegularGrid<2>::Create();
        QLength spacing = 100_um; //um
        p_grid->SetSpacing(spacing);
        c_vector<double, 3> dimensions;
        dimensions[0] = 7; // num x
        dimensions[1] = 5; // num_y
        dimensions[2] = 1; // num_z
        p_grid->SetDimensions(dimensions);

        auto p_grid_calc = GridCalculator<2>::Create();
        p_grid_calc->SetGrid(p_grid);

        // Make a vessel
        auto p_node1 = VesselNode<2>::Create(0.0_m, 2.0*spacing);
        auto p_node2 = VesselNode<2>::Create(spacing, 2.0*spacing);
        p_node2->SetIsMigrating(true);
        auto p_vessel = Vessel<2>::Create(p_node1, p_node2);
        auto p_network = VesselNetwork<2>::Create();
        p_network->AddVessel(p_vessel);
        p_grid_calc->SetVesselNetwork(p_network);

        // Set up a vegf field
        std::shared_ptr<FunctionMap<2> > p_funciton_map = FunctionMap<2>::Create();
        p_funciton_map->SetGrid(p_grid);
        std::vector<QConcentration> vegf_field =
                std::vector<QConcentration >(dimensions[0]*dimensions[1], 0_M);

        QConcentration max_vegf(0.2e-6_nM);
        for(unsigned idx=0; idx<p_grid_calc->GetGrid()->GetNumberOfPoints(); idx++)
        {
            vegf_field[idx] = double(p_grid->GetPoint(idx).Convert(1_um)[0] / (double(dimensions[0]) * spacing) )*max_vegf;
        }
        p_funciton_map->UpdateSolution(vegf_field);

        // Set up the migration rule
        auto p_migration_rule = Owen2011MigrationRule<2>::Create();
        p_migration_rule->SetGridCalculator(p_grid_calc);
        p_migration_rule->SetMovementProbability(0.001);
        p_migration_rule->SetNetwork(p_network);
        p_migration_rule->SetDiscreteContinuumSolver(p_funciton_map);

        // Test that we move into the correct locations and that sometimes, but not always, we don't move.
        // Also check that we mostly move in the direction of the vegf gradient, but not always
        RandomNumberGenerator::Instance()->Reseed(522525);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.e-6, 1);

        unsigned not_moved = 0;
        unsigned num_right = 0;
        for(unsigned idx=0; idx<100; idx++)
        {
            std::vector<int> indices = p_migration_rule->GetIndices(std::vector<VesselNodePtr<2> > (1, p_node2));
            if (indices[0] == -1)
            {
                not_moved++;
            }
            else
            {
                TS_ASSERT(indices[0] == 16 or indices[0] == 22 or indices[0] == 8)
                if(indices[0] == 16)
                {
                    num_right ++;
                }
            }
        }
        TS_ASSERT(not_moved>0)
        TS_ASSERT(not_moved<100)
        TS_ASSERT(num_right>0)
        TS_ASSERT(num_right<100)
    }

    void TestSproutingAndMigrationWithFlow()
    {
        // Make a network

        // Set the grid to move on
        auto p_grid = RegularGrid<2>::Create();
        QLength spacing = 100.0_um; //um
        p_grid->SetSpacing(spacing);
        c_vector<double, 3> dimensions;
        dimensions[0] = 7; // num x
        dimensions[1] = 5; // num_y
        dimensions[2] = 1; // num_z
        p_grid->SetDimensions(dimensions);

        std::shared_ptr<GridCalculator<2> > p_grid_calc = GridCalculator<2>::Create();
        p_grid_calc->SetGrid(p_grid);

        std::shared_ptr<VesselNode<2> > p_node1 = VesselNode<2>::Create(0.0_m, 2.0*spacing);
        std::shared_ptr<VesselNode<2> > p_node2 = VesselNode<2>::Create(spacing, 2.0*spacing);
        std::shared_ptr<VesselNode<2> > p_node3 = VesselNode<2>::Create(2.0*spacing, 2.0*spacing);
        std::shared_ptr<VesselNode<2> > p_node4 = VesselNode<2>::Create(3.0*spacing, 2.0*spacing);
        std::shared_ptr<VesselNode<2> > p_node5 = VesselNode<2>::Create(4.0*spacing, 2.0*spacing);

        p_node3->SetIsMigrating(true);
        auto p_vessel1 = Vessel<2>::Create(p_node1, p_node2);
        auto p_vessel2 = Vessel<2>::Create(p_node2, p_node3);
        auto p_vessel3 = Vessel<2>::Create(p_node3, p_node4);
        auto p_vessel4 = Vessel<2>::Create(p_node4, p_node5);
        auto p_network = VesselNetwork<2>::Create();
        p_network->AddVessel(p_vessel1);
        p_network->AddVessel(p_vessel2);
        p_network->AddVessel(p_vessel3);
        p_network->AddVessel(p_vessel4);
        p_grid_calc->SetVesselNetwork(p_network);

        p_node1->GetFlowProperties()->SetIsInputNode(true);
        p_node1->GetFlowProperties()->SetPressure(3000.0_Pa);
        p_node5->GetFlowProperties()->SetIsOutputNode(true);
        p_node5->GetFlowProperties()->SetPressure(1000.0_Pa);

        VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, 10.0_um);
        std::vector<std::shared_ptr<VesselSegment<2> > > segments = p_network->GetVesselSegments();
        for(unsigned idx=0; idx<segments.size(); idx++)
        {
            segments[idx]->GetFlowProperties()->SetViscosity(1.e-3 * unit::poiseuille);
        }

        auto p_migration_rule = LatticeBasedMigrationRule<2>::Create();
        p_migration_rule->SetGridCalculator(p_grid_calc);
        p_migration_rule->SetMovementProbability(0.5);
        p_migration_rule->SetNetwork(p_network);

        VesselImpedanceCalculator<2> calculator;
        calculator.SetVesselNetwork(p_network);
        calculator.Calculate();

        FlowSolver<2> flow_solver;
        flow_solver.SetVesselNetwork(p_network);
        flow_solver.SetUp();
        flow_solver.Solve();

        // Grow the vessel
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(30, 30);
        AngiogenesisSolver<2> angiogenesis_solver;
        angiogenesis_solver.SetVesselNetwork(p_network);
        angiogenesis_solver.SetMigrationRule(p_migration_rule);
        angiogenesis_solver.SetVesselGridCalculator(p_grid_calc);

        auto p_handler = std::make_shared<OutputFileHandler>("TestLatticeBasedMigrationRulesWithFlow");
        angiogenesis_solver.SetOutputFileHandler(p_handler);
        angiogenesis_solver.Run(true);

        flow_solver.SetUp();
        flow_solver.Solve();
        p_network->Write(p_handler->GetOutputDirectoryFullPath() + "/network_flow_after.vtp");
    }
};

#endif // TESTLATTICEBASEDMIGRATIONRULES_HPP
