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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Vedang's Notes (vedang.narain@msdtc.ox.ac.uk)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*

This script contains tests that are used to simulate blood flow and oxygenation in different vascular architectures. 

Outputs can usually be found in /tmp/<home_directory>/testoutput. The pointwise data for the field can be obtained by opening spreadsheet view in Paraview.

H = Haematocrit
BC = boundary condition
RT = Radiotherapy
PQ = Perfusion Quotient
Fig. = Figure
# = number

TO DO LIST
 
Currently, we manually set no-flow vessels to have zero H, but maybe we should either put that bit of code in the solvers. Or set all relevant parameters (like viscosity) to zero too.

In the homogeneous cases, the Pries & Memory solvers are breaking the Voronoi simulation. The solvers are currently replaced by Const. H. The Fung solver breaks it with some layouts.

The Pries w. Memory solvers have been ignored for hex. and Voronoi since the network doesn't have directionality.

In the heterogeneous case, the Yang, Pries Memory, and Pries solvers are breaking the dichotomous simulation when alpha <6. When alpha<5, all but Yang work (which breaks at alpha=2). The Fung and Pries solvers break the Voronoi sometimes. The Fung, Pries, and Yang (sometimes) solvers don't converge in the hex network. The broken solvers are currently replaced by Const. H

The output node extends ~0.4 um beyond the x-extent of the simulation. It's because of the grid spacing resolution.

I've marked unknown functionality with a '???'. Need to read up on it later.

For some reason, you need to have an extra node for an unpruned network to have flow. However, single-path networks don't flow if you leave that node in. WHY? I've added a node at a distance of 100 picometres from the last one to circumvent the problem for now.

I should move my equilateral network build to the Network Generator function, once I figure out why that extra node is needed.

Link the inlet haematocrit (0.45) directly with Owen11 parameters.

The y-extent for the hexagonal network is currently set manually based on the simulations. Find a way to make it appropriately cover the hexagonal network automatically.

13/6/23

*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialisation
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef TESTHAEMATOCRITSOLVERS_HPP_
#define TESTHAEMATOCRITSOLVERS_HPP_
#define _BACKWARD_BACKWARD_WARNING_H 1  // Cut out the VTK deprecated warning

// Essential functionality
#include <boost/lexical_cast.hpp>
#include <cxxtest/TestSuite.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>  // needed for exp function
#include <sstream>
#include <string>
#include <time.h>
#include <vector>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"
#include "SimulationTime.hpp"
#include "SmartPointers.hpp"

// Geometry tools
#include "Part.hpp"

// Dimensional analysis
#include "BaseUnits.hpp"
#include "UnitCollection.hpp"
#include "GenericParameters.hpp"
#include "ParameterCollection.hpp"
#include "Owen11Parameters.hpp"
#include "Secomb04Parameters.hpp"
#include "Vertex.hpp"

// Grids and PDEs
#include "CellBasedDiscreteSource.hpp"
#include "CellStateDependentDiscreteSource.hpp"
#include "SimpleLinearEllipticFiniteDifferenceSolver.hpp"
#include "SimpleLinearEllipticFiniteElementSolver.hpp"
#include "DiscreteContinuumBoundaryCondition.hpp"
#include "DiscreteContinuumLinearEllipticPde.hpp"
#include "RegularGrid.hpp"
#include "GridCalculator.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "EllipticGrowingDomainPdeModifier.hpp"
#include "VesselBasedDiscreteSource.hpp"

// Vessel networks
#include "Vessel.hpp"
#include "VesselNetwork.hpp"
#include "VesselNetworkGenerator.hpp"
#include "VesselNetworkGeometryCalculator.hpp"
#include "VesselNetworkPropertyManager.hpp"
#include "VesselNode.hpp"
#include "VesselSegment.hpp"

// Flow
#include "NodeFlowProperties.hpp"
#include "VesselImpedanceCalculator.hpp"
#include "FlowSolver.hpp"
#include "WallShearStressCalculator.hpp"
#include "MechanicalStimulusCalculator.hpp"
#include "MetabolicStimulusCalculator.hpp"
#include "ShrinkingStimulusCalculator.hpp"
#include "StructuralAdaptationSolver.hpp"
#include "ViscosityCalculator.hpp"

// Haematocrit
// #include "AlarconHaematocritSolver.hpp"
#include "BetteridgeHaematocritSolver.hpp"
#include "ConstantHaematocritSolver.hpp"
// #include "LinnengarHaematocritSolver.hpp"
#include "YangHaematocritSolver.hpp"
#include "PriesHaematocritSolver.hpp"
#include "PriesWithMemoryHaematocritSolver.hpp"

// Cells
// #include "CancerCellMutationState.hpp"
// #include "StalkCellMutationState.hpp"
// #include "QuiescentCancerCellMutationState.hpp"
// #include "WildTypeCellMutationState.hpp"
// #include "Owen11CellPopulationGenerator.hpp"
// #include "Owen2011TrackingModifier.hpp"
// #include "CaBasedCellPopulation.hpp"
// #include "ApoptoticCellKiller.hpp"
// #include "SimpleOxygenBasedCellCycleModel.hpp"
// #include "StemCellProliferativeType.hpp"

// Forces
#include "GeneralisedLinearSpringForce.hpp"

// Vessel regression solver
#include "WallShearStressBasedRegressionSolver.hpp"

// General solver to collect all the flows
#include "MicrovesselSolver.hpp"
#include "MicrovesselSimulationModifier.hpp"
#include "OnLatticeSimulation.hpp"
#include "OffLatticeSimulation.hpp"

// Visualisation
#include "MicrovesselVtkScene.hpp"
#include "VtkSceneMicrovesselModifier.hpp"

// Keep this last
#include "PetscAndVtkSetupAndFinalize.hpp"

// ???
using namespace std;

// Make a test class
class TestHaematocritSolvers : public CxxTest::TestSuite
{

public:

    // // Start the execution timer
    // clock_t clkStart = clock();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Tests
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define key parameters here, i.e., things that are constant across architectures
    // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));  // = 40 um
    QLength grid_spacing = 10.0_um;  // the simulation time gets quite long if you reduce the resolution further
    // Note: 10.0_um is a rather low resolution for the tiny equilateral network.
    unsigned max_alpha = 3;  // set the max. heterog. parameter (alpha = 1+(max_alpha*0.1), i.e., 5 is 1.5) 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Equilateral Network
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Change resolution of equilateral network for better visualisation
    // QLength equilateral_grid_spacing = grid_spacing;
    QLength equilateral_grid_spacing = 0.1_um;

    // Make an equilateral network on a PDE grid that acts as a Dirichlet BC in 2D 
    void xTestEquilateralNetworkLineSource2D()
    {
        // Run the simulation with different heterogeneities
        for (unsigned n_alpha=0; n_alpha<=max_alpha; n_alpha++)
        {
            double alpha = 1.0+(double)n_alpha*0.1;  // alpha determines the relative radius of the left vessel 
        
            // Set up the reference length for the simulation
            QLength reference_length(1.0_um);
            BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

            // Set key vessel parameters
            QLength vessel_length(100.0*unit::microns);
            QLength vessel_height = (pow(2,0.5)*vessel_length)*0.5;

            // Set up the domain parameters
            QLength domain_x = vessel_height+vessel_height+vessel_length+vessel_length+0.0001_um;  // this should extend to the x-position of outlet node
            QLength domain_y = domain_x;  // should be the same as domain_x to make square domain
            QLength mid_domain_y = domain_y*0.5;

            // Set nodes based on an equilateral network
            std::shared_ptr<VesselNode<2> > p_node_1 = VesselNode<2>::Create(0.0_um, mid_domain_y);
            std::shared_ptr<VesselNode<2> > p_node_2 = VesselNode<2>::Create(vessel_length, mid_domain_y);
            std::shared_ptr<VesselNode<2> > p_node_3 = VesselNode<2>::Create(vessel_height+vessel_length, vessel_height+mid_domain_y);
            std::shared_ptr<VesselNode<2> > p_node_4 = VesselNode<2>::Create(vessel_height+vessel_length, -vessel_height+mid_domain_y);
            std::shared_ptr<VesselNode<2> > p_node_5 = VesselNode<2>::Create(vessel_height+vessel_height+vessel_length, mid_domain_y);
            std::shared_ptr<VesselNode<2> > p_node_6 = VesselNode<2>::Create(vessel_height+vessel_height+vessel_length+vessel_length, mid_domain_y);
            std::shared_ptr<VesselNode<2> > p_node_7 = VesselNode<2>::Create(domain_x, mid_domain_y);  // add a node at a distance of 100 picometres from the last node to work around glitch of disconnected flow when setting p_node_6 as the output node

            // Make segments 
            std::shared_ptr<VesselSegment<2> > p_segment_1 = VesselSegment<2>::Create(p_node_1, p_node_2);
            std::shared_ptr<VesselSegment<2> > p_segment_2 = VesselSegment<2>::Create(p_node_2, p_node_3);
            std::shared_ptr<VesselSegment<2> > p_segment_3 = VesselSegment<2>::Create(p_node_2, p_node_4);
            std::shared_ptr<VesselSegment<2> > p_segment_4 = VesselSegment<2>::Create(p_node_3, p_node_5);
            std::shared_ptr<VesselSegment<2> > p_segment_5 = VesselSegment<2>::Create(p_node_4, p_node_5);
            std::shared_ptr<VesselSegment<2> > p_segment_6 = VesselSegment<2>::Create(p_node_5, p_node_6);
            std::shared_ptr<VesselSegment<2> > p_segment_7 = VesselSegment<2>::Create(p_node_6, p_node_7);

            // Make vessels
            std::shared_ptr<Vessel<2> > p_vessel_1 = Vessel<2>::Create(p_segment_1);
            std::shared_ptr<Vessel<2> > p_vessel_2 = Vessel<2>::Create(p_segment_2);
            std::shared_ptr<Vessel<2> > p_vessel_3 = Vessel<2>::Create(p_segment_3);
            std::shared_ptr<Vessel<2> > p_vessel_4 = Vessel<2>::Create(p_segment_4);
            std::shared_ptr<Vessel<2> > p_vessel_5 = Vessel<2>::Create(p_segment_5);
            std::shared_ptr<Vessel<2> > p_vessel_6 = Vessel<2>::Create(p_segment_6);
            std::shared_ptr<Vessel<2> > p_vessel_7 = Vessel<2>::Create(p_segment_7);

            // Add the vessels to a vessel network
            std::shared_ptr<VesselNetwork<2> > p_network = VesselNetwork<2>::Create();
            p_network->AddVessel(p_vessel_1);
            p_network->AddVessel(p_vessel_2);
            p_network->AddVessel(p_vessel_3);
            p_network->AddVessel(p_vessel_4);
            p_network->AddVessel(p_vessel_5);
            p_network->AddVessel(p_vessel_6);
            p_network->AddVessel(p_vessel_7);

            // No pruning
            std::ostringstream strs;
            strs << std::fixed << std::setprecision( 1 );
            strs << "TestEquilateralNetwork/DirichletHaematocrit/Alpha" << alpha << "/NoPruning/";
            std::string str_directory_name = strs.str();
            auto p_output_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);         

            // Prune lower path (literally just removes vessels)
            // auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestEquilateralNetwork/LineSource2D/VTK/PruneUpperPath");   
            // p_network->RemoveVessel(p_vessel_2,true);
            // p_network->RemoveVessel(p_vessel_4,true);
            // p_network->RemoveVessel(p_vessel_7,true);

            // Prune upper path (literally just removes vessels)
            // auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestEquilateralNetwork/LineSource2D/VTK/PruneLowerPath");  
            // p_network->RemoveVessel(p_vessel_3,true);
            // p_network->RemoveVessel(p_vessel_5,true);
            // p_network->RemoveVessel(p_vessel_7,true);

            // Set up the grid for the finite difference solver
            auto p_grid = RegularGrid<2>::Create();
            p_grid->SetSpacing(equilateral_grid_spacing);
            c_vector<unsigned, 3> dimensions;
            dimensions[0] = unsigned((domain_x)/(equilateral_grid_spacing))+1; // num x
            dimensions[1] = unsigned((domain_y)/(equilateral_grid_spacing))+1; // num_y
            dimensions[2] = 1;
            p_grid->SetDimensions(dimensions);

            // Choose the PDE
            std::shared_ptr<DiscreteContinuumLinearEllipticPde<2> > p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();

            // Set the diffusivity and decay terms
            p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
            p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));
    
            // Set the O2 concentration value
            QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") * GenericParameters::mpGasConcentrationAtStp->GetValue("User");
            QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp * Owen11Parameters::mpReferencePartialPressure->GetValue("User");
            std::cout << vessel_oxygen_concentration << std::endl; 

            // Set the boundary condition to be the network acting as a Dirichlet source
            std::shared_ptr<DiscreteContinuumBoundaryCondition<2> > p_vessel_boundary_condition = DiscreteContinuumBoundaryCondition<2>::Create();
            p_vessel_boundary_condition->SetValue(vessel_oxygen_concentration);  // set the O2 concentration at the source
            p_vessel_boundary_condition->SetType(BoundaryConditionType::VESSEL_LINE);
            p_vessel_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);

            // Set up the finite difference solver for oxygen (which handles everything)
            SimpleLinearEllipticFiniteDifferenceSolver<2> solver;
            solver.SetGrid(p_grid);
            solver.SetPde(p_oxygen_pde);
            solver.AddBoundaryCondition(p_vessel_boundary_condition);
            solver.SetVesselNetwork(p_network);
            solver.SetLabel("oxygen");
            solver.SetFileName("oxygen_solution_0");

            // Run the solver and write the output file (visualise with Paraview; set Filters->Alphabetical->Tube)
            solver.SetFileHandler(p_output_file_handler);
            solver.SetWriteSolution(true);
            solver.Solve();
            p_network->Write(p_output_file_handler->GetOutputDirectoryFullPath() + "equilateral_network_results.vtp");
    
            // Dump our parameter collection to an xml file and, importantly, clear it for the next test
            ParameterCollection::Instance()->DumpToFile(p_output_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
            ParameterCollection::Instance()->Destroy();
            BaseUnits::Instance()->Destroy();
        }
    }

    // Make an equilateral network on a PDE grid that has flow with different haematocrit splitting rules
    void xTestEquilateralNetworkWithFlow2D()
    {
        // Run the simulation with different heterogeneities
        for (unsigned n_alpha=4; n_alpha<=4; n_alpha++)
        {
            double alpha = 1.0+(double)n_alpha*0.1;  // alpha determines the relative radius of the left vessel 
            
            // Run the simulation with different solvers of interest
            for (unsigned h_solver=1; h_solver<=1; h_solver++)
            {    
                // Set file name based on haematocrit solver
                std::ostringstream strs;
                strs << std::fixed << std::setprecision( 1 );
                if (h_solver==1)
                {
                    strs << "TestEquilateralNetwork/ConstantHaematocrit/Alpha" << alpha << "/NoPruning/";
                }
                else if (h_solver==2)
                {
                    strs << "TestEquilateralNetwork/PriesHaematocrit/Alpha" << alpha << "/NoPruning/";
                    
                }
                else if (h_solver==3)
                {
                    strs << "TestEquilateralNetwork/MemoryHaematocrit/Alpha" << alpha << "/NoPruning/";
                }
                else if (h_solver==4)
                {
                    strs << "TestEquilateralNetwork/FungHaematocrit/Alpha" << alpha << "/NoPruning/";
                }
                else if (h_solver==5)
                {
                    strs << "TestEquilateralNetwork/YangHaematocrit/Alpha" << alpha << "/NoPruning/";
                }
                std::string str_directory_name = strs.str();
                auto p_output_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                // Set up the reference length for the simulation
                QLength reference_length(1.0_um);
                BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

                // Set key vessel parameters
                QLength vessel_length(10.0_um);
                QLength vessel_height = (pow(2,0.5)*vessel_length)*0.5;

                // Set up the domain parameters
                QLength domain_x = vessel_height+vessel_height+vessel_length+vessel_length+0.0001_um;  // this should be the x-position of outlet node 
                QLength domain_y = domain_x;  // should be the same as domain_x to make square domain
                QLength mid_domain_y = domain_y*0.5;

                // Set nodes based on an equilateral network
                std::shared_ptr<VesselNode<2> > p_node_1 = VesselNode<2>::Create(0.0_um, mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_2 = VesselNode<2>::Create(vessel_length, mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_3 = VesselNode<2>::Create(vessel_height+vessel_length, vessel_height+mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_4 = VesselNode<2>::Create(vessel_height+vessel_length, -vessel_height+mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_5 = VesselNode<2>::Create(vessel_height+vessel_height+vessel_length, mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_6 = VesselNode<2>::Create(vessel_height+vessel_height+vessel_length+vessel_length, mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_7 = VesselNode<2>::Create(domain_x, mid_domain_y);  // add a node at a distance of 100 picometres from the last node to work around glitch of disconnected flow when setting p_node_6 as the output node

                // Make segments 
                std::shared_ptr<VesselSegment<2> > p_segment_1 = VesselSegment<2>::Create(p_node_1, p_node_2);
                std::shared_ptr<VesselSegment<2> > p_segment_2 = VesselSegment<2>::Create(p_node_2, p_node_3);
                std::shared_ptr<VesselSegment<2> > p_segment_3 = VesselSegment<2>::Create(p_node_2, p_node_4);
                std::shared_ptr<VesselSegment<2> > p_segment_4 = VesselSegment<2>::Create(p_node_3, p_node_5);
                std::shared_ptr<VesselSegment<2> > p_segment_5 = VesselSegment<2>::Create(p_node_4, p_node_5);
                std::shared_ptr<VesselSegment<2> > p_segment_6 = VesselSegment<2>::Create(p_node_5, p_node_6);
                std::shared_ptr<VesselSegment<2> > p_segment_7 = VesselSegment<2>::Create(p_node_6, p_node_7);

                // Make vessels
                std::shared_ptr<Vessel<2> > p_vessel_1 = Vessel<2>::Create(p_segment_1);
                std::shared_ptr<Vessel<2> > p_vessel_2 = Vessel<2>::Create(p_segment_2);
                std::shared_ptr<Vessel<2> > p_vessel_3 = Vessel<2>::Create(p_segment_3);
                std::shared_ptr<Vessel<2> > p_vessel_4 = Vessel<2>::Create(p_segment_4);
                std::shared_ptr<Vessel<2> > p_vessel_5 = Vessel<2>::Create(p_segment_5);
                std::shared_ptr<Vessel<2> > p_vessel_6 = Vessel<2>::Create(p_segment_6);
                std::shared_ptr<Vessel<2> > p_vessel_7 = Vessel<2>::Create(p_segment_7);

                // Add the vessels to a vessel network
                std::shared_ptr<VesselNetwork<2> > p_network = VesselNetwork<2>::Create();
                p_network->AddVessel(p_vessel_1);
                p_network->AddVessel(p_vessel_2);
                p_network->AddVessel(p_vessel_3);
                p_network->AddVessel(p_vessel_4);
                p_network->AddVessel(p_vessel_5);
                p_network->AddVessel(p_vessel_6);
                p_network->AddVessel(p_vessel_7);

                // No pruning
                // auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestEquilateralNetwork/VTK/NoPruning/");            

                // Prune lower path (literally just removes vessels)
                // auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestEquilateralNetwork/VTK/PruneUpperPath");   
                // p_network->RemoveVessel(p_vessel_2,true);
                // p_network->RemoveVessel(p_vessel_4,true);
                // p_network->RemoveVessel(p_vessel_7,true);

                // Prune upper path (literally just removes vessels)
                // auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestEquilateralNetwork/VTK/PruneLowerPath");  
                // p_network->RemoveVessel(p_vessel_3,true);
                // p_network->RemoveVessel(p_vessel_5,true);
                // p_network->RemoveVessel(p_vessel_7,true);
                
                // Specify which nodes are the inlets and outlets
                p_network->GetNode(0)->GetFlowProperties()->SetIsInputNode(true);
                p_network->GetNode(0)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));
                p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetIsOutputNode(true);
                p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));

                // Set segment radii values
                QLength vessel_radius(1.0 *GenericParameters::mpCapillaryRadius->GetValue());
                VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);

                // Set heterogeneous radii for upper path
                QLength alpha_radius(alpha*vessel_radius);
                p_vessel_2->SetRadius(alpha_radius);
                p_vessel_4->SetRadius(alpha_radius);

                // Set segment viscosity values
                QDynamicViscosity viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();
                auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                p_viscosity_calculator->SetVesselNetwork(p_network);
                p_viscosity_calculator->Calculate();
                // VesselNetworkPropertyManager<2>::SetSegmentViscosity(p_network, viscosity);

                // Set up the impedance calculator
                VesselImpedanceCalculator<2> impedance_calculator = VesselImpedanceCalculator<2>();
                impedance_calculator.SetVesselNetwork(p_network);
                impedance_calculator.Calculate();

                // Set up the flow solver
                FlowSolver<2> flow_solver = FlowSolver<2>();
                flow_solver.SetVesselNetwork(p_network);
                flow_solver.SetUp();
                flow_solver.Solve();

                // Set the haematocrit for all vessels
                // ConstantHaematocritSolver<2> haematocrit_solver = ConstantHaematocritSolver<2>();
                // haematocrit_solver.SetVesselNetwork(p_network);
                // haematocrit_solver.SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                // haematocrit_solver.Calculate();

                // Set the haematocrit solver
                std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                if (h_solver==1)
                {
                    std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;       
                }
                else if (h_solver==2)
                {
                    std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;                    
                }
                else if (h_solver==3)
                {
                    std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;      
                }
                else if (h_solver==4)
                {
                    std::cout << "Now using FungHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;      
                }
                else if (h_solver==5)
                {
                    std::cout << "Now using YangHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;      
                }
                // p_abstract_haematocrit_solver->Calculate();

                // // Set up the grid
                // std::shared_ptr<Part<2> > p_domain = Part<2>::Create();
                // p_domain->AddRectangle(domain_x, domain_y);
                // std::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
                // // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
                // // QLength grid_spacing(10_um);  // the simulation time gets quite long if you reduce the resolution further
                // p_grid->GenerateFromPart(p_domain, equilateral_grid_spacing);  // set domain and grid spacing

                // Set up the grid for the finite difference solver
                auto p_grid = RegularGrid<2>::Create();
                p_grid->SetSpacing(equilateral_grid_spacing);
                c_vector<unsigned, 3> dimensions;
                dimensions[0] = unsigned((domain_x)/(equilateral_grid_spacing))+1; // num x
                dimensions[1] = unsigned((domain_y)/(equilateral_grid_spacing))+1; // num_y
                dimensions[2] = 1;
                p_grid->SetDimensions(dimensions);
                
                // Choose the PDE
                std::shared_ptr<DiscreteContinuumLinearEllipticPde<2> > p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                
                // Set the diffusivity and decay terms
                p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                // Set up the discrete source
                std::shared_ptr<VesselBasedDiscreteSource<2> > p_vessel_source = VesselBasedDiscreteSource<2>::Create();
                QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                    GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                    Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                p_vessel_source->SetReferenceConcentration(vessel_oxygen_concentration);
                p_vessel_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_vessel_source->SetVesselPermeability(1.0*Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                p_oxygen_pde->AddDiscreteSource(p_vessel_source);

                // Set up the finite difference solver for oxygen (which handles everything)
                auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                p_oxygen_solver->SetPde(p_oxygen_pde);
                p_oxygen_solver->SetLabel("oxygen");
                p_oxygen_solver->SetGrid(p_grid);

                // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                double initial_haematocrit = 1;
                unsigned max_iter = 1000; 
                double tolerance2 = 1.e-10;
                std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                for(unsigned idx=0;idx<max_iter;idx++)
                {
                    // Run the solvers
                    impedance_calculator.Calculate();
                    flow_solver.SetUp();
                    flow_solver.Solve();
                    p_abstract_haematocrit_solver->Calculate();
                    p_viscosity_calculator->Calculate();

                    // Get the residual
                    double max_difference = 0.0;
                    double h_for_max = 0.0;
                    double prev_for_max = 0.0;
                    for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                    {
                        double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                        double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                        if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                        {
                            max_difference = difference;
                            h_for_max = current_haematocrit;
                            prev_for_max = previous_haematocrit[jdx];
                        }
                        previous_haematocrit[jdx] = current_haematocrit;
                    }
                    std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                    // Print the final or intermediary convergence results
                    if(max_difference<=tolerance2)  
                    {
                        std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                        break;
                    }
                    else
                    {
                        if(idx%1==0)
                        {
                            std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                            // std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                            // std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                            // p_network->Write(output_file);
                        }
                    }

                    // If there is no convergence after all the iterations, print the error message.
                    if(idx==max_iter-1)
                    {
                        std::cout << "Problem encountered in Equilateral Network with h_solver = " << h_solver << " using alpha = " << n_alpha << std::endl;
                        EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                    }
                }
                
                // Run the simulation 
                SimulationTime::Instance()->SetStartTime(0.0);
                SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                p_microvessel_solver->SetVesselNetwork(p_network);
                p_microvessel_solver->SetOutputFileHandler(p_output_file_handler);
                p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                p_microvessel_solver->Run();

                // Print the average oxygenation
                std::vector<double> solution = p_oxygen_solver->GetSolution();
                double average_oxygen = 0.0;
                for(unsigned jdx=0;jdx<solution.size();jdx++)
                {
                    average_oxygen += solution[jdx];
                }
                average_oxygen /= double(solution.size());
                std::cout << "Average oxygen: " << average_oxygen << std::endl;
                
                // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                std::string output_file = p_output_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                p_network->Write(output_file);

                // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                ParameterCollection::Instance()->DumpToFile(p_output_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                ParameterCollection::Instance()->Destroy();
                BaseUnits::Instance()->Destroy();
                SimulationTime::Instance()->Destroy();

                // // Run the solver and write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                // solver.SetFileHandler(p_output_file_handler);
                // solver.SetWriteSolution(true);
                // solver.Solve();
                // p_network->Write(p_output_file_handler->GetOutputDirectoryFullPath() + "equilateral_network_results.vtp");

                // // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                // ParameterCollection::Instance()->DumpToFile(p_output_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                // ParameterCollection::Instance()->Destroy();
                // BaseUnits::Instance()->Destroy();
            }
        }
    }    

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Debugging Networks
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Change resolution of custom network for better visualisation
    QLength irregular_grid_spacing = 1.0_um;

    // Make a forking network with one daughter vessel longer than the other
    void xTestLongVesselNetworkWithFlow2D()
    {
        // Change max. thickness heterogeneity
        unsigned max_alpha = 0;
   
        // Run the simulation with different heterogeneities
        for (unsigned n_alpha=0; n_alpha<=max_alpha; n_alpha++)
        {
            double alpha = 1.0+(double)n_alpha*0.1;  // alpha determines the relative radius of thicker vessels
            
            // Run the simulation with different solvers of interest
            for (unsigned h_solver=1; h_solver<6; h_solver++)
            {    
                // Set file name based on haematocrit solver
                std::ostringstream strs;
                strs << std::fixed << std::setprecision( 1 );
                strs << "TestDebuggingNetworks/LongVesselNetwork/";
                if (h_solver==1)
                {
                    strs << "ConstantHaematocrit";
                }
                else if (h_solver==2)
                {
                    strs << "PriesHaematocrit";
                    
                }
                else if (h_solver==3)
                {
                    strs << "MemoryHaematocrit";
                }
                else if (h_solver==4)
                {
                    strs << "FungHaematocrit";
                }
                else if (h_solver==5)
                {
                    strs << "YangHaematocrit";
                }
                strs << "/Alpha" << alpha << "/NoPruning/";
                std::string str_directory_name = strs.str();
                auto p_output_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                // Set up the reference length for the simulation
                QLength reference_length(1.0_um);
                BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

                // Set key vessel parameters
                QLength vessel_length(100.0_um);
                QLength vessel_height = (pow(2,0.5)*vessel_length)*0.5;

                // Set up the domain parameters
                QLength domain_x = vessel_height+vessel_height+vessel_length+vessel_length+0.0001_um;  // this should be the x-position of outlet node 
                QLength domain_y = domain_x;  // should be the same as domain_x to make square domain
                QLength mid_domain_y = domain_y*0.5;

                // Set nodes based on desired network
                std::shared_ptr<VesselNode<2> > p_node_1 = VesselNode<2>::Create(0.0_um, mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_2 = VesselNode<2>::Create(vessel_length, mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_3 = VesselNode<2>::Create(vessel_height+vessel_length, vessel_height+mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_4 = VesselNode<2>::Create(vessel_height+vessel_length, -vessel_height+mid_domain_y-25.0_um);
                std::shared_ptr<VesselNode<2> > p_node_5 = VesselNode<2>::Create(vessel_height+vessel_height+vessel_length, mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_6 = VesselNode<2>::Create(vessel_height+vessel_height+vessel_length+vessel_length, mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_7 = VesselNode<2>::Create(domain_x, mid_domain_y);  // add a node at a distance of 100 picometres from the last node to work around glitch of disconnected flow when setting p_node_6 as the output node

                // Make segments 
                std::shared_ptr<VesselSegment<2> > p_segment_1 = VesselSegment<2>::Create(p_node_1, p_node_2);
                std::shared_ptr<VesselSegment<2> > p_segment_2 = VesselSegment<2>::Create(p_node_2, p_node_3);
                std::shared_ptr<VesselSegment<2> > p_segment_3 = VesselSegment<2>::Create(p_node_2, p_node_4);
                std::shared_ptr<VesselSegment<2> > p_segment_4 = VesselSegment<2>::Create(p_node_3, p_node_5);
                std::shared_ptr<VesselSegment<2> > p_segment_5 = VesselSegment<2>::Create(p_node_4, p_node_5);
                std::shared_ptr<VesselSegment<2> > p_segment_6 = VesselSegment<2>::Create(p_node_5, p_node_6);
                std::shared_ptr<VesselSegment<2> > p_segment_7 = VesselSegment<2>::Create(p_node_6, p_node_7);

                // Make vessels
                std::shared_ptr<Vessel<2> > p_vessel_1 = Vessel<2>::Create(p_segment_1);
                std::shared_ptr<Vessel<2> > p_vessel_2 = Vessel<2>::Create(p_segment_2);
                std::shared_ptr<Vessel<2> > p_vessel_3 = Vessel<2>::Create(p_segment_3);
                std::shared_ptr<Vessel<2> > p_vessel_4 = Vessel<2>::Create(p_segment_4);
                std::shared_ptr<Vessel<2> > p_vessel_5 = Vessel<2>::Create(p_segment_5);
                std::shared_ptr<Vessel<2> > p_vessel_6 = Vessel<2>::Create(p_segment_6);
                std::shared_ptr<Vessel<2> > p_vessel_7 = Vessel<2>::Create(p_segment_7);

                // Add the vessels to a vessel network
                std::shared_ptr<VesselNetwork<2> > p_network = VesselNetwork<2>::Create();
                p_network->AddVessel(p_vessel_1);
                p_network->AddVessel(p_vessel_2);
                p_network->AddVessel(p_vessel_3);
                p_network->AddVessel(p_vessel_4);
                p_network->AddVessel(p_vessel_5);
                p_network->AddVessel(p_vessel_6);
                p_network->AddVessel(p_vessel_7);

                // No pruning
                // auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestEquilateralNetwork/VTK/NoPruning/");            

                // Prune lower path (literally just removes vessels)
                // auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestEquilateralNetwork/VTK/PruneUpperPath");   
                // p_network->RemoveVessel(p_vessel_2,true);
                // p_network->RemoveVessel(p_vessel_4,true);
                // p_network->RemoveVessel(p_vessel_7,true);

                // Prune upper path (literally just removes vessels)
                // auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestEquilateralNetwork/VTK/PruneLowerPath");  
                // p_network->RemoveVessel(p_vessel_3,true);
                // p_network->RemoveVessel(p_vessel_5,true);
                // p_network->RemoveVessel(p_vessel_7,true);
                
                // Specify which nodes are the inlets and outlets
                p_network->GetNode(0)->GetFlowProperties()->SetIsInputNode(true);
                p_network->GetNode(0)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));
                p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetIsOutputNode(true);
                p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));

                // Set segment radii values
                QLength vessel_radius(1.0 *GenericParameters::mpCapillaryRadius->GetValue());
                VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);

                // Set heterogeneous radii for desired vessels
                // QLength alpha_radius(alpha*vessel_radius);
                // p_vessel_2->SetRadius(alpha_radius);
                // p_vessel_4->SetRadius(alpha_radius);

                // Set segment viscosity values
                QDynamicViscosity viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();
                auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                p_viscosity_calculator->SetVesselNetwork(p_network);
                p_viscosity_calculator->Calculate();
                // VesselNetworkPropertyManager<2>::SetSegmentViscosity(p_network, viscosity);

                // Set up the impedance calculator
                VesselImpedanceCalculator<2> impedance_calculator = VesselImpedanceCalculator<2>();
                impedance_calculator.SetVesselNetwork(p_network);
                impedance_calculator.Calculate();

                // Set up the flow solver
                FlowSolver<2> flow_solver = FlowSolver<2>();
                flow_solver.SetVesselNetwork(p_network);
                flow_solver.SetUp();
                flow_solver.Solve();

                // Set the haematocrit for all vessels
                // ConstantHaematocritSolver<2> haematocrit_solver = ConstantHaematocritSolver<2>();
                // haematocrit_solver.SetVesselNetwork(p_network);
                // haematocrit_solver.SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                // haematocrit_solver.Calculate();

                // Set the haematocrit solver
                std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                if (h_solver==1)
                {
                    std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;       
                }
                else if (h_solver==2)
                {
                    std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;                    
                }
                else if (h_solver==3)
                {
                    std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;      
                }
                else if (h_solver==4)
                {
                    std::cout << "Now using FungHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;      
                }
                else if (h_solver==5)
                {
                    std::cout << "Now using YangHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;      
                }

                // Set up the grid for the finite difference solver
                auto p_grid = RegularGrid<2>::Create();
                p_grid->SetSpacing(irregular_grid_spacing);
                c_vector<unsigned, 3> dimensions;
                dimensions[0] = unsigned((domain_x)/(irregular_grid_spacing))+1; // num x
                dimensions[1] = unsigned((domain_y)/(irregular_grid_spacing))+1; // num_y
                dimensions[2] = 1;
                p_grid->SetDimensions(dimensions);
                
                // Choose the PDE
                std::shared_ptr<DiscreteContinuumLinearEllipticPde<2> > p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                
                // Set the diffusivity and decay terms
                p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                // Set up the discrete source
                std::shared_ptr<VesselBasedDiscreteSource<2> > p_vessel_source = VesselBasedDiscreteSource<2>::Create();
                QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                    GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                    Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                p_vessel_source->SetReferenceConcentration(vessel_oxygen_concentration);
                p_vessel_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_vessel_source->SetVesselPermeability(1.0*Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                p_oxygen_pde->AddDiscreteSource(p_vessel_source);

                // Set up the finite difference solver for oxygen (which handles everything)
                auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                p_oxygen_solver->SetPde(p_oxygen_pde);
                p_oxygen_solver->SetLabel("oxygen");
                p_oxygen_solver->SetGrid(p_grid);

                // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                double initial_haematocrit = 0.45;
                unsigned max_iter = 1000; 
                double tolerance2 = 1.e-10;
                std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                for(unsigned idx=0;idx<max_iter;idx++)
                {
                    // Run the solvers
                    impedance_calculator.Calculate();
                    flow_solver.SetUp();
                    flow_solver.Solve();
                    p_abstract_haematocrit_solver->Calculate();
                    p_viscosity_calculator->Calculate();

                    // Get the residual
                    double max_difference = 0.0;
                    double h_for_max = 0.0;
                    double prev_for_max = 0.0;
                    for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                    {
                        double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                        double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                        if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                        {
                            max_difference = difference;
                            h_for_max = current_haematocrit;
                            prev_for_max = previous_haematocrit[jdx];
                        }
                        previous_haematocrit[jdx] = current_haematocrit;
                    }
                    std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                    // Print the final or intermediary convergence results
                    if(max_difference<=tolerance2)  
                    {
                        std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                        break;
                    }
                    else
                    {
                        if(idx%1==0)
                        {
                            std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                            // std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                            // std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                            // p_network->Write(output_file);
                        }
                    }

                    // If there is no convergence after all the iterations, print the error message.
                    if(idx==max_iter-1)
                    {
                        std::cout << "Problem encountered in Irregular Network with h_solver = " << h_solver << " using alpha = " << n_alpha << std::endl;
                        EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                    }
                }
                
                // Run the simulation 
                SimulationTime::Instance()->SetStartTime(0.0);
                SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                p_microvessel_solver->SetVesselNetwork(p_network);
                p_microvessel_solver->SetOutputFileHandler(p_output_file_handler);
                p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                p_microvessel_solver->Run();

                // Print the average oxygenation
                std::vector<double> solution = p_oxygen_solver->GetSolution();
                double average_oxygen = 0.0;
                for(unsigned jdx=0;jdx<solution.size();jdx++)
                {
                    average_oxygen += solution[jdx];
                }
                average_oxygen /= double(solution.size());
                std::cout << "Average oxygen: " << average_oxygen << std::endl;
                
                // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                std::string output_file = p_output_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                p_network->Write(output_file);

                // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                ParameterCollection::Instance()->DumpToFile(p_output_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                ParameterCollection::Instance()->Destroy();
                BaseUnits::Instance()->Destroy();
                SimulationTime::Instance()->Destroy();
            }
        }
    }    

    // Make a 2D hexagonal network on a PDE grid with flow and H-splitting (single inlet and outlet)
    void xTestSingleFeedHexagonalUnitWithFlow2D()
    {
        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Define the key parameters
        QLength large_vessel_radius = 10_um;
        double dimless_domain_size_x = 400.0;  // x-coordinate of output node
        double dimless_domain_size_y = 346.41 + 173.205;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
        unsigned dimless_vessel_length = 100.0;
        std::vector<std::vector<unsigned> > Order;
        std::shared_ptr<VesselNetwork<2> > p_network;
        std::vector<std::shared_ptr<Vessel<2> > > vessels;
        string line2;

        // Initialise the simulation space
        QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
        QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
        std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
        std::ofstream outfileMean;

        // Define the key parameters
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        double initial_haematocrit = 0.45;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";

        // Read the network from a file
        VesselNetworkGenerator<2> network_generator;
        std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix_SingleFeed.txt");
        std::vector<std::vector<double> > rEdgesMatrix;
        string line;
        
        // Break down the rows in the file into column values
        while (std::getline(in, line)) 
        {
            rEdgesMatrix.push_back(std::vector<double>());
            std::stringstream split(line);
            double value;
            while (split >> value)
            {
                rEdgesMatrix.back().push_back(value);
            }
        }

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<6; h_solver++)
        {
            // Set file name
            std::ostringstream strs;
            strs << std::fixed << std::setprecision( 2 );
            strs << "TestDebuggingNetworks/SingleFeedHexagonalUnit/NoPruning/";        

            // std::ostringstream strs;
            if (h_solver==1)
            {
                strs << "ConstantHaematocrit";
            }
            else if (h_solver==2)
            {
                strs << "PriesHaematocrit";
            }
            else if (h_solver==3)
            {
                strs << "MemoryHaematocrit";
            }
            else if (h_solver==4)
            {
                strs << "FungHaematocrit";
            }
            else if (h_solver==5)
            {
                strs << "YangHaematocrit";
            }
            std::string str_directory_name = strs.str();
            auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

            // Generate the network
            p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);
            
            // Set inlet and outlet nodes
            auto p_segment = p_network->GetVesselSegments()[0];
            vessels = p_network->GetVessels();
            for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
            {
                if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
                {
                    if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)
                    {
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
                    }
                    if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
                    {
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
                    }
                }
                if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)   	
                {
                    if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)
                    {
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
                    }
                    if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
                    {
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
                    }
                }
            }
            p_segment = p_network->GetVesselSegments()[0];
            p_segment->SetRadius(large_vessel_radius);
            VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   

            // Set the haematocrit solver
            std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
            if (h_solver==1)
            {
                std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;       
            }
            else if (h_solver==2)
            {
                std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;                    
            }
            else if (h_solver==3)
            {
                std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;      
            }
            else if (h_solver==4)
            {
                std::cout << "Now using FungHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;      
            }
            else if (h_solver==5)
            {
                std::cout << "Now using YangHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;      
            }

            // Set up the viscosity solver
            auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
            p_viscosity_calculator->SetPlasmaViscosity(viscosity);
            p_viscosity_calculator->SetVesselNetwork(p_network);
            p_viscosity_calculator->Calculate();
            
            // Set up the impedance solver
            auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
            p_impedance_calculator->SetVesselNetwork(p_network);
            p_impedance_calculator->Calculate();

            // Set up the flow solver 
            FlowSolver<2> flow_solver;
            flow_solver.SetVesselNetwork(p_network);
            // flow_solver.SetUp();
            flow_solver.SetUseDirectSolver(true);
            // flow_solver.Solve();

            // Set up the grid for the finite difference solver
            auto p_grid = RegularGrid<2>::Create();
            p_grid->SetSpacing(irregular_grid_spacing);
            c_vector<unsigned, 3> dimensions;
            dimensions[0] = unsigned((domain_side_length_x)/(irregular_grid_spacing))+1; // num x
            dimensions[1] = unsigned((domain_side_length_y)/(irregular_grid_spacing))+1; // num_y
            dimensions[2] = 1;
            p_grid->SetDimensions(dimensions);

            // Choose the PDE
            auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
            
            // Set the diffusivity and decay terms
            p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
            p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

            // Set up the discrete source
            auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
            QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                    GenericParameters::mpGasConcentrationAtStp->GetValue("User");
            QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                    Owen11Parameters::mpReferencePartialPressure->GetValue("User");
            p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
            p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
            p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
            p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

            // Set up the finite difference solver for oxygen (which handles everything)
            auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
            p_oxygen_solver->SetPde(p_oxygen_pde);
            p_oxygen_solver->SetLabel("oxygen");
            p_oxygen_solver->SetGrid(p_grid);
            // solver.SetFileName("oxygen_solution_0");

            // Set up flag for broken solver
            unsigned broken_solver = 0;

            // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
            unsigned max_iter = 1000;  // 1000 
            double tolerance2 = 1.e-10;
            std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
            std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
            for(unsigned idx=0;idx<max_iter;idx++)
            {
                // Run the solvers
                p_impedance_calculator->Calculate();
                flow_solver.SetUp();
                flow_solver.Solve();
                p_abstract_haematocrit_solver->Calculate();
                p_viscosity_calculator->Calculate();

                // Get the residual
                double max_difference = 0.0;
                double h_for_max = 0.0;
                double prev_for_max = 0.0;
                for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                {
                    double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                    double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                    if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                    {
                        max_difference = difference;
                        h_for_max = current_haematocrit;
                        prev_for_max = previous_haematocrit[jdx];
                    }
                    previous_haematocrit[jdx] = current_haematocrit;
                }
                std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                // Print the final or intermediary convergence results
                if(max_difference<=tolerance2)  
                {
                    std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                    broken_solver = 0;
                    break;
                }
                else
                {
                    if(idx%1==0)
                    {
                        std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                        std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                        p_network->Write(output_file);
                    }
                }

                // If there is no convergence after all the iterations, print the error message.
                if(idx==max_iter-1)
                {
                    std::cout << "Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << std::endl;
                    error_log << "\n Problem encountered in hexagonal network with h_solver = " << to_string(h_solver); 
                    broken_solver = 1;
                    break;
                    // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                }
            }

            // If solver doesn't converge, move on to next one
            if (broken_solver == 1)
            {
                continue;
            }

            // Run the simulation 
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
            auto p_microvessel_solver = MicrovesselSolver<2>::Create();
            p_microvessel_solver->SetVesselNetwork(p_network);
            p_microvessel_solver->SetOutputFileHandler(p_file_handler);
            p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
            p_microvessel_solver->Run();

            // Print the average oxygenation
            std::vector<double> solution = p_oxygen_solver->GetSolution();
            double average_oxygen = 0.0;
            for(unsigned jdx=0;jdx<solution.size();jdx++)
            {
                average_oxygen += solution[jdx];
            }
            average_oxygen /= double(solution.size());
            std::cout << "Average oxygen: " << average_oxygen << std::endl;
            
            // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
            std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
            p_network->Write(output_file);

            // Dump our parameter collection to an xml file and, importantly, clear it for the next test
            ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
            ParameterCollection::Instance()->Destroy();
            BaseUnits::Instance()->Destroy();
            SimulationTime::Instance()->Destroy();
        }
         
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a 2D hexagonal network on a PDE grid with flow and H-splitting (multiple inlets and outlets)
    void xTestMultiFeedHexagonalUnitWithFlow2D()
    {
        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Define the key parameters
        QLength large_vessel_radius = 10_um;
        double dimless_domain_size_x = 400.0;  // x-coordinate of output node
        double dimless_domain_size_y = 433.013 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
        unsigned dimless_vessel_length = 100.0;
        std::vector<std::vector<unsigned> > Order;
        std::shared_ptr<VesselNetwork<2> > p_network;
        std::vector<std::shared_ptr<Vessel<2> > > vessels;
        string line2;

        // Initialise the simulation space
        QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
        QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
        std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
        std::ofstream outfileMean;

        // Define the key parameters
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        double initial_haematocrit = 0.45;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";

        // Read the network from a file
        VesselNetworkGenerator<2> network_generator;
        std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix_MultiFeed.txt");
        std::vector<std::vector<double> > rEdgesMatrix;
        string line;
        
        // Break down the rows in the file into column values
        while (std::getline(in, line)) 
        {
            rEdgesMatrix.push_back(std::vector<double>());
            std::stringstream split(line);
            double value;
            while (split >> value)
            {
                rEdgesMatrix.back().push_back(value);
            }
        }

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<6; h_solver++)
        {
            // Set file name
            std::ostringstream strs;
            strs << std::fixed << std::setprecision( 2 );
            strs << "TestDebuggingNetworks/MultiFeedHexagonalUnit/NoPruning/";        

            // std::ostringstream strs;
            if (h_solver==1)
            {
                strs << "ConstantHaematocrit";
            }
            else if (h_solver==2)
            {
                strs << "PriesHaematocrit";
            }
            else if (h_solver==3)
            {
                strs << "MemoryHaematocrit";
            }
            else if (h_solver==4)
            {
                strs << "FungHaematocrit";
            }
            else if (h_solver==5)
            {
                strs << "YangHaematocrit";
            }
            std::string str_directory_name = strs.str();
            auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

            // Generate the network
            p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);
            
            // Set inlet and outlet nodes
            auto p_segment = p_network->GetVesselSegments()[0];
            vessels = p_network->GetVessels();
            for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
            {
                if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
                {
                    if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)
                    {
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
                    }
                    if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
                    {
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
                    }
                }
                if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)   	
                {
                    if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)
                    {
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
                    }
                    if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
                    {
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
                    }
                }
            }
            p_segment = p_network->GetVesselSegments()[0];
            p_segment->SetRadius(large_vessel_radius);
            VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   

            // Set the haematocrit solver
            std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
            if (h_solver==1)
            {
                std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;       
            }
            else if (h_solver==2)
            {
                std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;                    
            }
            else if (h_solver==3)
            {
                std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;      
            }
            else if (h_solver==4)
            {
                std::cout << "Now using FungHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;      
            }
            else if (h_solver==5)
            {
                std::cout << "Now using YangHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;      
            }

            // Set up the viscosity solver
            auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
            p_viscosity_calculator->SetPlasmaViscosity(viscosity);
            p_viscosity_calculator->SetVesselNetwork(p_network);
            p_viscosity_calculator->Calculate();
            
            // Set up the impedance solver
            auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
            p_impedance_calculator->SetVesselNetwork(p_network);
            p_impedance_calculator->Calculate();

            // Set up the flow solver 
            FlowSolver<2> flow_solver;
            flow_solver.SetVesselNetwork(p_network);
            // flow_solver.SetUp();
            flow_solver.SetUseDirectSolver(true);
            // flow_solver.Solve();

            // Set up the grid for the finite difference solver
            auto p_grid = RegularGrid<2>::Create();
            p_grid->SetSpacing(irregular_grid_spacing);
            c_vector<unsigned, 3> dimensions;
            dimensions[0] = unsigned((domain_side_length_x)/(irregular_grid_spacing))+1; // num x
            dimensions[1] = unsigned((domain_side_length_y)/(irregular_grid_spacing))+1; // num_y
            dimensions[2] = 1;
            p_grid->SetDimensions(dimensions);

            // Choose the PDE
            auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
            
            // Set the diffusivity and decay terms
            p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
            p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

            // Set up the discrete source
            auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
            QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                    GenericParameters::mpGasConcentrationAtStp->GetValue("User");
            QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                    Owen11Parameters::mpReferencePartialPressure->GetValue("User");
            p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
            p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
            p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
            p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

            // Set up the finite difference solver for oxygen (which handles everything)
            auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
            p_oxygen_solver->SetPde(p_oxygen_pde);
            p_oxygen_solver->SetLabel("oxygen");
            p_oxygen_solver->SetGrid(p_grid);
            // solver.SetFileName("oxygen_solution_0");

            // Set up flag for broken solver
            unsigned broken_solver = 0;

            // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
            unsigned max_iter = 1000;  // 1000 
            double tolerance2 = 1.e-10;
            std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
            std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
            for(unsigned idx=0;idx<max_iter;idx++)
            {
                // Run the solvers
                p_impedance_calculator->Calculate();
                flow_solver.SetUp();
                flow_solver.Solve();
                p_abstract_haematocrit_solver->Calculate();
                p_viscosity_calculator->Calculate();

                // Get the residual
                double max_difference = 0.0;
                double h_for_max = 0.0;
                double prev_for_max = 0.0;
                for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                {
                    double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                    double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                    if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                    {
                        max_difference = difference;
                        h_for_max = current_haematocrit;
                        prev_for_max = previous_haematocrit[jdx];
                    }
                    previous_haematocrit[jdx] = current_haematocrit;
                }
                std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                // Print the final or intermediary convergence results
                if(max_difference<=tolerance2)  
                {
                    std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;                    
                    broken_solver = 0;
                    break;
                }
                else
                {
                    if(idx%1==0)
                    {
                        std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                        std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                        p_network->Write(output_file);
                    }
                }

                // If there is no convergence after all the iterations, print the error message.
                if(idx==max_iter-1)
                {
                    std::cout << "Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << std::endl;
                    error_log << "\n Problem encountered in hexagonal network with h_solver = " << to_string(h_solver); 
                    broken_solver = 1;
                    break;
                    // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                }
            }

            // If solver doesn't converge, move on to next one
            if (broken_solver == 1)
            {
                continue;
            }

            // Run the simulation 
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
            auto p_microvessel_solver = MicrovesselSolver<2>::Create();
            p_microvessel_solver->SetVesselNetwork(p_network);
            p_microvessel_solver->SetOutputFileHandler(p_file_handler);
            p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
            p_microvessel_solver->Run();

            // Print the average oxygenation
            std::vector<double> solution = p_oxygen_solver->GetSolution();
            double average_oxygen = 0.0;
            for(unsigned jdx=0;jdx<solution.size();jdx++)
            {
                average_oxygen += solution[jdx];
            }
            average_oxygen /= double(solution.size());
            std::cout << "Average oxygen: " << average_oxygen << std::endl;
            
            // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
            std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
            p_network->Write(output_file);

            // Dump our parameter collection to an xml file and, importantly, clear it for the next test
            ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
            ParameterCollection::Instance()->Destroy();
            BaseUnits::Instance()->Destroy();
            SimulationTime::Instance()->Destroy();
        }
         
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make an equilateral network on a PDE grid that has flow with different haematocrit splitting rules and radii scaled by Murray's law
    void xTestMurrayNetworkWithFlow2D()
    {
        // Run the simulation with different heterogeneities
        for (unsigned n_alpha=0; n_alpha<=max_alpha; n_alpha++)
        {
            double alpha = 1.0+(double)n_alpha*0.1;  // alpha determines the relative radius of the left vessel 
            
            // Run the simulation with different solvers of interest
            for (unsigned h_solver=5; h_solver<6; h_solver++)
            {    
                // Set file name based on haematocrit solver
                std::ostringstream strs;
                strs << std::fixed << std::setprecision( 1 );
                if (h_solver==1)
                {
                    strs << "TestMurrayNetwork/ConstantHaematocrit/Alpha" << alpha << "/NoPruning/";
                }
                else if (h_solver==2)
                {
                    strs << "TestMurrayNetwork/PriesHaematocrit/Alpha" << alpha << "/NoPruning/";
                    
                }
                else if (h_solver==3)
                {
                    strs << "TestMurrayNetwork/MemoryHaematocrit/Alpha" << alpha << "/NoPruning/";
                }
                else if (h_solver==4)
                {
                    strs << "TestMurrayNetwork/FungHaematocrit/Alpha" << alpha << "/NoPruning/";
                }
                else if (h_solver==5)
                {
                    strs << "TestMurrayNetwork/YangHaematocrit/Alpha" << alpha << "/NoPruning/";
                }
                std::string str_directory_name = strs.str();
                auto p_output_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                // Set up the reference length for the simulation
                QLength reference_length(1.0_um);
                BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

                // Set key vessel parameters
                QLength vessel_length(100.0_um);
                QLength vessel_height = (pow(2,0.5)*vessel_length)*0.5;

                // Set up the domain parameters
                QLength domain_x = vessel_height+vessel_height+vessel_length+vessel_length+0.0001_um;  // this should be the x-position of outlet node 
                QLength domain_y = domain_x;  // should be the same as domain_x to make square domain
                QLength mid_domain_y = domain_y*0.5;

                // Set nodes based on an equilateral network
                std::shared_ptr<VesselNode<2> > p_node_1 = VesselNode<2>::Create(0.0_um, mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_2 = VesselNode<2>::Create(vessel_length, mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_3 = VesselNode<2>::Create(vessel_height+vessel_length, vessel_height+mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_4 = VesselNode<2>::Create(vessel_height+vessel_length, -vessel_height+mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_5 = VesselNode<2>::Create(vessel_height+vessel_height+vessel_length, mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_6 = VesselNode<2>::Create(vessel_height+vessel_height+vessel_length+vessel_length, mid_domain_y);
                std::shared_ptr<VesselNode<2> > p_node_7 = VesselNode<2>::Create(domain_x, mid_domain_y);  // add a node at a distance of 100 picometres from the last node to work around glitch of disconnected flow when setting p_node_6 as the output node

                // Make segments 
                std::shared_ptr<VesselSegment<2> > p_segment_1 = VesselSegment<2>::Create(p_node_1, p_node_2);
                std::shared_ptr<VesselSegment<2> > p_segment_2 = VesselSegment<2>::Create(p_node_2, p_node_3);
                std::shared_ptr<VesselSegment<2> > p_segment_3 = VesselSegment<2>::Create(p_node_2, p_node_4);
                std::shared_ptr<VesselSegment<2> > p_segment_4 = VesselSegment<2>::Create(p_node_3, p_node_5);
                std::shared_ptr<VesselSegment<2> > p_segment_5 = VesselSegment<2>::Create(p_node_4, p_node_5);
                std::shared_ptr<VesselSegment<2> > p_segment_6 = VesselSegment<2>::Create(p_node_5, p_node_6);
                std::shared_ptr<VesselSegment<2> > p_segment_7 = VesselSegment<2>::Create(p_node_6, p_node_7);

                // Make vessels
                std::shared_ptr<Vessel<2> > p_vessel_1 = Vessel<2>::Create(p_segment_1);
                std::shared_ptr<Vessel<2> > p_vessel_2 = Vessel<2>::Create(p_segment_2);
                std::shared_ptr<Vessel<2> > p_vessel_3 = Vessel<2>::Create(p_segment_3);
                std::shared_ptr<Vessel<2> > p_vessel_4 = Vessel<2>::Create(p_segment_4);
                std::shared_ptr<Vessel<2> > p_vessel_5 = Vessel<2>::Create(p_segment_5);
                std::shared_ptr<Vessel<2> > p_vessel_6 = Vessel<2>::Create(p_segment_6);
                std::shared_ptr<Vessel<2> > p_vessel_7 = Vessel<2>::Create(p_segment_7);

                // Add the vessels to a vessel network
                std::shared_ptr<VesselNetwork<2> > p_network = VesselNetwork<2>::Create();
                p_network->AddVessel(p_vessel_1);
                p_network->AddVessel(p_vessel_2);
                p_network->AddVessel(p_vessel_3);
                p_network->AddVessel(p_vessel_4);
                p_network->AddVessel(p_vessel_5);
                p_network->AddVessel(p_vessel_6);
                p_network->AddVessel(p_vessel_7);

                // No pruning
                // auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestEquilateralNetwork/VTK/NoPruning/");            

                // Prune lower path (literally just removes vessels)
                // auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestEquilateralNetwork/VTK/PruneUpperPath");   
                // p_network->RemoveVessel(p_vessel_2,true);
                // p_network->RemoveVessel(p_vessel_4,true);
                // p_network->RemoveVessel(p_vessel_7,true);

                // Prune upper path (literally just removes vessels)
                // auto p_output_file_handler = std::make_shared<OutputFileHandler>("TestEquilateralNetwork/VTK/PruneLowerPath");  
                // p_network->RemoveVessel(p_vessel_3,true);
                // p_network->RemoveVessel(p_vessel_5,true);
                // p_network->RemoveVessel(p_vessel_7,true);
                
                // Specify which nodes are the inlets and outlets
                p_network->GetNode(0)->GetFlowProperties()->SetIsInputNode(true);
                p_network->GetNode(0)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));
                p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetIsOutputNode(true);
                p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));

                // Set daughter segment radii values
                QLength vessel_radius(39.685_um);
                VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);

                // Set inlet parent radius value
                QLength inlet_radius(50.0_um);
                p_vessel_1->SetRadius(inlet_radius);   
                p_vessel_6->SetRadius(inlet_radius);     
                p_vessel_7->SetRadius(inlet_radius);      


                // Set heterogeneous radii for upper path
                QLength alpha_radius(alpha*vessel_radius);
                p_vessel_2->SetRadius(alpha_radius);
                p_vessel_4->SetRadius(alpha_radius);

                // Set segment viscosity values
                QDynamicViscosity viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();
                auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                p_viscosity_calculator->SetVesselNetwork(p_network);
                p_viscosity_calculator->Calculate();
                // VesselNetworkPropertyManager<2>::SetSegmentViscosity(p_network, viscosity);

                // Set up the impedance calculator
                VesselImpedanceCalculator<2> impedance_calculator = VesselImpedanceCalculator<2>();
                impedance_calculator.SetVesselNetwork(p_network);
                impedance_calculator.Calculate();

                // Set up the flow solver
                FlowSolver<2> flow_solver = FlowSolver<2>();
                flow_solver.SetVesselNetwork(p_network);
                flow_solver.SetUp();
                flow_solver.Solve();

                // Set the haematocrit for all vessels
                // ConstantHaematocritSolver<2> haematocrit_solver = ConstantHaematocritSolver<2>();
                // haematocrit_solver.SetVesselNetwork(p_network);
                // haematocrit_solver.SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                // haematocrit_solver.Calculate();

                // Set the haematocrit solver
                std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                if (h_solver==1)
                {
                    std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;       
                }
                else if (h_solver==2)
                {
                    std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;                    
                }
                else if (h_solver==3)
                {
                    std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;      
                }
                else if (h_solver==4)
                {
                    std::cout << "Now using FungHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;      
                }
                else if (h_solver==5)
                {
                    std::cout << "Now using YangHaematocritSolver..." << std::endl;
                    auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                    p_haematocrit_solver->SetVesselNetwork(p_network);
                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_abstract_haematocrit_solver = p_haematocrit_solver;      
                }
                // p_abstract_haematocrit_solver->Calculate();

                // // Set up the grid
                // std::shared_ptr<Part<2> > p_domain = Part<2>::Create();
                // p_domain->AddRectangle(domain_x, domain_y);
                // std::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
                // // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
                // // QLength grid_spacing(10_um);  // the simulation time gets quite long if you reduce the resolution further
                // p_grid->GenerateFromPart(p_domain, equilateral_grid_spacing);  // set domain and grid spacing

                // Set up the grid for the finite difference solver
                auto p_grid = RegularGrid<2>::Create();
                p_grid->SetSpacing(equilateral_grid_spacing);
                c_vector<unsigned, 3> dimensions;
                dimensions[0] = unsigned((domain_x)/(equilateral_grid_spacing))+1; // num x
                dimensions[1] = unsigned((domain_y)/(equilateral_grid_spacing))+1; // num_y
                dimensions[2] = 1;
                p_grid->SetDimensions(dimensions);
                
                // Choose the PDE
                std::shared_ptr<DiscreteContinuumLinearEllipticPde<2> > p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                
                // Set the diffusivity and decay terms
                p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                // Set up the discrete source
                std::shared_ptr<VesselBasedDiscreteSource<2> > p_vessel_source = VesselBasedDiscreteSource<2>::Create();
                QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                    GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                    Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                p_vessel_source->SetReferenceConcentration(vessel_oxygen_concentration);
                p_vessel_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_vessel_source->SetVesselPermeability(1.0*Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                p_oxygen_pde->AddDiscreteSource(p_vessel_source);

                // Set up the finite difference solver for oxygen (which handles everything)
                auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                p_oxygen_solver->SetPde(p_oxygen_pde);
                p_oxygen_solver->SetLabel("oxygen");
                p_oxygen_solver->SetGrid(p_grid);

                // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                double initial_haematocrit = 0.45;
                unsigned max_iter = 1000; 
                double tolerance2 = 1.e-6;  // Paraview only displays 6 decimal places anyway (too small and you succumb to floating point errors)
                std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                for(unsigned idx=0;idx<max_iter;idx++)
                {
                    // Run the solvers
                    impedance_calculator.Calculate();
                    flow_solver.SetUp();
                    flow_solver.Solve();
                    p_abstract_haematocrit_solver->Calculate();
                    p_viscosity_calculator->Calculate();

                    // Get the residual
                    double max_difference = 0.0;
                    double h_for_max = 0.0;
                    double prev_for_max = 0.0;
                    for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                    {
                        double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                        double difference = 0.0;
                        difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                        std::cout << "Current: " << current_haematocrit << "- Previous: " << previous_haematocrit[jdx] << " = Difference: " << difference << std::endl;
                        if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                        {
                            max_difference = difference;
                            h_for_max = current_haematocrit;
                            prev_for_max = previous_haematocrit[jdx];
                        }
                        previous_haematocrit[jdx] = current_haematocrit;
                    }
                    std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                    // Print the final or intermediary convergence results
                    if(max_difference<=tolerance2)  
                    {
                        std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                        break;
                    }
                    else
                    {
                        if(idx%1==0)
                        {
                            std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                            std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                            std::string output_file = p_output_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                            p_network->Write(output_file);
                        }
                    }

                    // If there is no convergence after all the iterations, print the error message.
                    if(idx==max_iter-1)
                    {
                        std::cout << "Problem encountered in Murray Network with h_solver = " << h_solver << " using alpha = " << n_alpha << std::endl;
                        EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                    }
                }
                
                // Run the simulation 
                SimulationTime::Instance()->SetStartTime(0.0);
                SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                p_microvessel_solver->SetVesselNetwork(p_network);
                p_microvessel_solver->SetOutputFileHandler(p_output_file_handler);
                p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                p_microvessel_solver->Run();

                // Print the average oxygenation
                std::vector<double> solution = p_oxygen_solver->GetSolution();
                double average_oxygen = 0.0;
                for(unsigned jdx=0;jdx<solution.size();jdx++)
                {
                    average_oxygen += solution[jdx];
                }
                average_oxygen /= double(solution.size());
                std::cout << "Average oxygen: " << average_oxygen << std::endl;
                
                // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                std::string output_file = p_output_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                p_network->Write(output_file);

                // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                ParameterCollection::Instance()->DumpToFile(p_output_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                ParameterCollection::Instance()->Destroy();
                BaseUnits::Instance()->Destroy();
                SimulationTime::Instance()->Destroy();

                // // Run the solver and write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                // solver.SetFileHandler(p_output_file_handler);
                // solver.SetWriteSolution(true);
                // solver.Solve();
                // p_network->Write(p_output_file_handler->GetOutputDirectoryFullPath() + "equilateral_network_results.vtp");

                // // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                // ParameterCollection::Instance()->DumpToFile(p_output_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                // ParameterCollection::Instance()->Destroy();
                // BaseUnits::Instance()->Destroy();
            }
        }
    }    

    // Make a 2D hexagonal network on a PDE grid with flow and H-splitting (with a neighbourhood of units)
    void xTestHexagonalNeighbourhoodWithFlow2D()
    {
        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Define the key parameters
        QLength large_vessel_radius = 10_um;
        double dimless_domain_size_x = 700.0;  // x-coordinate of output node
        double dimless_domain_size_y = 606.218 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
        unsigned dimless_vessel_length = 100.0;
        std::vector<std::vector<unsigned> > Order;
        std::shared_ptr<VesselNetwork<2> > p_network;
        std::vector<std::shared_ptr<Vessel<2> > > vessels;
        string line2;

        // Initialise the simulation space
        QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
        QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
        std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
        std::ofstream outfileMean;

        // Define the key parameters
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        double initial_haematocrit = 0.45;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";

        // Read the network from a file
        VesselNetworkGenerator<2> network_generator;
        std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix_Neighbourhood.txt");
        std::vector<std::vector<double> > rEdgesMatrix;
        string line;
        
        // Break down the rows in the file into column values
        while (std::getline(in, line)) 
        {
            rEdgesMatrix.push_back(std::vector<double>());
            std::stringstream split(line);
            double value;
            while (split >> value)
            {
                rEdgesMatrix.back().push_back(value);
            }
        }

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<6; h_solver++)
        {
            // Set file name
            std::ostringstream strs;
            strs << std::fixed << std::setprecision( 2 );
            strs << "TestDebuggingNetworks/HexagonalNeighbourhood/NoPruning/";        

            // std::ostringstream strs;
            if (h_solver==1)
            {
                strs << "ConstantHaematocrit";
            }
            else if (h_solver==2)
            {
                strs << "PriesHaematocrit";
            }
            else if (h_solver==3)
            {
                strs << "MemoryHaematocrit";
            }
            else if (h_solver==4)
            {
                strs << "FungHaematocrit";
            }
            else if (h_solver==5)
            {
                strs << "YangHaematocrit";
            }
            std::string str_directory_name = strs.str();
            auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

            // Generate the network
            p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);
            
            // Set inlet and outlet nodes
            auto p_segment = p_network->GetVesselSegments()[0];
            vessels = p_network->GetVessels();
            for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
            {
                if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
                {
                    if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)
                    {
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
                    }
                    if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
                    {
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
                        (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
                    }
                }
                if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)   	
                {
                    if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)
                    {
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
                    }
                    if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
                    {
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
                        (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
                    }
                }
            }
            p_segment = p_network->GetVesselSegments()[0];
            p_segment->SetRadius(large_vessel_radius);
            VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   

            // Set the haematocrit solver
            std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
            if (h_solver==1)
            {
                std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;       
            }
            else if (h_solver==2)
            {
                std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;                    
            }
            else if (h_solver==3)
            {
                std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;      
            }
            else if (h_solver==4)
            {
                std::cout << "Now using FungHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;      
            }
            else if (h_solver==5)
            {
                std::cout << "Now using YangHaematocritSolver..." << std::endl;
                auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                p_haematocrit_solver->SetVesselNetwork(p_network);
                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                p_abstract_haematocrit_solver = p_haematocrit_solver;      
            }

            // Set up the viscosity solver
            auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
            p_viscosity_calculator->SetPlasmaViscosity(viscosity);
            p_viscosity_calculator->SetVesselNetwork(p_network);
            p_viscosity_calculator->Calculate();
            
            // Set up the impedance solver
            auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
            p_impedance_calculator->SetVesselNetwork(p_network);
            p_impedance_calculator->Calculate();

            // Set up the flow solver 
            FlowSolver<2> flow_solver;
            flow_solver.SetVesselNetwork(p_network);
            // flow_solver.SetUp();
            flow_solver.SetUseDirectSolver(true);
            // flow_solver.Solve();

            // Set up the grid for the finite difference solver
            auto p_grid = RegularGrid<2>::Create();
            p_grid->SetSpacing(irregular_grid_spacing);
            c_vector<unsigned, 3> dimensions;
            dimensions[0] = unsigned((domain_side_length_x)/(irregular_grid_spacing))+1; // num x
            dimensions[1] = unsigned((domain_side_length_y)/(irregular_grid_spacing))+1; // num_y
            dimensions[2] = 1;
            p_grid->SetDimensions(dimensions);

            // Choose the PDE
            auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
            
            // Set the diffusivity and decay terms
            p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
            p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

            // Set up the discrete source
            auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
            QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                    GenericParameters::mpGasConcentrationAtStp->GetValue("User");
            QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                    Owen11Parameters::mpReferencePartialPressure->GetValue("User");
            p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
            p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
            p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
            p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

            // Set up the finite difference solver for oxygen (which handles everything)
            auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
            p_oxygen_solver->SetPde(p_oxygen_pde);
            p_oxygen_solver->SetLabel("oxygen");
            p_oxygen_solver->SetGrid(p_grid);
            // solver.SetFileName("oxygen_solution_0");

            // Set up flag for broken solver
            unsigned broken_solver = 0;

            // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
            unsigned max_iter = 1000;  // 1000 
            double tolerance2 = 1.e-10;
            std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
            std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
            for(unsigned idx=0;idx<max_iter;idx++)
            {
                // Run the solvers
                p_impedance_calculator->Calculate();
                flow_solver.SetUp();
                flow_solver.Solve();
                p_abstract_haematocrit_solver->Calculate();
                p_viscosity_calculator->Calculate();

                // Get the residual
                double max_difference = 0.0;
                double h_for_max = 0.0;
                double prev_for_max = 0.0;
                for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                {
                    double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                    double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                    if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                    {
                        max_difference = difference;
                        h_for_max = current_haematocrit;
                        prev_for_max = previous_haematocrit[jdx];
                    }
                    previous_haematocrit[jdx] = current_haematocrit;
                }
                std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                // Print the final or intermediary convergence results
                if(max_difference<=tolerance2)  
                {
                    std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;                    
                    broken_solver = 0;
                    break;
                }
                else
                {
                    if(idx%1==0)
                    {
                        std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                        std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                        p_network->Write(output_file);
                    }
                }

                // If there is no convergence after all the iterations, print the error message.
                if(idx==max_iter-1)
                {
                    std::cout << "Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << std::endl;
                    error_log << "\n Problem encountered in hexagonal network with h_solver = " << to_string(h_solver); 
                    broken_solver = 1;
                    break;
                    // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                }
            }

            // If solver doesn't converge, move on to next one
            if (broken_solver == 1)
            {
                continue;
            }

            // Run the simulation 
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
            auto p_microvessel_solver = MicrovesselSolver<2>::Create();
            p_microvessel_solver->SetVesselNetwork(p_network);
            p_microvessel_solver->SetOutputFileHandler(p_file_handler);
            p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
            p_microvessel_solver->Run();

            // Print the average oxygenation
            std::vector<double> solution = p_oxygen_solver->GetSolution();
            double average_oxygen = 0.0;
            for(unsigned jdx=0;jdx<solution.size();jdx++)
            {
                average_oxygen += solution[jdx];
            }
            average_oxygen /= double(solution.size());
            std::cout << "Average oxygen: " << average_oxygen << std::endl;
            
            // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
            std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
            p_network->Write(output_file);

            // Dump our parameter collection to an xml file and, importantly, clear it for the next test
            ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
            ParameterCollection::Instance()->Destroy();
            BaseUnits::Instance()->Destroy();
            SimulationTime::Instance()->Destroy();
        }
         
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a parallel network on a PDE grid that acts as a Dirichlet BC in 2D 
    void xTestParallelNetworkLineSource2D()
    {
        double alpha = 1.0;  // alpha determines the relative radius of the left vessel 
    
        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Set key vessel parameters
        QLength vessel_length(100.0*unit::microns);
        QLength vessel_height = (pow(2,0.5)*vessel_length)*0.5;

        // Set up the domain parameters
        QLength domain_x = vessel_height+vessel_height+vessel_length+vessel_length+0.0001_um;  // this should extend to the x-position of outlet node
        QLength domain_y = domain_x;  // should be the same as domain_x to make square domain
        QLength mid_domain_y = domain_y*0.5;

        // Set nodes based on an equilateral network
        std::shared_ptr<VesselNode<2> > p_node_1 = VesselNode<2>::Create(0.0_um, mid_domain_y+50.0_um);
        std::shared_ptr<VesselNode<2> > p_node_2 = VesselNode<2>::Create(domain_x, mid_domain_y+50.0_um);

        std::shared_ptr<VesselNode<2> > p_node_3 = VesselNode<2>::Create(0.0_um, mid_domain_y);
        std::shared_ptr<VesselNode<2> > p_node_4 = VesselNode<2>::Create(domain_x, mid_domain_y);

        std::shared_ptr<VesselNode<2> > p_node_5 = VesselNode<2>::Create(0.0_um, mid_domain_y-50.0_um);
        std::shared_ptr<VesselNode<2> > p_node_6 = VesselNode<2>::Create(domain_x, mid_domain_y-50.0_um);

        std::shared_ptr<VesselNode<2> > p_node_7 = VesselNode<2>::Create(0.0_um, mid_domain_y-60.0_um);
        std::shared_ptr<VesselNode<2> > p_node_8 = VesselNode<2>::Create(domain_x, mid_domain_y-60.0_um);

        std::shared_ptr<VesselNode<2> > p_node_9 = VesselNode<2>::Create(0.0_um, mid_domain_y-70.0_um);
        std::shared_ptr<VesselNode<2> > p_node_10 = VesselNode<2>::Create(domain_x, mid_domain_y-70.0_um);

        // std::shared_ptr<VesselNode<2> > p_node_7 = VesselNode<2>::Create(0.0_um, mid_domain_y);  // add a node at a distance of 100 picometres from the last node to work around glitch of disconnected flow when setting p_node_6 as the output node

        // Make segments 
        std::shared_ptr<VesselSegment<2> > p_segment_1 = VesselSegment<2>::Create(p_node_1, p_node_2);
        std::shared_ptr<VesselSegment<2> > p_segment_2 = VesselSegment<2>::Create(p_node_3, p_node_4);
        std::shared_ptr<VesselSegment<2> > p_segment_3 = VesselSegment<2>::Create(p_node_5, p_node_6);
        std::shared_ptr<VesselSegment<2> > p_segment_4 = VesselSegment<2>::Create(p_node_7, p_node_8);
        std::shared_ptr<VesselSegment<2> > p_segment_5 = VesselSegment<2>::Create(p_node_9, p_node_10);

        // Make vessels
        std::shared_ptr<Vessel<2> > p_vessel_1 = Vessel<2>::Create(p_segment_1);
        std::shared_ptr<Vessel<2> > p_vessel_2 = Vessel<2>::Create(p_segment_2);
        std::shared_ptr<Vessel<2> > p_vessel_3 = Vessel<2>::Create(p_segment_3);
        std::shared_ptr<Vessel<2> > p_vessel_4 = Vessel<2>::Create(p_segment_4);
        std::shared_ptr<Vessel<2> > p_vessel_5 = Vessel<2>::Create(p_segment_5);

        // Add the vessels to a vessel network
        std::shared_ptr<VesselNetwork<2> > p_network = VesselNetwork<2>::Create();
        p_network->AddVessel(p_vessel_1);
        p_network->AddVessel(p_vessel_2);
        p_network->AddVessel(p_vessel_3);
        p_network->AddVessel(p_vessel_4);
        p_network->AddVessel(p_vessel_5);

        // No pruning
        std::ostringstream strs;
        strs << std::fixed << std::setprecision( 1 );
        strs << "TestEquilateralNetwork/DirichletHaematocrit/Alpha" << alpha << "/NoPruning/";
        std::string str_directory_name = strs.str();
        auto p_output_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);         
        // p_network->RemoveVessel(p_vessel_7,true);

        // Set up the grid for the finite difference solver
        auto p_grid = RegularGrid<2>::Create();
        p_grid->SetSpacing(equilateral_grid_spacing);
        c_vector<unsigned, 3> dimensions;
        dimensions[0] = unsigned((domain_x)/(equilateral_grid_spacing))+1; // num x
        dimensions[1] = unsigned((domain_y)/(equilateral_grid_spacing))+1; // num_y
        dimensions[2] = 1;
        p_grid->SetDimensions(dimensions);

        // Choose the PDE
        std::shared_ptr<DiscreteContinuumLinearEllipticPde<2> > p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();

        // Set the diffusivity and decay terms
        p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
        p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

        // Set the O2 concentration value
        QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") * GenericParameters::mpGasConcentrationAtStp->GetValue("User");
        QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp * Owen11Parameters::mpReferencePartialPressure->GetValue("User");
        std::cout << vessel_oxygen_concentration << std::endl; 

        // Set the boundary condition to be the network acting as a Dirichlet source
        std::shared_ptr<DiscreteContinuumBoundaryCondition<2> > p_vessel_boundary_condition = DiscreteContinuumBoundaryCondition<2>::Create();
        p_vessel_boundary_condition->SetValue(vessel_oxygen_concentration);  // set the O2 concentration at the source
        p_vessel_boundary_condition->SetType(BoundaryConditionType::VESSEL_LINE);
        p_vessel_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);

        // Set up the finite difference solver for oxygen (which handles everything)
        SimpleLinearEllipticFiniteDifferenceSolver<2> solver;
        solver.SetGrid(p_grid);
        solver.SetPde(p_oxygen_pde);
        solver.AddBoundaryCondition(p_vessel_boundary_condition);
        solver.SetVesselNetwork(p_network);
        solver.SetLabel("oxygen");
        solver.SetFileName("oxygen_solution_0");

        // Run the solver and write the output file (visualise with Paraview; set Filters->Alphabetical->Tube)
        solver.SetFileHandler(p_output_file_handler);
        solver.SetWriteSolution(true);
        solver.Solve();
        p_network->Write(p_output_file_handler->GetOutputDirectoryFullPath() + "equilateral_network_results.vtp");

        // Dump our parameter collection to an xml file and, importantly, clear it for the next test
        ParameterCollection::Instance()->DumpToFile(p_output_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
        ParameterCollection::Instance()->Destroy();
        BaseUnits::Instance()->Destroy();
        
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Archival Forking Networks
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Choose the number of bifurcating generations (inlet and outlet vessels don't count)
    unsigned order = 7;

    // Make a multi-generation forking network on a PDE grid as a Dirichlet BC in 2D 
    void xTestDichotomousNetworkLineSource2D()
    {
        // Run the simulation with different heterogeneities
        for (unsigned n_alpha=0; n_alpha<=max_alpha; n_alpha++)
        {
            double alpha = 1.0+(double)n_alpha*0.1;  // alpha determines the relative radius of the left vessel 
            
            // Set up the reference length for the simulation
            QLength reference_length(1.0_um);
            BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

            // Set key vessel parameters
            double dimless_length = 1.0;  

            // ???
            for(unsigned i_aux=1; i_aux<order+1; i_aux++)
            {
                dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
            }

            // Set input radius
            QLength input_radius(10.0*GenericParameters::mpCapillaryRadius->GetValue());  // *5_um
            // VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);

            // Generate the networks
            VesselNetworkGenerator<2> network_generator;
            double lambda;  // lambda = length/diameter
            double twicelambda;  // used as an input parameter
            // for (unsigned k_aux=1; k_aux<5; k_aux++)  // generate the network for various lambdas
            for (unsigned k_aux=1; k_aux<2; k_aux++)  // generate the network for lambda = 4
            {
                lambda = 2.0+double(k_aux)*2.0;
                twicelambda = 2.0*lambda;

                // Height of of first-order vessels
                QLength main_vert_length = 0.9*twicelambda*input_radius*pow(2.0,-1.0/3.0);
                
                // Height of the domain
                QLength domain_side_length_y = 4.0*main_vert_length;

                // Length of the domain
                QLength domain_side_length_x = dimless_length*2.0*twicelambda*input_radius;

                // Generate the name of the directory based on the value of lambda
                std::ostringstream strs;
                strs << std::fixed << std::setprecision( 1 );
                strs << "TestDichotomousNetwork/DirichletHaematocrit/Alpha" << alpha << "/NoPruning/LambdaEquals" << lambda;
                std::string str_directory_name = strs.str();
                auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                // Generate the network
                std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapy(order, main_vert_length, input_radius, alpha, twicelambda);

                // Set up the grid for the finite difference solver
                auto p_grid = RegularGrid<2>::Create();
                // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
                // QLength grid_spacing(1.0*1_um);  // the simulation time gets quite long if you reduce the resolution further
                // QLength grid_spacing = 10_um;
                p_grid->SetSpacing(grid_spacing);
                c_vector<unsigned, 3> dimensions;
                dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                dimensions[2] = 1;
                p_grid->SetDimensions(dimensions);
                // std::shared_ptr<Part<2> > p_domain = Part<2>::Create();
                // p_domain->AddRectangle(domain_x, domain_y);
                // std::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
                // p_grid->GenerateFromPart(p_domain, grid_spacing);  // set domain and grid spacing

                // Choose the PDE
                auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                
                // Set the diffusivity and decay terms
                p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                // Set the O2 concentration value
                QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") * GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp * Owen11Parameters::mpReferencePartialPressure->GetValue("User");

                // Set the boundary condition to be the network acting as a Dirichlet source
                std::shared_ptr<DiscreteContinuumBoundaryCondition<2> > p_vessel_boundary_condition = DiscreteContinuumBoundaryCondition<2>::Create();
                p_vessel_boundary_condition->SetValue(vessel_oxygen_concentration);  // set the O2 concentration at the source
                p_vessel_boundary_condition->SetType(BoundaryConditionType::VESSEL_LINE);
                p_vessel_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);
                    
                // Set up the finite difference solver for oxygen (which handles everything)
                SimpleLinearEllipticFiniteDifferenceSolver<2> solver;
                solver.SetGrid(p_grid);
                solver.SetPde(p_oxygen_pde);
                solver.AddBoundaryCondition(p_vessel_boundary_condition);
                solver.SetVesselNetwork(p_network);
                solver.SetLabel("oxygen");
                solver.SetFileName("oxygen_solution_0");

                // Run the solver and write the output file (visualise with Paraview; set Filters->Alphabetical->Tube)
                solver.SetFileHandler(p_file_handler);
                solver.SetWriteSolution(true);
                solver.Solve();
                p_network->Write(p_file_handler->GetOutputDirectoryFullPath() + "FinalHaematocrit.vtp");
        
                // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                ParameterCollection::Instance()->Destroy();
                BaseUnits::Instance()->Destroy();
            }
        }
    }

    // Make a multi-generation forking network on a PDE grid with different haematocrit splitting rules
    void xTestDichotomousNetworkWithFlow2D()
    {
        // Run the simulation with different heterogeneities
        for (unsigned n_alpha=0; n_alpha<=0; n_alpha++)
        {
            double alpha = 1.0+(double)n_alpha*0.1;  // alpha determines the relative radius of the left vessel 
            
            // Run the simulation with different solvers of interest
            for (unsigned h_solver=0; h_solver<1; h_solver++)
            {      
                // Set up the reference length for the simulation
                QLength reference_length(1.0_um);
                BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

                // Set key vessel parameters
                // QLength vessel_length(100.0_um);
                // QLength vessel_height = (pow(2,0.5)*vessel_length)*0.5;
                double dimless_length = 1.0;  

                // ???
                for(unsigned i_aux=1; i_aux<order+1; i_aux++)
                {
                    dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
                }

                // Set input radius
                QLength input_radius(5.0 *GenericParameters::mpCapillaryRadius->GetValue());  // = 50_um
                // VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);

                // Set the viscosity
                QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
                // QDynamicViscosity viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();  // = 0.0012

                // Set the inlet and initial haematocrit
                // p_haematocrit_calculator->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));  // = 0.45
                // double inlet_haematocrit = 0.45;
                double initial_haematocrit = 0.45;

                // Generate the networks
                VesselNetworkGenerator<2> network_generator;
                double lambda;  // lambda = length/diameter
                double twicelambda;  // used as an input parameter
                // for (unsigned k_aux=1; k_aux<5; k_aux++)  // generate the network for various lambdas
                for (unsigned k_aux=1; k_aux<2; k_aux++)  // generate the network for lambda = 4
                {
                    lambda = 2.0+double(k_aux)*2.0;
                    twicelambda = 2.0*lambda;

                    // Height of of first-order vessels
                    QLength main_vert_length = 0.9*twicelambda*input_radius*pow(2.0,-1.0/3.0);
                    
                    // Height of the domain
                    // QLength domain_side_length_y = 4.0*main_vert_length;

                    // Length of the domain
                    QLength domain_side_length_x = dimless_length*2.0*twicelambda*input_radius;

                    // Set file name based on haematocrit solver
                    std::ostringstream strs;
                    strs << std::fixed << std::setprecision( 1 );
                    if (h_solver==1)
                    {
                        strs << "TestDichotomousNetwork/ConstantHaematocrit/Alpha" << alpha << "/NoPruning/GridSpacingEquals" << grid_spacing;
                    }
                    else if (h_solver==2)
                    {
                        strs << "TestDichotomousNetwork/PriesHaematocrit/Alpha" << alpha << "/NoPruning/GridSpacingEquals" << grid_spacing;
                    }
                    else if (h_solver==3)
                    {
                        strs << "TestDichotomousNetwork/MemoryHaematocrit/Alpha" << alpha << "/NoPruning/GridSpacingEquals" << grid_spacing;
                    }
                    else if (h_solver==4)
                    {
                        strs << "TestDichotomousNetwork/FungHaematocrit/Alpha" << alpha << "/NoPruning/GridSpacingEquals" << grid_spacing;
                    }
                    else if (h_solver==5)
                    {
                        strs << "TestDichotomousNetwork/YangHaematocrit/Alpha" << alpha << "/NoPruning/GridSpacingEquals" << grid_spacing;
                    }
                    std::string str_directory_name = strs.str();
                    auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                    // Generate the network
                    std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapy(order, main_vert_length, input_radius, alpha, twicelambda);

                    // Identify input and output nodes and assign them properties
                    VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(0.0_um,2.0*main_vert_length));
                    VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(domain_side_length_x, 2.0*main_vert_length));
                    p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
                    // p_network->GetNode(0)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));  // = 3333.06
                    p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
                    p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
                    // p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));  // = 1999.84
                    p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);

                    // Set up the grid for the finite difference solver
                    // auto p_grid = RegularGrid<2>::Create();
                    // // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
                    // // QLength grid_spacing(1.0*1_um);  // the simulation time gets quite long if you reduce the resolution further
                    // QLength grid_spacing = 10_um;
                    // p_grid->SetSpacing(grid_spacing);
                    // c_vector<unsigned, 3> dimensions;
                    // dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                    // dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                    // dimensions[2] = 1;
                    // p_grid->SetDimensions(dimensions);
                    // // std::shared_ptr<Part<2> > p_domain = Part<2>::Create();
                    // // p_domain->AddRectangle(domain_x, domain_y);
                    // // std::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
                    // // p_grid->GenerateFromPart(p_domain, grid_spacing);  // set domain and grid spacing

                    // // Choose the PDE
                    // auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                    
                    // // Set the diffusivity and decay terms
                    // p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                    // p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                    // // Set up the discrete source
                    // auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                    // QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                    //         GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                    // QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                    //         Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                    // p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                    // p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                    // p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    // p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                    // // Set up the finite difference solver for oxygen (which handles everything)
                    // auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                    // p_oxygen_solver->SetPde(p_oxygen_pde);
                    // p_oxygen_solver->SetLabel("oxygen");
                    // p_oxygen_solver->SetGrid(p_grid);

                    // Set the haematocrit solver
                    std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                    if (h_solver==1)
                    {
                        std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;       
                    }
                    else if (h_solver==2)
                    {
                        std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                    
                    }
                    else if (h_solver==3)
                    {
                        std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;      
                    }
                    else if (h_solver==4)
                    {
                        std::cout << "Now using FungHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;      
                    }
                    else if (h_solver==5)
                    {
                        std::cout << "Now using YangHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;      
                    }

                    // Set up the viscosity calculator
                    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                    p_viscosity_calculator->SetVesselNetwork(p_network);
                    p_viscosity_calculator->Calculate();
                    
                    // Set up the impedance calculator
                    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                    p_impedance_calculator->SetVesselNetwork(p_network);
                    p_impedance_calculator->Calculate();

                    // Set up the flow solver
                    FlowSolver<2> flow_solver;
                    flow_solver.SetVesselNetwork(p_network);
                    flow_solver.SetUp();

                    // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance. The loop runs the solvers, calculates the max. difference in the network between the previous H of a vessel and the new H in the vessel, and repeats this until all vessels converge on a value within a certain tolerance.
                    // unsigned max_iter = 1000;  // simulation fails if it doesn't converge in these many iterations
                    // double tolerance2 = 1.e-10;
                    std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                    std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                    // for(unsigned iteration=0;iteration<max_iter;iteration++)
                    // {
                        // Run the solvers
                        p_impedance_calculator->Calculate();
                        flow_solver.SetUp();
                        flow_solver.Solve();
                        p_abstract_haematocrit_solver->Calculate();
                        p_viscosity_calculator->Calculate();

                        // Check for convergence 
                        // double max_difference = 0.0;
                        // double h_for_max = 0.0;
                        // double prev_for_max = 0.0;
                        for(unsigned segment_index=0;segment_index<segments.size();segment_index++)  // for all the segments in the network
                        {
                            // // Set segments with no flow to have no haematocrit
                            if (fabs(segments[segment_index]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                            {
                                segments[segment_index]->GetFlowProperties()->SetHaematocrit(0.0);
                            }

                            // double current_haematocrit = segments[segment_index]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                            // double difference = std::abs(current_haematocrit - previous_haematocrit[segment_index]);  // difference in new H in segment and old H
                            
                            // // Log the max. difference calculated
                            // if(difference>max_difference)
                            // {
                            //     max_difference = difference;
                            //     h_for_max = current_haematocrit;
                            //     prev_for_max = previous_haematocrit[segment_index];
                            // }
                            // previous_haematocrit[segment_index] = current_haematocrit;
                        }

                    //     // Print the max. difference of the iteration
                    //     std::cout << "The maximum difference calculated in iteration " << iteration << " is " << h_for_max << " - " << prev_for_max << " = " << max_difference << std::endl;

                    //     if(max_difference<=tolerance2)  
                    //     {
                    //         std::cout << "Converged after: " << iteration << " iterations." <<  std::endl;
                    //         break;
                    //     }
                    //     else
                    //     {
                    //         if(iteration%1==0)
                    //         {
                    //             // std::cout << "Max Difference at iter: " << iteration << " is " << max_difference << std::endl;
                    //             std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(iteration) + ".vtp";
                    //             std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                    //             p_network->Write(output_file);
                    //         }
                    //     }
                    //     if(iteration==max_iter-1)  // if there is no convergence after all the iterations, print the error message
                    //     {
                    //         std::cout << "Problem encountered in Dichotomous Network with h_solver = " << h_solver << " using alpha = " << n_alpha << " and lambda = " << lambda << std::endl;
                    //         EXCEPTION("Did not converge after " + std::to_string(iteration) + " iterations.");
                    //     }
                    // }
                
                    // // Run the simulation 
                    // SimulationTime::Instance()->SetStartTime(0.0);
                    // SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                    // auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                    // p_microvessel_solver->SetVesselNetwork(p_network);
                    // p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                    // p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                    // p_microvessel_solver->Run();

                    // // Print the average oxygenation
                    // std::vector<double> solution = p_oxygen_solver->GetSolution();
                    // double average_oxygen = 0.0;
                    // for(unsigned jdx=0;jdx<solution.size();jdx++)
                    // {
                    //     average_oxygen += solution[jdx];
                    // }
                    // average_oxygen /= double(solution.size());
                    // std::cout << "Average oxygen: " << average_oxygen << std::endl;
                    
                    // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                    std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                    p_network->Write(output_file);

                    // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                    ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                    ParameterCollection::Instance()->Destroy();
                    BaseUnits::Instance()->Destroy();
                    SimulationTime::Instance()->Destroy();
                }
            }
        }
    }

    // Make a multi-generation forking network on a PDE grid with different haematocrit splitting rules and radius threshold pruning
    void xTestDichotomousNetworkWithPruningAndFlow2D()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";

        // Set network in filename
        std::string network_name = "TestDichotomousNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("forking_radius_threshold_perfusion_quotients.txt");
        outfile.close();

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=1; h_solver++)
        {   
            // Generate the network for various lambdas
            for (unsigned k_aux=1; k_aux<2; k_aux++)  // generate the network for lambda = 4
            {
            // // Run the simulation with different heterogeneities
            // for (unsigned n_alpha=0; n_alpha<=max_alpha; n_alpha++)
            // { 
                // Set key vessel parameters
                double dimless_length = 1.0;  

                // ??? Non-dimensionalising the length
                for(unsigned i_aux=1; i_aux<order+1; i_aux++)
                {
                    dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
                }
   
                // Set up the reference length for the simulation
                QLength reference_length(1.0_um);
                BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

                // Set input radius
                QLength input_radius(5.0 *GenericParameters::mpCapillaryRadius->GetValue());  // = 50_um
                // VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);

                // Set the viscosity
                QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
                // QDynamicViscosity viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();  // = 0.0012

                // Set the inlet and initial haematocrit
                double initial_haematocrit = 0.45;

                // Get pre-RT Perfusion Quotient
		        QFlowRate threshold = 1.e-12*unit::metre_cubed_per_second;  // set flow threshold for what counts as a perfused vessel

                // Generate the networks
                VesselNetworkGenerator<2> network_generator;

                // Set lambda
                double lambda;  // lambda = length/diameter
                double twicelambda;  // used as an input parameter
                lambda = 2.0+double(k_aux)*2.0;
                twicelambda = 2.0*lambda;

                // Run the simulation with different heterogeneities
                for (unsigned n_alpha=0; n_alpha<=max_alpha; n_alpha++)
                { 
                    // Set up flag for broken solver
                    unsigned broken_solver = 0;

                    // Set alpha
                    double alpha = 1.0+(double)n_alpha*0.1;  // alpha determines the relative radius of the left vessel 

                    // Height of of first-order vessels
                    QLength main_vert_length = 0.9*twicelambda*input_radius*pow(2.0,-1.0/3.0);
                    
                    // Height of the domain
                    QLength domain_side_length_y = 4.0*main_vert_length;

                    // Length of the domain
                    QLength domain_side_length_x = dimless_length*2.0*twicelambda*input_radius;

                    // Generate the network
                    std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapy(order, main_vert_length, input_radius, alpha, twicelambda);

                    // Identify input and output nodes and assign them properties
                    VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(0.0_um,2.0*main_vert_length));
                    VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(domain_side_length_x, 2.0*main_vert_length));
                    p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
                    // p_network->GetNode(0)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));  // = 3333.06
                    p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
                    p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
                    // p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));  // = 1999.84
                    p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);

                    // Set the h-solver
                    std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                    std::string solver_name;
                    if (h_solver==1)
                    {
                        solver_name = "ConstantHaematocrit/"; 
                        std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==2)
                    {
                        solver_name = "PriesHaematocrit/"; 
                        std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                
                    }
                    else if (h_solver==3)
                    {
                        solver_name = "MemoryHaematocrit/"; 
                        std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==4)
                    {
                        solver_name = "FungHaematocrit/"; 
                        std::cout << "Now using FungHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==5)
                    {
                        solver_name = "YangHaematocrit/"; 
                        std::cout << "Now using YangHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();   
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                                    
                    // Store the lambda values for file name
                    std::stringstream lambda_stream;
                    lambda_stream << std::fixed << std::setprecision(0) << lambda;
                    std::string lambda_string = lambda_stream.str();
                    std::string lambda_directory = network_name + solver_name + "/Lambda" + lambda_string; 

                    // Set up the grid for the finite difference solver
                    auto p_grid = RegularGrid<2>::Create();
                    // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
                    // QLength grid_spacing(1.0*1_um);  // the simulation time gets quite long if you reduce the resolution further
                    // QLength grid_spacing = 10_um;
                    p_grid->SetSpacing(grid_spacing);
                    c_vector<unsigned, 3> dimensions;
                    dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                    dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                    dimensions[2] = 1;
                    p_grid->SetDimensions(dimensions);

                    // Choose the PDE
                    auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                    
                    // Set the diffusivity and decay terms
                    p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                    p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                    // Set up the discrete source
                    auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                    QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                            GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                    QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                            Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                    p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                    p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                    p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                    // Set up the finite difference solver for oxygen (which handles everything)
                    auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                    p_oxygen_solver->SetPde(p_oxygen_pde);
                    p_oxygen_solver->SetLabel("oxygen");
                    p_oxygen_solver->SetGrid(p_grid);

                    // Set up the viscosity calculator
                    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                    p_viscosity_calculator->SetVesselNetwork(p_network);
                    p_viscosity_calculator->Calculate();
                    
                    // Set up the impedance calculator
                    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                    p_impedance_calculator->SetVesselNetwork(p_network);
                    p_impedance_calculator->Calculate();

                    // Set up the flow solver
                    FlowSolver<2> flow_solver;
                    flow_solver.SetVesselNetwork(p_network);
                    // flow_solver.SetUp();
                    // // flow_solver.SetUseDirectSolver(false);
                    // flow_solver.Solve();

                    // Set up pruning parameters
                    QLength max_radius_to_kill = 31_um;  // set threshold up to which vessels should be pruned
                    QLength radius_step = 1_um;
                    QLength current_radius_threshold = 0_um;

                    // Prunes vessels from thinnest up to simulate RT
                    while (current_radius_threshold <= max_radius_to_kill)
                    {         
                        // Display status message
                        double current_radius_threshold_um = current_radius_threshold*1000000;  // in micrometres
                        std::cout << "Now pruning up to vessels with radius = " << current_radius_threshold_um << " um" << std::endl;

                        // Remove vessels under radius threshold
                        p_network->RemoveThinVessels(current_radius_threshold, false);  

                        // Set filename
                        std::stringstream alpha_stream;
                        std::stringstream threshold_stream;
                        alpha_stream << std::fixed << std::setprecision(2) << alpha;
                        threshold_stream << std::fixed << std::setprecision(0) << current_radius_threshold_um;
                        std::string alpha_string = alpha_stream.str();
                        std::string threshold_string = threshold_stream.str();
                        std::string file_name = lambda_directory + "/Alpha" + alpha_string + "/RadiusThreshold" + threshold_string;
                        std::string str_directory_name = file_name;
                        auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                        // outfile << alpha << " " << current_radius_threshold << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 

                        // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance. The loop runs the solvers, calculates the max. difference in the network between the previous H of a vessel and the new H in the vessel, and repeats this until all vessels converge on a value within a certain tolerance.
                        unsigned max_iter = 1000;  // simulation fails if it doesn't converge in these many iterations
                        double tolerance2 = 0.001;
                        std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                        std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                        for (unsigned idx=0;idx<max_iter;idx++)
                        {
                            // Run the solvers
                            p_abstract_haematocrit_solver->Calculate();
                            p_viscosity_calculator->Calculate();
                            p_impedance_calculator->Calculate();
                            flow_solver.SetUp();
                            flow_solver.Solve();

                            // Check for convergence 
                            double max_difference = 0.0;
                            double h_for_max = 0.0;
                            double prev_for_max = 0.0;
                            for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)  // for all the segments 
                            {
                                // Set segments with no flow to have no haematocrit
                                if (fabs(segments[segment_index]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                {
                                    segments[segment_index]->GetFlowProperties()->SetHaematocrit(0.0);
                                }
                                
                                // Check for convergence
                                double current_haematocrit = segments[segment_index]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                                double difference = std::abs(current_haematocrit - previous_haematocrit[segment_index]);  // difference in new H in segment and old H
                                if(difference>max_difference)  // if the difference is greater than previous max.
                                {
                                    max_difference = difference;
                                    h_for_max = current_haematocrit;
                                    prev_for_max = previous_haematocrit[segment_index];
                                }
                                previous_haematocrit[segment_index] = current_haematocrit;
                            }
                            std::cout << "Segment H. at max difference: " << h_for_max << ", Prev. segment H at max difference:" << prev_for_max << std::endl;
                            if(max_difference<=tolerance2)  
                            {
                                std::cout << "Converged after: " << idx << " iterations." <<  std::endl;
                                broken_solver = 0;
                                break;
                            }
                            else
                            {
                                if(idx%1==0)
                                {
                                    std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                                    std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                                    std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                                    p_network->Write(output_file);
                                }
                            }
                            if(idx==max_iter-1)  // if there is no convergence after all the iterations, print the error message
                            {
                                std::cout << "Problem encountered in Dichotomous Network with h_solver = " << h_solver << " using n_alpha = " << n_alpha << " and lambda = " << lambda << " and radius threshold = " << current_radius_threshold << std::endl;
                                // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                                error_log << "\n Problem encountered in Dichotomous Network with h_solver = " << h_solver << " using n_alpha = " << n_alpha << " and lambda = " << lambda << " and radius threshold = " << current_radius_threshold; 
                                broken_solver = 1;
                                break;
                            }
                        }

                        // If solver doesn't converge, move on to next one
                        if (broken_solver == 1)
                        {
                            continue;
                        }

                        // Run the simulation 
                        SimulationTime::Instance()->SetStartTime(0.0);
                        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                        auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                        p_microvessel_solver->SetVesselNetwork(p_network);
                        p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                        p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                        p_microvessel_solver->Run();

                        // Print the average oxygenation
                        std::vector<double> solution = p_oxygen_solver->GetSolution();
                        double average_oxygen = 0.0;
                        for(unsigned jdx=0;jdx<solution.size();jdx++)
                        {
                            average_oxygen += solution[jdx];
                        }
                        average_oxygen /= double(solution.size());
                        std::cout << "Average oxygen: " << average_oxygen << std::endl;

                        // Write the PQs
                        outfile.open("forking_radius_threshold_perfusion_quotients.txt", std::ios_base::app);
                        outfile << network_name << " " << solver_name << " " << lambda_string << " " << alpha_string << " " << threshold_string << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                        outfile.close();
                        
                        // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                        p_network->Write(output_file);

                        // // ??? Assuming radius doesn't exceed threshold, write the network to the output file
                        // if(current_radius_threshold>12.9_um && current_radius_threshold<13.1_um)
                        // {
                        //     p_network->Write(output_file_final_RT);
                        // }
                        current_radius_threshold = current_radius_threshold + radius_step;  // move onto next radius threshold for next step

                        // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                        ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                        ParameterCollection::Instance()->Destroy();
                        BaseUnits::Instance()->Destroy();
                        SimulationTime::Instance()->Destroy();
                    }
                }
            }
        }
        // Close the output file

        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a multi-generation forking network on a PDE grid with different h-splitting rules and stochastic pruning
    void xTestDichotomousNetworkWithStochasticPruningAndFlow2D()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";

        // Set network in filename
        std::string network_name = "TestDichotomousNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("forking_stochastic_pruning_perfusion_quotients.txt", std::ios_base::app);
        outfile.close();

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=1; h_solver++)
        {   
            // Generate the network for various lambdas
            for (unsigned k_aux=1; k_aux<2; k_aux++)  // generate the network for lambda = 4
            {
                // Set key vessel parameters
                double dimless_length = 1.0;  

                // ??? Non-dimensionalising the length
                for(unsigned i_aux=1; i_aux<order+1; i_aux++)
                {
                    dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
                }
   
                // Set up the reference length for the simulation
                QLength reference_length(1.0_um);
                BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

                // Set input radius
                QLength input_radius(5.0 *GenericParameters::mpCapillaryRadius->GetValue());  // = 50_um
                // VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);

                // Set the viscosity
                QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
                // QDynamicViscosity viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();  // = 0.0012

                // Set the inlet and initial haematocrit
                double initial_haematocrit = 0.45;

                // Set threshold for perfusion quotient
		        QFlowRate threshold = 1.e-12*unit::metre_cubed_per_second;  // set flow threshold for what counts as a perfused vessel

                // Generate the networks
                VesselNetworkGenerator<2> network_generator;

                // Set lambda
                double lambda;  // lambda = length/diameter
                double twicelambda;  // used as an input parameter
                lambda = 2.0+double(k_aux)*2.0;
                twicelambda = 2.0*lambda;

                // Run the simulation with different heterogeneities
                for (unsigned n_alpha=0; n_alpha<=4; n_alpha++)
                { 
                    // Set up flag for broken solver
                    unsigned broken_solver = 0;

                    // Set alpha
                    double alpha = 1.0+(double)n_alpha*0.1;  // alpha determines the relative radius of the left vessel 

                    // Height of of first-order vessels
                    QLength main_vert_length = 0.9*twicelambda*input_radius*pow(2.0,-1.0/3.0);
                    
                    // Height of the domain
                    QLength domain_side_length_y = 4.0*main_vert_length;

                    // Length of the domain
                    QLength domain_side_length_x = dimless_length*2.0*twicelambda*input_radius;

                    // Set the number of trials and the maximum beta value
                    int max_trials = 100;
                    int max_beta = 35;
                    double beta_step = 0.5;  // the number of values on the x-axis per beta, i.e., x-axis resolution

                    // Iterate over beta values
                    double beta = 0.5;
                    while (beta <= max_beta)
                    {   
                        // Run the trials                    
                        int n_trial = 1;
                        while (n_trial <= max_trials)
                        {  
                            // Generate the network
                            std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapy(order, main_vert_length, input_radius, alpha, twicelambda);

                            // Identify input and output nodes and assign them properties
                            VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(0.0_um,2.0*main_vert_length));
                            VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(domain_side_length_x, 2.0*main_vert_length));
                            p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
                            // p_network->GetNode(0)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));  // = 3333.06
                            p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
                            p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
                            // p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));  // = 1999.84
                            p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);

                            // Set the h-solver
                            std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                            std::string solver_name;
                            if (h_solver==1)
                            {
                                solver_name = "ConstantHaematocrit/"; 
                                std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                                auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                                p_haematocrit_solver->SetVesselNetwork(p_network);
                                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                                p_abstract_haematocrit_solver = p_haematocrit_solver;     
                            }
                            else if (h_solver==2)
                            {
                                solver_name = "PriesHaematocrit/"; 
                                std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                                auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                                p_haematocrit_solver->SetVesselNetwork(p_network);
                                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                                p_abstract_haematocrit_solver = p_haematocrit_solver;                
                            }
                            else if (h_solver==3)
                            {
                                solver_name = "MemoryHaematocrit/"; 
                                std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                                auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                                p_haematocrit_solver->SetVesselNetwork(p_network);
                                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                                p_abstract_haematocrit_solver = p_haematocrit_solver;     
                            }
                            else if (h_solver==4)
                            {
                                solver_name = "FungHaematocrit/"; 
                                std::cout << "Now using FungHaematocritSolver..." << std::endl;
                                auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                                p_haematocrit_solver->SetVesselNetwork(p_network);
                                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                                p_abstract_haematocrit_solver = p_haematocrit_solver;     
                            }
                            else if (h_solver==5)
                            {
                                solver_name = "YangHaematocrit/"; 
                                std::cout << "Now using YangHaematocritSolver..." << std::endl;
                                auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();   
                                p_haematocrit_solver->SetVesselNetwork(p_network);
                                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                                p_abstract_haematocrit_solver = p_haematocrit_solver;     
                            }
                                            
                            // Store the lambda values for file name
                            std::stringstream lambda_stream;
                            lambda_stream << std::fixed << std::setprecision(0) << lambda;
                            std::string lambda_string = lambda_stream.str();
                            std::string lambda_directory = network_name + solver_name + "/Lambda" + lambda_string; 

                            // Set up the grid for the finite difference solver
                            auto p_grid = RegularGrid<2>::Create();
                            // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
                            // QLength grid_spacing(1.0*1_um);  // the simulation time gets quite long if you reduce the resolution further
                            // QLength grid_spacing = 10_um;
                            p_grid->SetSpacing(grid_spacing);
                            c_vector<unsigned, 3> dimensions;
                            dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                            dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                            dimensions[2] = 1;
                            p_grid->SetDimensions(dimensions);

                            // Choose the PDE
                            auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                            
                            // Set the diffusivity and decay terms
                            p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                            p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                            // Set up the discrete source
                            auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                            QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                                    GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                            QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                                    Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                            p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                            p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                            p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                            p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                            // Set up the finite difference solver for oxygen (which handles everything)
                            auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                            p_oxygen_solver->SetPde(p_oxygen_pde);
                            p_oxygen_solver->SetLabel("oxygen");
                            p_oxygen_solver->SetGrid(p_grid);

                            // Set up the viscosity calculator
                            auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                            p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                            p_viscosity_calculator->SetVesselNetwork(p_network);
                            p_viscosity_calculator->Calculate();
                            
                            // Set up the impedance calculator
                            auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                            p_impedance_calculator->SetVesselNetwork(p_network);
                            p_impedance_calculator->Calculate();

                            // Set up the flow solver
                            FlowSolver<2> flow_solver;
                            flow_solver.SetVesselNetwork(p_network);
                            flow_solver.SetUp();
   
                            // Conduct stochastic pruning (thinner vessels are more likely to be pruned)
                            p_network->StochasticPruning(beta, false);  

                            // Set filename
                            std::stringstream alpha_stream;
                            std::stringstream trial_stream;
                            std::stringstream beta_stream;
                            alpha_stream << std::fixed << std::setprecision(2) << alpha;
                            trial_stream << std::fixed << std::setprecision(0) << n_trial;
                            beta_stream << std::fixed << std::setprecision(2) << beta;
                            std::string alpha_string = alpha_stream.str();
                            std::string trial_string = trial_stream.str();
                            std::string beta_string = beta_stream.str();
                            std::string file_name = lambda_directory + "/Alpha" + alpha_string + "/Beta" + beta_string + "/Trial" + trial_string;
                            std::string str_directory_name = file_name;
                            auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                            // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance. The loop runs the solvers, calculates the max. difference in the network between the previous H of a vessel and the new H in the vessel, and repeats this until all vessels converge on a value within a certain tolerance.
                            unsigned max_iter = 100;  // simulation fails if it doesn't converge in these many iterations
                            double tolerance2 = 0.01; //1.e-10
                            std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                            std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                            for (unsigned idx=0;idx<max_iter;idx++)
                            {
                                // Run the solvers
                                p_impedance_calculator->Calculate();
                                flow_solver.SetUp();
                                flow_solver.Solve();
                                p_abstract_haematocrit_solver->Calculate();
                                p_viscosity_calculator->Calculate();

                                // Check for convergence 
                                double max_difference = 0.0;
                                double h_for_max = 0.0;
                                double prev_for_max = 0.0;
                                for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)  // for all the segments 
                                {
                                    // Set segments with no flow to have no haematocrit
                                    if (fabs(segments[segment_index]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                    {
                                        segments[segment_index]->GetFlowProperties()->SetHaematocrit(0.0);
                                    }                                    
                                    double current_haematocrit = segments[segment_index]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                                    double difference = std::abs(current_haematocrit - previous_haematocrit[segment_index]);  // difference in new H in segment and old H
                                    if(difference>max_difference)  // if the difference is greater than previous max.
                                    {
                                        max_difference = difference;
                                        h_for_max = current_haematocrit;
                                        prev_for_max = previous_haematocrit[segment_index];
                                    }
                                    previous_haematocrit[segment_index] = current_haematocrit;
                                }
                                std::cout << "Segment H. at max difference: " << h_for_max << ", Prev. segment H at max difference:" << prev_for_max << std::endl;
                                if(max_difference<=tolerance2)  
                                {
                                    std::cout << "Converged after: " << idx << " iterations." <<  std::endl;
                                    broken_solver = 0;
                                    break;
                                }
                                else
                                {
                                    if(idx%1==0)
                                    {
                                        std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                                        // std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                                        // std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                                        // p_network->Write(output_file);
                                    }
                                }
                                if(idx==max_iter-1)  // if there is no convergence after all the iterations, print the error message
                                {
                                    std::cout << "Problem encountered in Dichotomous Network with h_solver = " << h_solver << " using n_alpha = " << n_alpha << " and lambda = " << lambda << " in trial = " << n_trial << std::endl;
                                    // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                                    error_log << "\n Problem encountered in Dichotomous Network with h_solver = " << h_solver << " using n_alpha = " << n_alpha << " and lambda = " << lambda << " in trial = " << n_trial; 
                                    broken_solver = 1;
                                    break;
                                }
                            }

                            // If solver doesn't converge, move on to next one
                            if (broken_solver == 1)
                            {
                                continue;
                            }

                            // Run the simulation 
                            SimulationTime::Instance()->SetStartTime(0.0);
                            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                            auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                            p_microvessel_solver->SetVesselNetwork(p_network);
                            p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                            p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                            p_microvessel_solver->Run();

                            // Print the average oxygenation
                            std::vector<double> solution = p_oxygen_solver->GetSolution();
                            double average_oxygen = 0.0;
                            for(unsigned jdx=0;jdx<solution.size();jdx++)
                            {
                                average_oxygen += solution[jdx];
                            }
                            average_oxygen /= double(solution.size());
                            std::cout << "Average oxygen: " << average_oxygen << std::endl;

                            // Write the PQs
                            outfile.open("forking_stochastic_pruning_perfusion_quotients.txt", std::ios_base::app);
                            outfile << network_name << " " << solver_name << " " << lambda_string << " " << alpha_string << " " << beta_string << " " << trial_string << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                            outfile.close();
                            
                            // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                            std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                            p_network->Write(output_file);

                            // // ??? Assuming radius doesn't exceed threshold, write the network to the output file
                            // if(current_radius_threshold>12.9_um && current_radius_threshold<13.1_um)
                            // {
                            //     p_network->Write(output_file_final_RT);
                            // }
                            n_trial = n_trial + 1;  // move onto next trial

                            // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                            ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                            ParameterCollection::Instance()->Destroy();
                            BaseUnits::Instance()->Destroy();
                            SimulationTime::Instance()->Destroy();
                        }
                        // Move on to next beta
                        beta = beta + beta_step;
                    }
                }
            }
        }
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a multi-generation forking network (with random daughter vessels being assigned greater radii) on a PDE grid with different haematocrit splitting rules and radius threshold pruning
    void xTestRandomDichotomousNetworkWithPruningAndFlow2D()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";

        // Set network in filename
        std::string network_name = "TestRandomDichotomousNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("forking_random_radius_threshold_pruning_perfusion_quotients.txt");
        outfile.close();

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=4; h_solver++)
        {   
            // Generate the network for various lambdas
            for (unsigned k_aux=1; k_aux<2; k_aux++)  // generate the network for lambda = 4
            {
                // // Run the simulation with different heterogeneities
                // for (unsigned n_alpha=0; n_alpha<=max_alpha; n_alpha++)
                // { 
                // Set key vessel parameters
                double dimless_length = 1.0;  

                // ??? Non-dimensionalising the length
                for(unsigned i_aux=1; i_aux<order+1; i_aux++)
                {
                    dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
                }

                // Set up the reference length for the simulation
                QLength reference_length(1.0_um);
                BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

                // Set input radius
                QLength input_radius(10.0 *GenericParameters::mpCapillaryRadius->GetValue());  // = 50_um
                // VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);

                // Set the viscosity
                QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
                // QDynamicViscosity viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();  // = 0.0012

                // Set the inlet and initial haematocrit
                double initial_haematocrit = 0.45;

                // Get pre-RT Perfusion Quotient
                QFlowRate threshold = 1.e-12*unit::metre_cubed_per_second;  // set flow threshold for what counts as a perfused vessel

                // Generate the networks
                VesselNetworkGenerator<2> network_generator;

                // Set lambda
                double lambda;  // lambda = length/diameter
                double twicelambda;  // used as an input parameter
                lambda = 2.0+double(k_aux)*2.0;
                twicelambda = 2.0*lambda;

                // Run the simulation with different heterogeneities
                for (unsigned n_alpha=0; n_alpha<=max_alpha; n_alpha++)
                { 
                    // Set the number of trials to generate the network
                    int max_trials = 100;

                    // Run the trials                    
                    int n_trial = 1;
                    while (n_trial <= max_trials)
                    {  
                        // Set up flag for broken solver
                        unsigned broken_trial = 0;

                        // Set alpha
                        double alpha = 1.0+(double)n_alpha*0.1;  // alpha determines the relative radius of the left vessel 

                        // Height of of first-order vessels
                        QLength main_vert_length = 0.9*twicelambda*input_radius*pow(2.0,-1.0/3.0);
                        
                        // Height of the domain
                        QLength domain_side_length_y = 4.0*main_vert_length;

                        // Length of the domain
                        QLength domain_side_length_x = dimless_length*2.0*twicelambda*input_radius;

                        // Generate the network
                        std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapyRandomHeterogeneity(order, main_vert_length, input_radius, alpha, twicelambda);

                        // Identify input and output nodes and assign them properties
                        VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(0.0_um,2.0*main_vert_length));
                        VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(domain_side_length_x, 2.0*main_vert_length));
                        p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
                        // p_network->GetNode(0)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));  // = 3333.06
                        p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
                        p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
                        // p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));  // = 1999.84
                        p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);

                        // Set the h-solver
                        std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                        std::string solver_name;
                        if (h_solver==1)
                        {
                            solver_name = "ConstantHaematocrit/"; 
                            std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                            auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                            p_haematocrit_solver->SetVesselNetwork(p_network);
                            p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                            p_abstract_haematocrit_solver = p_haematocrit_solver;     
                        }
                        else if (h_solver==2)
                        {
                            solver_name = "PriesHaematocrit/"; 
                            std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                            auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                            // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                            p_haematocrit_solver->SetVesselNetwork(p_network);
                            p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                            p_abstract_haematocrit_solver = p_haematocrit_solver;                
                        }
                        else if (h_solver==3)
                        {
                            solver_name = "MemoryHaematocrit/"; 
                            std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                            auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                            // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                            p_haematocrit_solver->SetVesselNetwork(p_network);
                            p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                            p_abstract_haematocrit_solver = p_haematocrit_solver;     
                        }
                        else if (h_solver==4)
                        {
                            solver_name = "FungHaematocrit/"; 
                            std::cout << "Now using FungHaematocritSolver..." << std::endl;
                            auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                            // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                            p_haematocrit_solver->SetVesselNetwork(p_network);
                            p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                            p_abstract_haematocrit_solver = p_haematocrit_solver;     
                        }
                        else if (h_solver==5)
                        {
                            solver_name = "YangHaematocrit/"; 
                            std::cout << "Now using YangHaematocritSolver..." << std::endl;
                            auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                            // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();   
                            p_haematocrit_solver->SetVesselNetwork(p_network);
                            p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                            p_abstract_haematocrit_solver = p_haematocrit_solver;     
                        }
                                        
                        // Store the lambda values for file name
                        std::stringstream lambda_stream;
                        lambda_stream << std::fixed << std::setprecision(0) << lambda;
                        std::string lambda_string = lambda_stream.str();
                        std::string lambda_directory = network_name + solver_name + "/Lambda" + lambda_string; 

                        // Set up the grid for the finite difference solver
                        auto p_grid = RegularGrid<2>::Create();
                        // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
                        // QLength grid_spacing(1.0*1_um);  // the simulation time gets quite long if you reduce the resolution further
                        // QLength grid_spacing = 10_um;
                        p_grid->SetSpacing(grid_spacing);
                        c_vector<unsigned, 3> dimensions;
                        dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                        dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                        dimensions[2] = 1;
                        p_grid->SetDimensions(dimensions);

                        // Choose the PDE
                        auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                        
                        // Set the diffusivity and decay terms
                        p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                        p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                        // Set up the discrete source
                        auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                        QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                                GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                        QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                                Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                        p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                        p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                        p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                        // Set up the finite difference solver for oxygen (which handles everything)
                        auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                        p_oxygen_solver->SetPde(p_oxygen_pde);
                        p_oxygen_solver->SetLabel("oxygen");
                        p_oxygen_solver->SetGrid(p_grid);

                        // Set up the viscosity calculator
                        auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                        p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                        p_viscosity_calculator->SetVesselNetwork(p_network);
                        p_viscosity_calculator->Calculate();
                        
                        // Set up the impedance calculator
                        auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                        p_impedance_calculator->SetVesselNetwork(p_network);
                        p_impedance_calculator->Calculate();

                        // Set up the flow solver
                        FlowSolver<2> flow_solver;
                        flow_solver.SetVesselNetwork(p_network);
                        flow_solver.SetUp();

                        // Set up pruning parameters
                        QLength max_radius_to_kill = 31_um;  // set threshold up to which vessels should be pruned
                        QLength radius_step = 1_um;
                        QLength current_radius_threshold = 0_um;

                        // Prunes vessels from thinnest up to simulate RT
                        while (current_radius_threshold <= max_radius_to_kill)
                        {         
                            // Display status message
                            double current_radius_threshold_um = current_radius_threshold*1000000;  // in micrometres
                            std::cout << "Now pruning up to vessels with radius = " << current_radius_threshold_um << " um" << std::endl;

                            // Remove vessels under radius threshold
                            p_network->RemoveThinVessels(current_radius_threshold, false);  

                            // Set filename
                            std::stringstream alpha_stream;
                            std::stringstream threshold_stream;
                            std::stringstream trial_stream;
                            alpha_stream << std::fixed << std::setprecision(2) << alpha; 
                            trial_stream << std::fixed << std::setprecision(0) << n_trial;
                            threshold_stream << std::fixed << std::setprecision(0) << current_radius_threshold_um;
                            std::string alpha_string = alpha_stream.str();                           
                            std::string trial_string = trial_stream.str();
                            std::string threshold_string = threshold_stream.str();
                            std::string file_name = lambda_directory + "/Alpha" + alpha_string + "/Trial" + trial_string+ "/RadiusThreshold" + threshold_string;
                            std::string str_directory_name = file_name;
                            auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                            // outfile << alpha << " " << current_radius_threshold << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 

                            // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance. The loop runs the solvers, calculates the max. difference in the network between the previous H of a vessel and the new H in the vessel, and repeats this until all vessels converge on a value within a certain tolerance.
                            unsigned max_iter = 100;  // simulation fails if it doesn't converge in these many iterations
                            double tolerance2 = 0.01;
                            std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                            std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                            for (unsigned idx=0;idx<max_iter;idx++)
                            {
                                // Run the solvers
                                p_impedance_calculator->Calculate();
                                flow_solver.SetUp();
                                flow_solver.Solve();
                                p_abstract_haematocrit_solver->Calculate();
                                p_viscosity_calculator->Calculate();

                                // Check for convergence 
                                double max_difference = 0.0;
                                double h_for_max = 0.0;
                                double prev_for_max = 0.0;
                                for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)  // for all the segments 
                                {
                                    // Set segments with no flow to have no haematocrit
                                    if (fabs(segments[segment_index]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                    {
                                        segments[segment_index]->GetFlowProperties()->SetHaematocrit(0.0);
                                    }
                                    
                                    // Check for convergence
                                    double current_haematocrit = segments[segment_index]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                                    double difference = std::abs(current_haematocrit - previous_haematocrit[segment_index]);  // difference in new H in segment and old H
                                    if(difference>max_difference)  // if the difference is greater than previous max.
                                    {
                                        max_difference = difference;
                                        h_for_max = current_haematocrit;
                                        prev_for_max = previous_haematocrit[segment_index];
                                    }
                                    previous_haematocrit[segment_index] = current_haematocrit;
                                }
                                std::cout << "Segment H. at max difference: " << h_for_max << ", Prev. segment H at max difference:" << prev_for_max << std::endl;
                                if(max_difference<=tolerance2)  
                                {
                                    std::cout << "Converged after: " << idx << " iterations." <<  std::endl;
                                    broken_trial = 0;
                                    break;
                                }
                                else
                                {
                                    if(idx%1==0)
                                    {
                                        std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                                        std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                                        p_network->Write(output_file);
                                    }
                                }
                                if(idx==max_iter-1)  // if there is no convergence after all the iterations, print the error message
                                {
                                    std::cout << "Problem encountered in Dichotomous Network with h_solver = " << h_solver << " using n_alpha = " << n_alpha << " and lambda = " << lambda << " in trial = " << n_trial << std::endl;
                                    // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                                    error_log << "\n Problem encountered in Dichotomous Network with h_solver = " << h_solver << " using n_alpha = " << n_alpha << " and lambda = " << lambda  << " in trial = " << n_trial << " and radius threshold = " << current_radius_threshold; 
                                    broken_trial = 1;
                                    break;
                                }
                            }

                            // If solver doesn't converge, restart trial
                            if (broken_trial == 1)
                            {
                                break;
                            }

                            // Run the simulation 
                            SimulationTime::Instance()->SetStartTime(0.0);
                            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                            auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                            p_microvessel_solver->SetVesselNetwork(p_network);
                            p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                            p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                            p_microvessel_solver->Run();

                            // Print the average oxygenation
                            std::vector<double> solution = p_oxygen_solver->GetSolution();
                            double average_oxygen = 0.0;
                            for(unsigned jdx=0;jdx<solution.size();jdx++)
                            {
                                average_oxygen += solution[jdx];
                            }
                            average_oxygen /= double(solution.size());
                            std::cout << "Average oxygen: " << average_oxygen << std::endl;

                            // Write the PQs
                            outfile.open("forking_random_radius_threshold_pruning_perfusion_quotients.txt", std::ios_base::app);
                            outfile << network_name << " " << solver_name << " " << lambda_string << " " << alpha_string << " " << trial_string << " " << threshold_string << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                            outfile.close();
                            
                            // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                            std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                            p_network->Write(output_file);

                            // // ??? Assuming radius doesn't exceed threshold, write the network to the output file
                            // if(current_radius_threshold>12.9_um && current_radius_threshold<13.1_um)
                            // {
                            //     p_network->Write(output_file_final_RT);
                            // }
                            current_radius_threshold = current_radius_threshold + radius_step;  // move onto next radius threshold for next step

                            // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                            ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                            ParameterCollection::Instance()->Destroy();
                            BaseUnits::Instance()->Destroy();
                            SimulationTime::Instance()->Destroy();
                        }
                        
                        // If trial broke, restart, otherwise move on
                        if (broken_trial == 0)
                        {
                            n_trial = n_trial + 1;
                        }
                    }
                }
            }
        }
        // Close the output file

        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a multi-generation forking network (with random daughter vessels being assigned greater radii) on a PDE grid with different h-splitting rules and stochastic pruning
    void xTestRandomDichotomousNetworkWithStochasticPruningAndFlow2D()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";

        // Set network in filename
        std::string network_name = "TestRandomDichotomousNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("forking_random_stochastic_pruning_perfusion_quotients.txt", std::ios_base::app);
        outfile.close();

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=3; h_solver<=3; h_solver++)
        {   
            // Generate the network for various lambdas
            for (unsigned k_aux=1; k_aux<2; k_aux++)  // generate the network for lambda = 4
            {
                // Set key vessel parameters
                double dimless_length = 1.0;  

                // ??? Non-dimensionalising the length
                for(unsigned i_aux=1; i_aux<order+1; i_aux++)
                {
                    dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
                }
   
                // Set up the reference length for the simulation
                QLength reference_length(1.0_um);
                BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

                // Set input radius
                QLength input_radius(10.0 *GenericParameters::mpCapillaryRadius->GetValue());  // = 50_um
                // VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);

                // Set the viscosity
                QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
                // QDynamicViscosity viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();  // = 0.0012

                // Set the inlet and initial haematocrit
                double initial_haematocrit = 0.45;

                // Set threshold for perfusion quotient
		        QFlowRate threshold = 1.e-12*unit::metre_cubed_per_second;  // set flow threshold for what counts as a perfused vessel

                // Generate the networks
                VesselNetworkGenerator<2> network_generator;

                // Set lambda
                double lambda;  // lambda = length/diameter
                double twicelambda;  // used as an input parameter
                lambda = 2.0+double(k_aux)*2.0;
                twicelambda = 2.0*lambda;

                // Run the simulation with different heterogeneities
                for (unsigned n_alpha=3; n_alpha<=max_alpha; n_alpha++)
                { 
                    // Set up flag for broken solver
                    unsigned broken_solver = 0;

                    // Set alpha
                    double alpha = 1.0+(double)n_alpha*0.1;  // alpha determines the relative radius of the left vessel 

                    // Height of of first-order vessels
                    QLength main_vert_length = 0.9*twicelambda*input_radius*pow(2.0,-1.0/3.0);
                    
                    // Height of the domain
                    QLength domain_side_length_y = 4.0*main_vert_length;

                    // Length of the domain
                    QLength domain_side_length_x = dimless_length*2.0*twicelambda*input_radius;

                    // Set the number of trials and the maximum beta value
                    int max_trials = 100;
                    int max_beta = 35;
                    
                    // Iterate over beta values
                    int beta = 1;
                    while (beta <= max_beta)
                    {
                        // Run the trials                    
                        int n_trial = 1;
                        while (n_trial <= max_trials)
                        {  
                            // Generate the network
                            std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapyRandomHeterogeneity(order, main_vert_length, input_radius, alpha, twicelambda);

                            // Identify input and output nodes and assign them properties
                            VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(0.0_um,2.0*main_vert_length));
                            VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(domain_side_length_x, 2.0*main_vert_length));
                            p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
                            // p_network->GetNode(0)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));  // = 3333.06
                            p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
                            p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
                            // p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));  // = 1999.84
                            p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);

                            // Set the h-solver
                            std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                            std::string solver_name;
                            if (h_solver==1)
                            {
                                solver_name = "ConstantHaematocrit/"; 
                                std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                                auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                                p_haematocrit_solver->SetVesselNetwork(p_network);
                                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                                p_abstract_haematocrit_solver = p_haematocrit_solver;     
                            }
                            else if (h_solver==2)
                            {
                                solver_name = "PriesHaematocrit/"; 
                                std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                                auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                                p_haematocrit_solver->SetVesselNetwork(p_network);
                                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                                p_abstract_haematocrit_solver = p_haematocrit_solver;                
                            }
                            else if (h_solver==3)
                            {
                                solver_name = "MemoryHaematocrit/"; 
                                std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                                auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                                p_haematocrit_solver->SetVesselNetwork(p_network);
                                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                                p_abstract_haematocrit_solver = p_haematocrit_solver;     
                            }
                            else if (h_solver==4)
                            {
                                solver_name = "FungHaematocrit/"; 
                                std::cout << "Now using FungHaematocritSolver..." << std::endl;
                                auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                                // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                                p_haematocrit_solver->SetVesselNetwork(p_network);
                                p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                                p_abstract_haematocrit_solver = p_haematocrit_solver;     
                            }
                                            
                            // Store the lambda values for file name
                            std::stringstream lambda_stream;
                            lambda_stream << std::fixed << std::setprecision(0) << lambda;
                            std::string lambda_string = lambda_stream.str();
                            std::string lambda_directory = network_name + solver_name + "/Lambda" + lambda_string; 

                            // Set up the grid for the finite difference solver
                            auto p_grid = RegularGrid<2>::Create();
                            // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
                            // QLength grid_spacing(1.0*1_um);  // the simulation time gets quite long if you reduce the resolution further
                            // QLength grid_spacing = 10_um;
                            p_grid->SetSpacing(grid_spacing);
                            c_vector<unsigned, 3> dimensions;
                            dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                            dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                            dimensions[2] = 1;
                            p_grid->SetDimensions(dimensions);

                            // Choose the PDE
                            auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                            
                            // Set the diffusivity and decay terms
                            p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                            p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                            // Set up the discrete source
                            auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                            QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                                    GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                            QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                                    Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                            p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                            p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                            p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                            p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                            // Set up the finite difference solver for oxygen (which handles everything)
                            auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                            p_oxygen_solver->SetPde(p_oxygen_pde);
                            p_oxygen_solver->SetLabel("oxygen");
                            p_oxygen_solver->SetGrid(p_grid);

                            // Set up the viscosity calculator
                            auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                            p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                            p_viscosity_calculator->SetVesselNetwork(p_network);
                            p_viscosity_calculator->Calculate();
                            
                            // Set up the impedance calculator
                            auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                            p_impedance_calculator->SetVesselNetwork(p_network);
                            p_impedance_calculator->Calculate();

                            // Set up the flow solver
                            FlowSolver<2> flow_solver;
                            flow_solver.SetVesselNetwork(p_network);
                            flow_solver.SetUp();
   
                            // Conduct stochastic pruning (thinner vessels are more likely to be pruned)
                            p_network->StochasticPruning(beta, false);  

                            // Set filename
                            std::stringstream alpha_stream;
                            std::stringstream trial_stream;
                            std::stringstream beta_stream;
                            alpha_stream << std::fixed << std::setprecision(2) << alpha;
                            trial_stream << std::fixed << std::setprecision(0) << n_trial;
                            beta_stream << std::fixed << std::setprecision(0) << beta;
                            std::string alpha_string = alpha_stream.str();
                            std::string trial_string = trial_stream.str();
                            std::string beta_string = beta_stream.str();
                            std::string file_name = lambda_directory + "/Alpha" + alpha_string + "/Beta" + beta_string + "/Trial" + trial_string;
                            std::string str_directory_name = file_name;
                            auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                            // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance. The loop runs the solvers, calculates the max. difference in the network between the previous H of a vessel and the new H in the vessel, and repeats this until all vessels converge on a value within a certain tolerance.
                            unsigned max_iter = 100;  // simulation fails if it doesn't converge in these many iterations
                            double tolerance2 = 0.01; //1.e-10
                            std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                            std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                            for (unsigned idx=0;idx<max_iter;idx++)
                            {
                                // Run the solvers
                                p_impedance_calculator->Calculate();
                                flow_solver.SetUp();
                                flow_solver.Solve();
                                p_abstract_haematocrit_solver->Calculate();
                                p_viscosity_calculator->Calculate();

                                // Check for convergence 
                                double max_difference = 0.0;
                                double h_for_max = 0.0;
                                double prev_for_max = 0.0;
                                for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)  // for all the segments 
                                {
                                    // Set segments with no flow to have no haematocrit
                                    if (fabs(segments[segment_index]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                    {
                                        segments[segment_index]->GetFlowProperties()->SetHaematocrit(0.0);
                                    }                                    
                                    double current_haematocrit = segments[segment_index]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                                    double difference = std::abs(current_haematocrit - previous_haematocrit[segment_index]);  // difference in new H in segment and old H
                                    if(difference>max_difference)  // if the difference is greater than previous max.
                                    {
                                        max_difference = difference;
                                        h_for_max = current_haematocrit;
                                        prev_for_max = previous_haematocrit[segment_index];
                                    }
                                    previous_haematocrit[segment_index] = current_haematocrit;
                                }
                                std::cout << "Segment H. at max difference: " << h_for_max << ", Prev. segment H at max difference:" << prev_for_max << std::endl;
                                if(max_difference<=tolerance2)  
                                {
                                    std::cout << "Converged after: " << idx << " iterations." <<  std::endl;
                                    broken_solver = 0;
                                    break;
                                }
                                else
                                {
                                    if(idx%1==0)
                                    {
                                        std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                                        // std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                                        // std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                                        // p_network->Write(output_file);
                                    }
                                }
                                if(idx==max_iter-1)  // if there is no convergence after all the iterations, print the error message
                                {
                                    std::cout << "Problem encountered in Dichotomous Network with h_solver = " << h_solver << " using n_alpha = " << n_alpha << " and lambda = " << lambda << " in trial = " << n_trial << std::endl;
                                    // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                                    error_log << "\n Problem encountered in Dichotomous Network with h_solver = " << h_solver << " using n_alpha = " << n_alpha << " and lambda = " << lambda << " in trial = " << n_trial; 
                                    broken_solver = 1;
                                    break;
                                }
                            }

                            // If solver doesn't converge, move on to next one
                            if (broken_solver == 1)
                            {
                                continue;
                            }

                            // Run the simulation 
                            SimulationTime::Instance()->SetStartTime(0.0);
                            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                            auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                            p_microvessel_solver->SetVesselNetwork(p_network);
                            p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                            p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                            p_microvessel_solver->Run();

                            // Print the average oxygenation
                            std::vector<double> solution = p_oxygen_solver->GetSolution();
                            double average_oxygen = 0.0;
                            for(unsigned jdx=0;jdx<solution.size();jdx++)
                            {
                                average_oxygen += solution[jdx];
                            }
                            average_oxygen /= double(solution.size());
                            std::cout << "Average oxygen: " << average_oxygen << std::endl;

                            // Write the PQs
                            outfile.open("forking_random_stochastic_pruning_perfusion_quotients.txt", std::ios_base::app);
                            outfile << network_name << " " << solver_name << " " << lambda_string << " " << alpha_string << " " << beta_string << " " << trial_string << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                            outfile.close();
                            
                            // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                            std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                            p_network->Write(output_file);

                            // // ??? Assuming radius doesn't exceed threshold, write the network to the output file
                            // if(current_radius_threshold>12.9_um && current_radius_threshold<13.1_um)
                            // {
                            //     p_network->Write(output_file_final_RT);
                            // }
                            n_trial = n_trial + 1;  // move onto next trial

                            // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                            ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                            ParameterCollection::Instance()->Destroy();
                            BaseUnits::Instance()->Destroy();
                            SimulationTime::Instance()->Destroy();
                        }
                        // Move on to next beta
                        beta = beta + 1;
                    }
                }
            }
        }
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a multi-generation forking network on a PDE grid with different h-splitting rules and stochastic pruning
    void xTestDichotomousNetworkWithIndividualPruningAndFlow2D()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";

        // Set network in filename
        std::string network_name = "TestDichotomousNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("forking_nonrandom_individual_pruning_perfusion_quotients.txt", std::ios_base::app);
        outfile.close();

        // Define the key pruning parameters
        unsigned n_vessels = 252;  // number of non-inlet/outlet vessels from which to select ones to make thin
        // double percToKill = 0.2;  // percentage of vessels to kill
        // unsigned ToBeKilled = (unsigned)(percToKill*n_vessels);  // number to kill
        unsigned ToBeKilled = n_vessels;  // number to kill

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=4; h_solver++)
        {   
            // Generate the network for various lambdas
            for (unsigned k_aux=1; k_aux<2; k_aux++)  // generate the network for lambda = 4
            {
                // Set key vessel parameters
                double dimless_length = 1.0;  

                // ??? Non-dimensionalising the length
                for(unsigned i_aux=1; i_aux<order+1; i_aux++)
                {
                    dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
                }
   
                // Set up the reference length for the simulation
                QLength reference_length(1.0_um);
                BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

                // Set input radius
                QLength input_radius(10.0 *GenericParameters::mpCapillaryRadius->GetValue());  // = 50_um
                // VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);

                // Set the viscosity
                QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
                // QDynamicViscosity viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();  // = 0.0012

                // Set the inlet and initial haematocrit
                double initial_haematocrit = 0.45;

                // Set threshold for perfusion quotient
		        QFlowRate threshold = 1.e-12*unit::metre_cubed_per_second;  // set flow threshold for what counts as a perfused vessel

                // Generate the networks
                VesselNetworkGenerator<2> network_generator;

                // Set lambda
                double lambda;  // lambda = length/diameter
                double twicelambda;  // used as an input parameter
                lambda = 2.0+double(k_aux)*2.0;
                twicelambda = 2.0*lambda;

                // Run the simulation with different heterogeneities
                for (unsigned n_alpha=0; n_alpha<=max_alpha; n_alpha++)
                { 
                    // Set up flag for broken solver
                    unsigned broken_solver = 0;

                    // Set alpha
                    double alpha = 1.0+(double)n_alpha*0.1;  // alpha determines the relative radius of the left vessel 

                    // Height of of first-order vessels
                    QLength main_vert_length = 0.9*twicelambda*input_radius*pow(2.0,-1.0/3.0);
                    
                    // Height of the domain
                    QLength domain_side_length_y = 4.0*main_vert_length;

                    // Length of the domain
                    QLength domain_side_length_x = dimless_length*2.0*twicelambda*input_radius;
            
                    // Generate the network
                    std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapy(order, main_vert_length, input_radius, alpha, twicelambda);

                    // Identify input and output nodes and assign them properties
                    VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(0.0_um,2.0*main_vert_length));
                    VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(domain_side_length_x, 2.0*main_vert_length));
                    p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
                    // p_network->GetNode(0)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));  // = 3333.06
                    p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
                    p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
                    // p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));  // = 1999.84
                    p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);

                    // Set the h-solver
                    std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                    std::string solver_name;
                    if (h_solver==1)
                    {
                        solver_name = "ConstantHaematocrit/"; 
                        std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==2)
                    {
                        solver_name = "PriesHaematocrit/"; 
                        std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                
                    }
                    else if (h_solver==3)
                    {
                        solver_name = "MemoryHaematocrit/"; 
                        std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==4)
                    {
                        solver_name = "FungHaematocrit/"; 
                        std::cout << "Now using FungHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                                            
                    // Store the lambda values for file name
                    std::stringstream lambda_stream;
                    lambda_stream << std::fixed << std::setprecision(0) << lambda;
                    std::string lambda_string = lambda_stream.str();
                    std::string lambda_directory = network_name + solver_name + "/Lambda" + lambda_string; 

                    // Set up the grid for the finite difference solver
                    auto p_grid = RegularGrid<2>::Create();
                    // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
                    // QLength grid_spacing(1.0*1_um);  // the simulation time gets quite long if you reduce the resolution further
                    // QLength grid_spacing = 10_um;
                    p_grid->SetSpacing(grid_spacing);
                    c_vector<unsigned, 3> dimensions;
                    dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                    dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                    dimensions[2] = 1;
                    p_grid->SetDimensions(dimensions);

                    // Choose the PDE
                    auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                    
                    // Set the diffusivity and decay terms
                    p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                    p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                    // Set up the discrete source
                    auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                    QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                            GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                    QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                            Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                    p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                    p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                    p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                    // Set up the finite difference solver for oxygen (which handles everything)
                    auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                    p_oxygen_solver->SetPde(p_oxygen_pde);
                    p_oxygen_solver->SetLabel("oxygen");
                    p_oxygen_solver->SetGrid(p_grid);

                    // Set up the viscosity calculator
                    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                    p_viscosity_calculator->SetVesselNetwork(p_network);
                    p_viscosity_calculator->Calculate();
                    
                    // Set up the impedance calculator
                    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                    p_impedance_calculator->SetVesselNetwork(p_network);
                    p_impedance_calculator->Calculate();

                    // Set up the flow solver
                    FlowSolver<2> flow_solver;
                    flow_solver.SetVesselNetwork(p_network);
                    flow_solver.SetUp();
   
                    // Prune all vessels up to specified dose 
                    for(unsigned KilledVessels=0; KilledVessels < ToBeKilled+1; KilledVessels++)
                    { 
                        // Display status message
                        std::cout << "Now killed " << KilledVessels << " vessels." << std::endl;  

                        // Set filename
                        std::stringstream alpha_stream;
                        std::stringstream kill_stream;
                        alpha_stream << std::fixed << std::setprecision(2) << alpha;
                        kill_stream << std::fixed << std::setprecision(0) << KilledVessels;
                        std::string alpha_string = alpha_stream.str();
                        std::string kill_string = kill_stream.str();
                        std::string file_name = lambda_directory + "/Alpha" + alpha_string + "/Kills" + kill_string;
                        std::string str_directory_name = file_name;
                        auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                        // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance. The loop runs the solvers, calculates the max. difference in the network between the previous H of a vessel and the new H in the vessel, and repeats this until all vessels converge on a value within a certain tolerance.
                        unsigned max_iter = 100;  // simulation fails if it doesn't converge in these many iterations
                        double tolerance2 = 0.01; //1.e-10
                        std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                        std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                        for (unsigned idx=0;idx<max_iter;idx++)
                        {
                            // Run the solvers
                            p_impedance_calculator->Calculate();
                            flow_solver.SetUp();
                            flow_solver.Solve();
                            p_abstract_haematocrit_solver->Calculate();
                            p_viscosity_calculator->Calculate();

                            // Check for convergence 
                            double max_difference = 0.0;
                            double h_for_max = 0.0;
                            double prev_for_max = 0.0;
                            for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)  // for all the segments 
                            {
                                // Set segments with no flow to have no haematocrit
                                if (fabs(segments[segment_index]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                {
                                    segments[segment_index]->GetFlowProperties()->SetHaematocrit(0.0);
                                }                                    
                                double current_haematocrit = segments[segment_index]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                                double difference = std::abs(current_haematocrit - previous_haematocrit[segment_index]);  // difference in new H in segment and old H
                                if(difference>max_difference)  // if the difference is greater than previous max.
                                {
                                    max_difference = difference;
                                    h_for_max = current_haematocrit;
                                    prev_for_max = previous_haematocrit[segment_index];
                                }
                                previous_haematocrit[segment_index] = current_haematocrit;
                            }
                            std::cout << "Segment H. at max difference: " << h_for_max << ", Prev. segment H at max difference:" << prev_for_max << std::endl;
                            if(max_difference<=tolerance2)  
                            {
                                std::cout << "Converged after: " << idx << " iterations." <<  std::endl;
                                broken_solver = 0;
                                break;
                            }
                            else
                            {
                                if(idx%1==0)
                                {
                                    std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                                    // std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                                    // std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                                    // p_network->Write(output_file);
                                }
                            }
                            if(idx==max_iter-1)  // if there is no convergence after all the iterations, print the error message
                            {
                                std::cout << "Problem encountered in Dichotomous Network with h_solver = " << h_solver << " using n_alpha = " << n_alpha << " and lambda = " << lambda << std::endl;
                                // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                                error_log << "\n Problem encountered in Dichotomous Network with h_solver = " << h_solver << " using n_alpha = " << n_alpha << " and lambda = " << lambda; 
                                broken_solver = 1;
                                break;
                            }
                        }

                        // If solver doesn't converge, move on to next one
                        if (broken_solver == 1)
                        {
                            continue;
                        }

                        // Run the simulation 
                        SimulationTime::Instance()->SetStartTime(0.0);
                        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                        auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                        p_microvessel_solver->SetVesselNetwork(p_network);
                        p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                        p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                        p_microvessel_solver->Run();

                        // Print the average oxygenation
                        std::vector<double> solution = p_oxygen_solver->GetSolution();
                        double average_oxygen = 0.0;
                        for(unsigned jdx=0;jdx<solution.size();jdx++)
                        {
                            average_oxygen += solution[jdx];
                        }
                        average_oxygen /= double(solution.size());
                        std::cout << "Average oxygen: " << average_oxygen << std::endl;

                        // Write the PQs
                        outfile.open("forking_nonrandom_individual_pruning_perfusion_quotients.txt", std::ios_base::app);
                        outfile << network_name << " " << solver_name << " " << lambda_string << " " << alpha_string << " " << KilledVessels << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                        outfile.close();
                        
                        // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                        p_network->Write(output_file);

                        // Remove the smallest vessel
                        QLength minimum_radius = input_radius;
                        unsigned int minimum_index = 0;
                        std::vector<std::shared_ptr<Vessel<2> > > vessels = p_network->GetVessels();
                        for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)  // for all the segments in the network
                        {
                            // Get the current segment's radius
                            QLength current_radius = vessels[vessel_index]->GetRadius();

                            // If the current radius is less than the minimum radius, record the new minimum
                            if (current_radius < minimum_radius)
                            {
                                minimum_radius = current_radius;
                                minimum_index = vessel_index;
                            }
                        }
                        p_network->RemoveVessel(vessels[minimum_index], true);  // remove the vessel

                        // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                        ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                        ParameterCollection::Instance()->Destroy();
                        BaseUnits::Instance()->Destroy();
                        SimulationTime::Instance()->Destroy();
                    }
                }
            }
        }
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a four-generation forking network on a PDE grid with different haematocrit splitting rules to compare with Hyakutake et al. (Microvascular Research, 2022)
    void xTestFourGenerationNetworkWithFlow2D()
    {
        // Run the simulation with different heterogeneities
        for (unsigned n_alpha=0; n_alpha<=0; n_alpha++)
        {
            double alpha = 1.0+(double)n_alpha*0.1;  // alpha determines the relative radius of the left vessel 
            
            // Run the simulation with different solvers of interest
            for (unsigned h_solver=3; h_solver<=3; h_solver++)
            {      
                // Set up the reference length for the simulation
                QLength reference_length(1.0_um);
                BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

                // Set key vessel parameters
                // QLength vessel_length(100.0_um);
                // QLength vessel_height = (pow(2,0.5)*vessel_length)*0.5;
                double dimless_length = 1.0;  

                // ???
                for(unsigned i_aux=1; i_aux<order+1; i_aux++)
                {
                    dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
                }

                // Set input radius
                QLength input_radius(2.0 *GenericParameters::mpCapillaryRadius->GetValue());  // = 5_um
                // QLength input_radius(10.0 *GenericParameters::mpCapillaryRadius->GetValue());  // = 50_um
                // VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);

                // Set the viscosity
                QDynamicViscosity viscosity = 5.e-3*unit::poiseuille;
                // QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
                // QDynamicViscosity viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();  // = 0.0012

                // Set the inlet and initial haematocrit
                // p_haematocrit_calculator->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));  // = 0.45
                // double initial_haematocrit = 0.45;
                double initial_haematocrit = 0.15;

                // Generate the networks
                VesselNetworkGenerator<2> network_generator;
                double lambda;  // lambda = length/diameter
                double twicelambda;  // used as an input parameter
                // for (unsigned k_aux=1; k_aux<5; k_aux++)  // generate the network for various lambdas
                for (unsigned k_aux=1; k_aux<2; k_aux++)  // generate the network for lambda = 4
                {
                    lambda = 6.0+double(k_aux)*2.0;
                    // lambda = 9.0;
                    twicelambda = 2.0*lambda;

                    // Height of of first-order vessels
                    QLength main_vert_length = 0.9*twicelambda*input_radius*pow(2.0,-1.0/3.0);
                    
                    // Height of the domain
                    QLength domain_side_length_y = 4.0*main_vert_length;

                    // Length of the domain
                    QLength domain_side_length_x = dimless_length*2.0*twicelambda*input_radius;

                    // Choose the number of bifurcating generations (inlet and outlet vessels don't count)
                    unsigned order = 4;

                    // Set file name based on haematocrit solver
                    std::ostringstream strs;
                    strs << std::fixed << std::setprecision( 1 );
                    if (h_solver==1)
                    {
                        strs << "TestDichotomousNetwork/ConstantHaematocrit/Alpha" << alpha << "/NoPruning/GridSpacingEquals" << grid_spacing;
                    }
                    else if (h_solver==2)
                    {
                        strs << "TestDichotomousNetwork/PriesHaematocrit/Alpha" << alpha << "/NoPruning/GridSpacingEquals" << grid_spacing;
                    }
                    else if (h_solver==3)
                    {
                        strs << "TestDichotomousNetwork/MemoryHaematocrit/Alpha" << alpha << "/NoPruning/GridSpacingEquals" << grid_spacing;
                    }
                    else if (h_solver==4)
                    {
                        strs << "TestDichotomousNetwork/FungHaematocrit/Alpha" << alpha << "/NoPruning/GridSpacingEquals" << grid_spacing;
                    }
                    else if (h_solver==5)
                    {
                        strs << "TestDichotomousNetwork/YangHaematocrit/Alpha" << alpha << "/NoPruning/GridSpacingEquals" << grid_spacing;
                    }
                    std::string str_directory_name = strs.str();
                    auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                    // Generate the network
                    std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapy(order, main_vert_length, input_radius, alpha, twicelambda);

                    // Identify input and output nodes and assign them properties
                    VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(0.0_um,2.0*main_vert_length));
                    VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(domain_side_length_x, 2.0*main_vert_length));
                    p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
                    // p_network->GetNode(0)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));  // = 3333.06
                    p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
                    p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
                    // p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));  // = 1999.84
                    p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);

                    // Set up the grid for the finite difference solver
                    auto p_grid = RegularGrid<2>::Create();
                    // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
                    // QLength grid_spacing(1.0*1_um);  // the simulation time gets quite long if you reduce the resolution further
                    QLength grid_spacing = 10_um;
                    p_grid->SetSpacing(grid_spacing);
                    c_vector<unsigned, 3> dimensions;
                    dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                    dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                    dimensions[2] = 1;
                    p_grid->SetDimensions(dimensions);
                    // std::shared_ptr<Part<2> > p_domain = Part<2>::Create();
                    // p_domain->AddRectangle(domain_x, domain_y);
                    // std::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
                    // p_grid->GenerateFromPart(p_domain, grid_spacing);  // set domain and grid spacing

                    // Choose the PDE
                    auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                    
                    // Set the diffusivity and decay terms
                    p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                    p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                    // Set up the discrete source
                    auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                    QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                            GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                    QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                            Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                    p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                    p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                    p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                    // Set up the finite difference solver for oxygen (which handles everything)
                    auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                    p_oxygen_solver->SetPde(p_oxygen_pde);
                    p_oxygen_solver->SetLabel("oxygen");
                    p_oxygen_solver->SetGrid(p_grid);

                    // Set the haematocrit solver
                    std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                    if (h_solver==1)
                    {
                        std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
                        p_abstract_haematocrit_solver = p_haematocrit_solver;       
                    }
                    else if (h_solver==2)
                    {
                        std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                    
                    }
                    else if (h_solver==3)
                    {
                        std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
                        p_abstract_haematocrit_solver = p_haematocrit_solver;      
                    }
                    else if (h_solver==4)
                    {
                        std::cout << "Now using FungHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
                        p_abstract_haematocrit_solver = p_haematocrit_solver;      
                    }

                    // Set up the viscosity calculator
                    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                    p_viscosity_calculator->SetVesselNetwork(p_network);
                    p_viscosity_calculator->Calculate();
                    
                    // Set up the impedance calculator
                    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                    p_impedance_calculator->SetVesselNetwork(p_network);
                    p_impedance_calculator->Calculate();

                    // Set up the flow solver
                    FlowSolver<2> flow_solver;
                    flow_solver.SetVesselNetwork(p_network);
                    flow_solver.SetUp();

                    // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance. The loop runs the solvers, calculates the max. difference in the network between the previous H of a vessel and the new H in the vessel, and repeats this until all vessels converge on a value within a certain tolerance.
                    unsigned max_iter = 1000;  // simulation fails if it doesn't converge in these many iterations
                    double tolerance2 = 1.e-10;
                    std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                    std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                    for(unsigned iteration=0;iteration<max_iter;iteration++)
                    {
                        // Run the solvers
                        p_impedance_calculator->Calculate();
                        flow_solver.SetUp();
                        flow_solver.Solve();
                        p_abstract_haematocrit_solver->Calculate();
                        p_viscosity_calculator->Calculate();

                        // Check for convergence 
                        double max_difference = 0.0;
                        double h_for_max = 0.0;
                        double prev_for_max = 0.0;
                        for(unsigned segment_index=0;segment_index<segments.size();segment_index++)  // for all the segments in the network
                        {
                            // // Set segments with no flow to have no haematocrit
                            // if (fabs(segments[segment_index]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                            // {
                            //     segments[segment_index]->GetFlowProperties()->SetHaematocrit(0.0);
                            // }

                            double current_haematocrit = segments[segment_index]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                            double difference = std::abs(current_haematocrit - previous_haematocrit[segment_index]);  // difference in new H in segment and old H
                            
                            // Log the max. difference calculated
                            if(difference>max_difference)
                            {
                                max_difference = difference;
                                h_for_max = current_haematocrit;
                                prev_for_max = previous_haematocrit[segment_index];
                            }
                            previous_haematocrit[segment_index] = current_haematocrit;
                        }

                        // Print the max. difference of the iteration
                        std::cout << "The maximum difference calculated in iteration " << iteration << " is " << h_for_max << " - " << prev_for_max << " = " << max_difference << std::endl;

                        if(max_difference<=tolerance2)  
                        {
                            std::cout << "Converged after: " << iteration << " iterations." <<  std::endl;
                            break;
                        }
                        else
                        {
                            if(iteration%1==0)
                            {
                                // std::cout << "Max Difference at iter: " << iteration << " is " << max_difference << std::endl;
                                std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(iteration) + ".vtp";
                                std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                                p_network->Write(output_file);
                            }
                        }
                        if(iteration==max_iter-1)  // if there is no convergence after all the iterations, print the error message
                        {
                            std::cout << "Problem encountered in Dichotomous Network with h_solver = " << h_solver << " using alpha = " << n_alpha << " and lambda = " << lambda << std::endl;
                            EXCEPTION("Did not converge after " + std::to_string(iteration) + " iterations.");
                        }
                    }
                
                    // Run the simulation 
                    SimulationTime::Instance()->SetStartTime(0.0);
                    SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                    auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                    p_microvessel_solver->SetVesselNetwork(p_network);
                    p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                    p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                    p_microvessel_solver->Run();

                    // Print the average oxygenation
                    std::vector<double> solution = p_oxygen_solver->GetSolution();
                    double average_oxygen = 0.0;
                    for(unsigned jdx=0;jdx<solution.size();jdx++)
                    {
                        average_oxygen += solution[jdx];
                    }
                    average_oxygen /= double(solution.size());
                    std::cout << "Average oxygen: " << average_oxygen << std::endl;
                    
                    // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                    std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                    p_network->Write(output_file);

                    // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                    ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                    ParameterCollection::Instance()->Destroy();
                    BaseUnits::Instance()->Destroy();
                    SimulationTime::Instance()->Destroy();
                }
            }
        }
    }

    // Make a multi-generation forking network with different h-splitting rules and individual pruning
    void xTestDichotomousNetworkWithIndividualPruningAndFlow2DAndVaryingMeansPaper1()
        {
            // Initialise error log
            std::ostringstream error_log;
            error_log << "\n The following simulations failed to converge: \n";

            // Set network in filename
            std::string network_name = "TestDichotomousNetwork/";

            // Create the output file for the PQs
            std::ofstream outfile;
            outfile.open("/tmp/narain/testoutput/TestDichotomousNetwork/forking_nonrandom_individual_pruning_perfusion_quotients.txt", std::ios_base::app);
            outfile.close();

            // Define the key pruning parameters
            // unsigned n_vessels = 250;  // number of non-inlet/outlet vessels from which to select ones to kill
            unsigned n_vessels = 508;  // number of non-inlet/outlet vessels from which to select ones to kill
            // unsigned n_vessels = 124;  // number of non-inlet/outlet vessels from which to select ones to kill
            // double percToKill = 0.2;  // percentage of vessels to kill
            // unsigned ToBeKilled = (unsigned)(percToKill*n_vessels);  // number to kill
            unsigned ToBeKilled = n_vessels;  // number to kill

            // Run the simulation with different solvers of interest
            for (unsigned h_solver=1; h_solver<=1; h_solver++)
            {   
                // Generate the network for various lambdas
                for (unsigned k_aux=1; k_aux<2; k_aux++)  // generate the network for lambda = 4
                {
                    // Set key vessel parameters
                    double dimless_length = 1.0;  

                    // Non-dimensionalising the length
                    for(unsigned i_aux=1; i_aux<order+1; i_aux++)
                    {
                        dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
                    }
    
                    // Set up the reference length for the simulation
                    QLength reference_length(1.0_um);
                    BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

                    // Set input radius
                    QLength input_radius(8.5 *GenericParameters::mpCapillaryRadius->GetValue());  // = *5_um (*10_um for diameter)
                    // VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);

                    // Set the viscosity
                    QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
                    // QDynamicViscosity viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();  // = 0.0012

                    // Set the inlet and initial haematocrit
                    double initial_haematocrit = 0.45;
                    // double initial_haematocrit = 0.36;
                    // double initial_haematocrit = 0.54;

                    // Set threshold for perfusion quotient
                    QFlowRate threshold = 3.e-12*unit::metre_cubed_per_second;  // set flow threshold for what counts as a perfused vessel

                    // Generate the networks
                    VesselNetworkGenerator<2> network_generator;

                    // Set lambda
                    double lambda;  // lambda = length/diameter
                    double twicelambda;  // used as an input parameter
                    lambda = 2.0+double(k_aux)*2.0;
                    twicelambda = 2.0*lambda;

                    // Set the number of different means
                    unsigned max_means = 3;

                    // Run the simulation with different heterogeneities
                    for (unsigned n_alpha=0; n_alpha<=max_alpha; n_alpha++)
                    { 
                        // Offset the radii depending on the mean and SD
                        for (unsigned n_mean=0; n_mean<max_means; n_mean++)
                        { 
                            // Set up flag for broken solver
                            unsigned broken_solver = 0;

                            // Set alpha
                            double alpha = 1.0+(double)n_alpha*0.1;  // alpha determines the relative radius of the left vessel 

                            // Height of of first-order vessels
                            QLength main_vert_length = 0.9*twicelambda*input_radius*pow(2.0,-1.0/3.0);
                            
                            // Height of the domain
                            // QLength domain_side_length_y = 4.0*main_vert_length;

                            // Length of the domain
                            QLength domain_side_length_x = dimless_length*2.0*twicelambda*input_radius;
                    
                            // Generate the network
                            std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapy(order, main_vert_length, input_radius, alpha, twicelambda);

                            // Identify input and output nodes and assign them properties
                            VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(0.0_um,2.0*main_vert_length));
                            VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(domain_side_length_x, 2.0*main_vert_length));
                            p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
                            // p_network->GetNode(0)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));  // = 3333.06
                            p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
                            p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
                            // p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));  // = 1999.84
                            p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);

                            // Specify the minimum offsets for all alphas (100 um, 6 gens)
                            // QLength min_array[] = {-9.95_um, -9.61_um, -8.75_um, -7.56_um, -6.21_um};
                            // QLength mean_array[] = {-4.21_um, -3.87_um, -3.01_um, -1.82_um, -0.47_um};
                            // QLength max_array[] = {0.93_um, 1.27_um, 2.13_um, 3.32_um, 4.67_um};

                            // Specify the minimum offsets for all alphas (75 um, 6 gens)
                            QLength min_array[] = {-1.77_um, -1.52_um, -0.87_um, 0.02_um, 1.04_um};
                            QLength mean_array[] = {3.97_um, 4.22_um, 4.87_um, 5.76_um, 6.78_um};
                            QLength max_array[] = {9.11_um, 9.36_um, 10.01_um, 10.9_um, 11.92_um};

                            // Specify the minimum offsets for all alphas (50 um, 6 gens)
                            // QLength min_array[] = {6.41_um, 6.58_um, 7.01_um, 7.6_um, 8.28_um};
                            // QLength mean_array[] = {12.15_um, 12.32_um, 12.75_um, 13.34_um, 14.02_um};
                            // QLength max_array[] = {17.29_um, 17.46_um, 17.89_um, 18.48_um, 19.16_um};
                        
                            // Specify the minimum offsets for all alphas (75 um, 5 gens)
                            // QLength min_array[] = {-7.64_um, -7.39_um, -6.74_um, -5.84_um, -4.8_um};
                            // QLength mean_array[] = {-1.9_um, -1.65_um, -1.0_um, -0.1_um, 0.94_um};
                            // QLength max_array[] = {3.24_um, 3.49_um, 4.14_um, 5.04_um, 6.08_um};

                            // Specify the minimum offsets for all alphas (50 um, 5 gens)
                            // QLength min_array[] = {2.49_um, 2.66_um, 3.09_um, 3.69_um, 4.38_um};
                            // QLength mean_array[] = {8.23_um, 8.4_um, 8.83_um, 9.43_um, 10.12_um};
                            // QLength max_array[] = {13.37_um, 13.54_um, 13.97_um, 14.57_um, 15.26_um};

                            // Initialise the offsets
                            std::vector<QLength> min_list(min_array, min_array+5);
                            std::vector<QLength> mean_list(mean_array, mean_array+5);
                            std::vector<QLength> max_list(max_array, max_array+5);

                            // Set the radius offset based on the alpha
                            std::vector<std::shared_ptr<Vessel<2> > > vessels;
                            vessels = p_network->GetVessels();
                            if (n_mean==0)
                                {
                                    QLength radius_offset = min_list[n_alpha];
                                    for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)
                                    {                            
                                        // Get the current segment's radius
                                        QLength current_radius = vessels[vessel_index]->GetRadius();

                                        // Calculate the new radius
                                        QLength new_vessel_radius = current_radius + (radius_offset*0.5);

                                        // Set the new radius
                                        vessels[vessel_index]->SetRadius(new_vessel_radius);
                                    }
                                }
                            else if (n_mean==1)
                            {
                                QLength radius_offset = mean_list[n_alpha];
                                for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)
                                {                            
                                    // Get the current segment's radius
                                    QLength current_radius = vessels[vessel_index]->GetRadius();

                                    // Calculate the new radius
                                    QLength new_vessel_radius = current_radius + (radius_offset*0.5);
    
                                    // Set the new radius
                                    vessels[vessel_index]->SetRadius(new_vessel_radius);
                                }
                            }
                            else if (n_mean==2)
                            {
                                QLength radius_offset = max_list[n_alpha];
                                for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)
                                {                            
                                    // Get the current segment's radius
                                    QLength current_radius = vessels[vessel_index]->GetRadius();

                                    // Calculate the new radius
                                    QLength new_vessel_radius = current_radius + (radius_offset*0.5);

                                    // Set the new radius
                                    vessels[vessel_index]->SetRadius(new_vessel_radius);
                                }
                            }
                            vessels = p_network->GetVessels();

                            // Set the h-solver
                            std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                            std::string solver_name;
                            if (h_solver==1)
                            {
                                solver_name = "ConstantHaematocrit/"; 
                                std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                                auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                                p_haematocrit_solver->SetVesselNetwork(p_network);
                                p_haematocrit_solver->SetHaematocrit(initial_haematocrit);
                                p_abstract_haematocrit_solver = p_haematocrit_solver;     
                            }
                                                    
                            // Store the lambda values for file name
                            std::stringstream lambda_stream;
                            lambda_stream << std::fixed << std::setprecision(0) << lambda;
                            std::string lambda_string = lambda_stream.str();
                            std::string lambda_directory = network_name + solver_name + "/Lambda" + lambda_string; 

                            // Store the n_mean values for file name
                            std::stringstream mean_stream;
                            mean_stream << std::fixed << std::setprecision(0) << n_mean;
                            std::string mean_string = mean_stream.str();

                            // Set up the viscosity calculator
                            auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                            p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                            p_viscosity_calculator->SetVesselNetwork(p_network);
                            p_viscosity_calculator->Calculate();
                            
                            // Set up the impedance calculator
                            auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                            p_impedance_calculator->SetVesselNetwork(p_network);
                            p_impedance_calculator->Calculate();

                            // Set up the flow solver
                            FlowSolver<2> flow_solver;
                            flow_solver.SetVesselNetwork(p_network);
                            flow_solver.SetUp();
        
                            // Prune all vessels up to specified dose 
                            for(unsigned KilledVessels=0; KilledVessels < ToBeKilled+1; KilledVessels++)
                            { 
                                // Display status message
                                std::cout << "Now killed " << KilledVessels << " vessels." << std::endl;  

                                // Set filename
                                std::stringstream alpha_stream;
                                std::stringstream kill_stream;
                                alpha_stream << std::fixed << std::setprecision(2) << alpha;
                                kill_stream << std::fixed << std::setprecision(0) << KilledVessels;
                                std::string alpha_string = alpha_stream.str();                            
                                std::string kill_string = kill_stream.str();
                                std::string file_name = lambda_directory + "/Alpha" + alpha_string + "/Mean" + mean_string + "/Kills" + kill_string;
                                std::string str_directory_name = file_name;
                                auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                                // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance. The loop runs the solvers, calculates the max. difference in the network between the previous H of a vessel and the new H in the vessel, and repeats this until all vessels converge on a value within a certain tolerance.
                                std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                                // Run the solvers
                                p_impedance_calculator->Calculate();
                                flow_solver.SetUp();
                                flow_solver.Solve();
                                p_abstract_haematocrit_solver->Calculate();
                                p_viscosity_calculator->Calculate();

                                // Check for convergence 
                                for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)  // for all the segments 
                                {
                                    // Set segments with no flow to have no haematocrit
                                    if (fabs(segments[segment_index]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                    {
                                        segments[segment_index]->GetFlowProperties()->SetHaematocrit(0.0);
                                    }                                    
                                }

                                // If solver doesn't converge, move on to next one
                                if (broken_solver == 1)
                                {
                                    continue;
                                }

                                // Write the PQs
                                outfile.open("/tmp/narain/testoutput/TestDichotomousNetwork/forking_nonrandom_individual_pruning_perfusion_quotients.txt", std::ios_base::app);
                                outfile << network_name << " " << solver_name << " " << lambda_string << " " << alpha_string << " " << mean_string << " " << KilledVessels << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                                outfile.close();
                                
                                // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                                std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                                p_network->Write(output_file);

                                // Remove the smallest vessel
                                QLength minimum_radius = input_radius;
                                unsigned int minimum_index = 0;
                                std::vector<std::shared_ptr<Vessel<2> > > vessels = p_network->GetVessels();
                                for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)  // for all the segments in the network
                                {
                                    // Get the current segment's radius
                                    QLength current_radius = vessels[vessel_index]->GetRadius();

                                    // If the current radius is less than the minimum radius, record the new minimum
                                    if (current_radius < minimum_radius)
                                    {
                                        minimum_radius = current_radius;
                                        minimum_index = vessel_index;
                                    }
                                }
                                p_network->RemoveVessel(vessels[minimum_index], true);  // remove the vessel

                                // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                                ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                                ParameterCollection::Instance()->Destroy();
                                BaseUnits::Instance()->Destroy();
                                SimulationTime::Instance()->Destroy();
                            }
                        }
                    }
                }
            }
            // Print the error log
            std::string error_message = error_log.str();
            std::cout << error_message << std::endl; 
        }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Current Forking Networks
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Make a multi-generation forking network with different h-splitting rules and individual pruning
    void xTestDichotomousNetworkWithIndividualPruningAndFlow2DAndVaryingMeans()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";

        // Set network in filename
        std::string network_name = "TestDichotomousNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("/tmp/narain/testoutput/TestDichotomousNetwork/forking_nonrandom_individual_pruning_perfusion_quotients.txt", std::ios_base::app);
        outfile.close();

        // Define the key pruning parameters
        unsigned n_vessels = 250;  // number of non-inlet/outlet vessels from which to select ones to kill
        // double percToKill = 0.2;  // percentage of vessels to kill
        // unsigned ToBeKilled = (unsigned)(percToKill*n_vessels);  // number to kill
        unsigned ToBeKilled = n_vessels;  // number to kill

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=1; h_solver++)
        {   
            // Generate the network for various lambdas
            for (unsigned k_aux=1; k_aux<2; k_aux++)  // generate the network for lambda = 4
            {
                // Set key vessel parameters
                double dimless_length = 1.0;  

                // ??? Non-dimensionalising the length
                for(unsigned i_aux=1; i_aux<order+1; i_aux++)
                {
                    dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
                }
   
                // Set up the reference length for the simulation
                QLength reference_length(1.0_um);
                BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

                // Set input radius
                QLength input_radius(7.5 *GenericParameters::mpCapillaryRadius->GetValue());  // = 37.5_um
                // VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);

                // Set the viscosity
                QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
                // QDynamicViscosity viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();  // = 0.0012

                // Set the inlet and initial haematocrit
                double initial_haematocrit = 0.45;

                // Set threshold for perfusion quotient
		        QFlowRate threshold = 3.e-12*unit::metre_cubed_per_second;  // set flow threshold for what counts as a perfused vessel

                // Generate the networks
                VesselNetworkGenerator<2> network_generator;

                // Set lambda
                double lambda;  // lambda = length/diameter
                double twicelambda;  // used as an input parameter
                lambda = 2.0+double(k_aux)*2.0;
                twicelambda = 2.0*lambda;

                // Set the number of different means
                unsigned max_means = 3;

                // Run the simulation with different heterogeneities
                for (unsigned n_alpha=0; n_alpha<=max_alpha; n_alpha++)
                { 
                    // Offset the radii depending on the mean and SD
                    for (unsigned n_mean=0; n_mean<max_means; n_mean++)
                    { 
                        // Set up flag for broken solver
                        unsigned broken_solver = 0;

                        // Set alpha
                        double alpha = 1.0+(double)n_alpha*0.1;  // alpha determines the relative radius of the left vessel 

                        // Height of of first-order vessels
                        QLength main_vert_length = 0.9*twicelambda*input_radius*pow(2.0,-1.0/3.0);
                        
                        // Height of the domain
                        QLength domain_side_length_y = 4.0*main_vert_length+10.0_um;  // 10 um allows for enough room on top and bottom of network (considering the grid size is 10 um)
                        std::cout << domain_side_length_y << std::endl;

                        // Length of the domain
                        QLength domain_side_length_x = dimless_length*2.0*twicelambda*input_radius;
                
                        // Generate the network
                        std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapy(order, main_vert_length, input_radius, alpha, twicelambda);

                        // Identify input and output nodes and assign them properties
                        VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(0.0_um,2.0*main_vert_length));
                        VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(domain_side_length_x, 2.0*main_vert_length));
                        p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
                        // p_network->GetNode(0)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));  // = 3333.06
                        p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
                        p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
                        // p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));  // = 1999.84
                        p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);

                        // Specify the minimum offsets for all alphas (100 um, 6 gens)
                        // QLength min_array[] = {-9.95_um, -9.61_um, -8.75_um, -7.56_um, -6.21_um};
                        // QLength mean_array[] = {-4.21_um, -3.87_um, -3.01_um, -1.82_um, -0.47_um};
                        // QLength max_array[] = {0.93_um, 1.27_um, 2.13_um, 3.32_um, 4.67_um};

                        // Specify the minimum offsets for all alphas (75 um, 6 gens)
                        QLength min_array[] = {-1.77_um, -1.52_um, -0.87_um, 0.02_um, 1.04_um};
                        QLength mean_array[] = {3.97_um, 4.22_um, 4.87_um, 5.76_um, 6.78_um};
                        QLength max_array[] = {9.11_um, 9.36_um, 10.01_um, 10.9_um, 11.92_um};

                        // Specify the minimum offsets for all alphas (50 um, 6 gens)
                        // QLength min_array[] = {6.41_um, 6.58_um, 7.01_um, 7.6_um, 8.28_um};
                        // QLength mean_array[] = {12.15_um, 12.32_um, 12.75_um, 13.34_um, 14.02_um};
                        // QLength max_array[] = {17.29_um, 17.46_um, 17.89_um, 18.48_um, 19.16_um};
                    
                        // Specify the minimum offsets for all alphas (75 um, 5 gens)
                        // QLength min_array[] = {-7.64_um, -7.39_um, -6.74_um, -5.84_um, -4.8_um};
                        // QLength mean_array[] = {-1.9_um, -1.65_um, -1.0_um, -0.1_um, 0.94_um};
                        // QLength max_array[] = {3.24_um, 3.49_um, 4.14_um, 5.04_um, 6.08_um};

                        // Specify the minimum offsets for all alphas (50 um, 5 gens)
                        // QLength min_array[] = {2.49_um, 2.66_um, 3.09_um, 3.69_um, 4.38_um};
                        // QLength mean_array[] = {8.23_um, 8.4_um, 8.83_um, 9.43_um, 10.12_um};
                        // QLength max_array[] = {13.37_um, 13.54_um, 13.97_um, 14.57_um, 15.26_um};

                        // Initialise the offsets
                        std::vector<QLength> min_list(min_array, min_array+5);
                        std::vector<QLength> mean_list(mean_array, mean_array+5);
                        std::vector<QLength> max_list(max_array, max_array+5);

                        // Set the radius offset based on the alpha
                        std::vector<std::shared_ptr<Vessel<2> > > vessels;
                        vessels = p_network->GetVessels();
                        if (n_mean==0)
                            {
                                QLength radius_offset = min_list[n_alpha];
                                for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)
                                {                            
                                    // Get the current segment's radius
                                    QLength current_radius = vessels[vessel_index]->GetRadius();

                                    // Calculate the new radius
                                    QLength new_vessel_radius = current_radius + (radius_offset*0.5);

                                    // Set the new radius
                                    vessels[vessel_index]->SetRadius(new_vessel_radius);
                                }
                            }
                        else if (n_mean==1)
                        {
                            QLength radius_offset = mean_list[n_alpha];
                            for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)
                            {                            
                                // Get the current segment's radius
                                QLength current_radius = vessels[vessel_index]->GetRadius();

                                // Calculate the new radius
                                QLength new_vessel_radius = current_radius + (radius_offset*0.5);
 
                                // Set the new radius
                                vessels[vessel_index]->SetRadius(new_vessel_radius);
                            }
                        }
                        else if (n_mean==2)
                        {
                            QLength radius_offset = max_list[n_alpha];
                            for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)
                            {                            
                                // Get the current segment's radius
                                QLength current_radius = vessels[vessel_index]->GetRadius();

                                // Calculate the new radius
                                QLength new_vessel_radius = current_radius + (radius_offset*0.5);

                                // Set the new radius
                                vessels[vessel_index]->SetRadius(new_vessel_radius);
                            }
                        }
                        vessels = p_network->GetVessels();

                        // Set the h-solver
                        std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                        std::string solver_name;
                        if (h_solver==1)
                        {
                            solver_name = "ConstantHaematocrit/"; 
                            std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                            auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                            p_haematocrit_solver->SetVesselNetwork(p_network);
                            p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                            p_abstract_haematocrit_solver = p_haematocrit_solver;     
                        }
                        else if (h_solver==2)
                        {
                            solver_name = "PriesHaematocrit/"; 
                            std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                            auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                            // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                            p_haematocrit_solver->SetVesselNetwork(p_network);
                            p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                            p_abstract_haematocrit_solver = p_haematocrit_solver;                
                        }
                        else if (h_solver==3)
                        {
                            solver_name = "MemoryHaematocrit/"; 
                            std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                            auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                            // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                            p_haematocrit_solver->SetVesselNetwork(p_network);
                            p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                            p_abstract_haematocrit_solver = p_haematocrit_solver;     
                        }
                        else if (h_solver==4)
                        {
                            solver_name = "FungHaematocrit/"; 
                            std::cout << "Now using FungHaematocritSolver..." << std::endl;
                            auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                            // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                            p_haematocrit_solver->SetVesselNetwork(p_network);
                            p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                            p_abstract_haematocrit_solver = p_haematocrit_solver;     
                        }
                                                
                        // Store the lambda values for file name
                        std::stringstream lambda_stream;
                        lambda_stream << std::fixed << std::setprecision(0) << lambda;
                        std::string lambda_string = lambda_stream.str();
                        std::string lambda_directory = network_name + solver_name + "/Lambda" + lambda_string; 

                        // Store the n_mean values for file name
                        std::stringstream mean_stream;
                        mean_stream << std::fixed << std::setprecision(0) << n_mean;
                        std::string mean_string = mean_stream.str();
                        // std::string lambda_directory = network_name + solver_name + "/Lambda" + lambda_string; 

                        // Set up the grid for the finite difference solver
                        auto p_grid = RegularGrid<2>::Create();
                        // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
                        // QLength grid_spacing(1.0*1_um);  // the simulation time gets quite long if you reduce the resolution further
                        // QLength grid_spacing = 10_um;
                        p_grid->SetSpacing(grid_spacing);
                        c_vector<unsigned, 3> dimensions;
                        dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                        dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                        dimensions[2] = 1;
                        p_grid->SetDimensions(dimensions);

                        // // Choose the PDE
                        auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                        
                        // // Set the diffusivity and decay terms
                        p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                        p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                        // Set up the discrete source
                        auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                        QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                                GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                        QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                                Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                        p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                        p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                        p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                        // // Set up the finite difference solver for oxygen (which handles everything)
                        auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                        p_oxygen_solver->SetPde(p_oxygen_pde);
                        p_oxygen_solver->SetLabel("oxygen");
                        p_oxygen_solver->SetGrid(p_grid);

                        // Set up the viscosity calculator
                        auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                        p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                        p_viscosity_calculator->SetVesselNetwork(p_network);
                        p_viscosity_calculator->Calculate();
                        
                        // Set up the impedance calculator
                        auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                        p_impedance_calculator->SetVesselNetwork(p_network);
                        p_impedance_calculator->Calculate();

                        // Set up the flow solver
                        FlowSolver<2> flow_solver;
                        flow_solver.SetVesselNetwork(p_network);
                        flow_solver.SetUp();
    
                        // Prune all vessels up to specified dose 
                        for(unsigned KilledVessels=0; KilledVessels < ToBeKilled+1; KilledVessels++)
                        { 
                            // Display status message
                            std::cout << "Now killed " << KilledVessels << " vessels." << std::endl;  

                            // Set filename
                            std::stringstream alpha_stream;
                            std::stringstream kill_stream;
                            alpha_stream << std::fixed << std::setprecision(2) << alpha;
                            kill_stream << std::fixed << std::setprecision(0) << KilledVessels;
                            std::string alpha_string = alpha_stream.str();                            
                            std::string kill_string = kill_stream.str();
                            std::string file_name = lambda_directory + "/Alpha" + alpha_string + "/Mean" + mean_string + "/Kills" + kill_string;
                            std::string str_directory_name = file_name;
                            auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                            // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance. The loop runs the solvers, calculates the max. difference in the network between the previous H of a vessel and the new H in the vessel, and repeats this until all vessels converge on a value within a certain tolerance.
                            unsigned max_iter = 100;  // simulation fails if it doesn't converge in these many iterations
                            double tolerance2 = 0.01; //1.e-10
                            std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                            std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                            for (unsigned idx=0;idx<max_iter;idx++)
                            {
                                // Run the solvers
                                p_impedance_calculator->Calculate();
                                flow_solver.SetUp();
                                flow_solver.Solve();
                                p_abstract_haematocrit_solver->Calculate();
                                p_viscosity_calculator->Calculate();

                                // Check for convergence 
                                double max_difference = 0.0;
                                double h_for_max = 0.0;
                                double prev_for_max = 0.0;
                                for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)  // for all the segments 
                                {
                                    // Set segments with no flow to have no haematocrit
                                    if (fabs(segments[segment_index]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                    {
                                        segments[segment_index]->GetFlowProperties()->SetHaematocrit(0.0);
                                    }                                    
                                    double current_haematocrit = segments[segment_index]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                                    double difference = std::abs(current_haematocrit - previous_haematocrit[segment_index]);  // difference in new H in segment and old H
                                    if(difference>max_difference)  // if the difference is greater than previous max.
                                    {
                                        max_difference = difference;
                                        h_for_max = current_haematocrit;
                                        prev_for_max = previous_haematocrit[segment_index];
                                    }
                                    previous_haematocrit[segment_index] = current_haematocrit;
                                }
                                std::cout << "Segment H. at max difference: " << h_for_max << ", Prev. segment H at max difference:" << prev_for_max << std::endl;
                                
                                if(max_difference<=tolerance2)  
                                {
                                    std::cout << "Converged after: " << idx << " iterations." <<  std::endl;
                                    broken_solver = 0;
                                    break;
                                }
                                else
                                {
                                    if(idx%1==0)
                                    {
                                        std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                                        // std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                                        // std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                                        // p_network->Write(output_file);
                                    }
                                }
                                
                                if(idx==max_iter-1)  // if there is no convergence after all the iterations, print the error message
                                {
                                    std::cout << "Problem encountered in Dichotomous Network with h_solver = " << h_solver << " using n_alpha = " << n_alpha << " and lambda = " << lambda << std::endl;
                                    // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                                    error_log << "\n Problem encountered in Dichotomous Network with h_solver = " << h_solver << " using n_alpha = " << n_alpha << " and lambda = " << lambda; 
                                    broken_solver = 1;
                                    break;
                                }
                            }

                            // If solver doesn't converge, move on to next one
                            if (broken_solver == 1)
                            {
                                continue;
                            }

                            // Run the simulation 
                            SimulationTime::Instance()->SetStartTime(0.0);
                            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                            auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                            p_microvessel_solver->SetVesselNetwork(p_network);
                            p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                            p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                            p_microvessel_solver->Run();

                            // Print the average oxygenation
                            std::vector<double> solution = p_oxygen_solver->GetSolution();
                            double average_oxygen = 0.0;
                            for(unsigned jdx=0;jdx<solution.size();jdx++)
                            {
                                average_oxygen += solution[jdx];
                            }
                            average_oxygen /= double(solution.size());
                            std::cout << "Average oxygen: " << average_oxygen << std::endl;
                            // }
                            // Write the PQs
                            outfile.open("/tmp/narain/testoutput/TestDichotomousNetwork/forking_nonrandom_individual_pruning_perfusion_quotients.txt", std::ios_base::app);
                            outfile << network_name << " " << solver_name << " " << lambda_string << " " << alpha_string << " " << mean_string << " " << KilledVessels << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                            outfile.close();
                            
                            // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                            std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                            p_network->Write(output_file);

                            // Remove the smallest vessel
                            QLength minimum_radius = input_radius;
                            unsigned int minimum_index = 0;
                            std::vector<std::shared_ptr<Vessel<2> > > vessels = p_network->GetVessels();
                            for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)  // for all the segments in the network
                            {
                                // Get the current segment's radius
                                QLength current_radius = vessels[vessel_index]->GetRadius();

                                // If the current radius is less than the minimum radius, record the new minimum
                                if (current_radius < minimum_radius)
                                {
                                    minimum_radius = current_radius;
                                    minimum_index = vessel_index;
                                }
                            }
                            p_network->RemoveVessel(vessels[minimum_index], true);  // remove the vessel

                            // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                            ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                            ParameterCollection::Instance()->Destroy();
                            BaseUnits::Instance()->Destroy();
                            SimulationTime::Instance()->Destroy();
                            // }
                        }
                    }
                }
            }
        }
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a multi-generation forking network with different h-splitting rules and radius thereshold pruning
    void xTestDichotomousNetworkWithThresholdPruningAndFlow2DAndVaryingMeans()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";

        // Set network in filename
        std::string network_name = "TestDichotomousNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("/tmp/narain/testoutput/TestDichotomousNetwork/forking_nonrandom_threshold_pruning_perfusion_quotients.txt", std::ios_base::app);
        outfile.close();

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=4; h_solver++)
        {   
            // Generate the network for various lambdas
            for (unsigned k_aux=1; k_aux<2; k_aux++)  // generate the network for lambda = 4
            {
                // Set key vessel parameters
                double dimless_length = 1.0;  

                // ??? Non-dimensionalising the length
                for(unsigned i_aux=1; i_aux<order+1; i_aux++)
                {
                    dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
                }
   
                // Set up the reference length for the simulation
                QLength reference_length(1.0_um);
                BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

                // Set input radius
                QLength input_radius(7.5 *GenericParameters::mpCapillaryRadius->GetValue());  // = 37.5_um
                // VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);

                // Set the viscosity
                QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
                // QDynamicViscosity viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();  // = 0.0012

                // Set the inlet and initial haematocrit
                double initial_haematocrit = 0.45;

                // Set threshold for perfusion quotient
		        QFlowRate threshold = 3.e-12*unit::metre_cubed_per_second;  // set flow threshold for what counts as a perfused vessel

                // Generate the networks
                VesselNetworkGenerator<2> network_generator;

                // Set lambda
                double lambda;  // lambda = length/diameter
                double twicelambda;  // used as an input parameter
                lambda = 2.0+double(k_aux)*2.0;
                twicelambda = 2.0*lambda;

                // Set the number of different means
                unsigned max_means = 3;

                // Run the simulation with different heterogeneities
                for (unsigned n_alpha=0; n_alpha<=max_alpha; n_alpha++)
                { 
                    // Offset the radii depending on the mean and SD
                    for (unsigned n_mean=0; n_mean<max_means; n_mean++)
                    { 
                        // Set up flag for broken solver
                        unsigned broken_solver = 0;

                        // Set alpha
                        double alpha = 1.0+(double)n_alpha*0.1;  // alpha determines the relative radius of the left vessel 

                        // Height of of first-order vessels
                        QLength main_vert_length = 0.9*twicelambda*input_radius*pow(2.0,-1.0/3.0);
                        
                        // Height of the domain
                        QLength domain_side_length_y = 4.0*main_vert_length+10.0_um;  // 10 um allows for enough room on top and bottom of network (considering the grid size is 10 um)

                        // Length of the domain
                        QLength domain_side_length_x = dimless_length*2.0*twicelambda*input_radius;
                
                        // Generate the network
                        std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapy(order, main_vert_length, input_radius, alpha, twicelambda);

                        // Identify input and output nodes and assign them properties
                        VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(0.0_um,2.0*main_vert_length));
                        VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(domain_side_length_x, 2.0*main_vert_length));
                        p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
                        // p_network->GetNode(0)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));  // = 3333.06
                        p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
                        p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
                        // p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));  // = 1999.84
                        p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);

                        // Specify the minimum offsets for all alphas (100 um, 6 gens)
                        // QLength min_array[] = {-9.95_um, -9.61_um, -8.75_um, -7.56_um, -6.21_um};
                        // QLength mean_array[] = {-4.21_um, -3.87_um, -3.01_um, -1.82_um, -0.47_um};
                        // QLength max_array[] = {0.93_um, 1.27_um, 2.13_um, 3.32_um, 4.67_um};

                        // Specify the minimum offsets for all alphas (75 um, 6 gens)
                        QLength min_array[] = {-1.77_um, -1.52_um, -0.87_um, 0.02_um, 1.04_um};
                        QLength mean_array[] = {3.97_um, 4.22_um, 4.87_um, 5.76_um, 6.78_um};
                        QLength max_array[] = {9.11_um, 9.36_um, 10.01_um, 10.9_um, 11.92_um};

                        // Specify the minimum offsets for all alphas (50 um, 6 gens)
                        // QLength min_array[] = {6.41_um, 6.58_um, 7.01_um, 7.6_um, 8.28_um};
                        // QLength mean_array[] = {12.15_um, 12.32_um, 12.75_um, 13.34_um, 14.02_um};
                        // QLength max_array[] = {17.29_um, 17.46_um, 17.89_um, 18.48_um, 19.16_um};
                    
                        // Specify the minimum offsets for all alphas (75 um, 5 gens)
                        // QLength min_array[] = {-7.64_um, -7.39_um, -6.74_um, -5.84_um, -4.8_um};
                        // QLength mean_array[] = {-1.9_um, -1.65_um, -1.0_um, -0.1_um, 0.94_um};
                        // QLength max_array[] = {3.24_um, 3.49_um, 4.14_um, 5.04_um, 6.08_um};

                        // Specify the minimum offsets for all alphas (50 um, 5 gens)
                        // QLength min_array[] = {2.49_um, 2.66_um, 3.09_um, 3.69_um, 4.38_um};
                        // QLength mean_array[] = {8.23_um, 8.4_um, 8.83_um, 9.43_um, 10.12_um};
                        // QLength max_array[] = {13.37_um, 13.54_um, 13.97_um, 14.57_um, 15.26_um};

                        // Specify no offsets for all alphas for the natural network
                        // QLength min_array[] = {0_um, 0_um, 0_um, 0_um, 0_um};
                        // QLength mean_array[] = {0_um, 0_um, 0_um, 0_um, 0_um};
                        // QLength max_array[] = {0_um, 0_um, 0_um, 0_um, 0_um};

                        // Initialise the offsets
                        std::vector<QLength> min_list(min_array, min_array+5);
                        std::vector<QLength> mean_list(mean_array, mean_array+5);
                        std::vector<QLength> max_list(max_array, max_array+5);

                        // Set the radius offset based on the alpha
                        std::vector<std::shared_ptr<Vessel<2> > > vessels;
                        vessels = p_network->GetVessels();
                        if (n_mean==0)
                            {
                                QLength radius_offset = min_list[n_alpha];
                                for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)
                                {                            
                                    // Get the current segment's radius
                                    QLength current_radius = vessels[vessel_index]->GetRadius();

                                    // Calculate the new radius
                                    QLength new_vessel_radius = current_radius + (radius_offset*0.5);

                                    // Set the new radius
                                    vessels[vessel_index]->SetRadius(new_vessel_radius);
                                }
                            }
                        else if (n_mean==1)
                        {
                            QLength radius_offset = mean_list[n_alpha];
                            for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)
                            {                            
                                // Get the current segment's radius
                                QLength current_radius = vessels[vessel_index]->GetRadius();

                                // Calculate the new radius
                                QLength new_vessel_radius = current_radius + (radius_offset*0.5);
 
                                // Set the new radius
                                vessels[vessel_index]->SetRadius(new_vessel_radius);
                            }
                        }
                        else if (n_mean==2)
                        {
                            QLength radius_offset = max_list[n_alpha];
                            for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)
                            {                            
                                // Get the current segment's radius
                                QLength current_radius = vessels[vessel_index]->GetRadius();

                                // Calculate the new radius
                                QLength new_vessel_radius = current_radius + (radius_offset*0.5);

                                // Set the new radius
                                vessels[vessel_index]->SetRadius(new_vessel_radius);
                            }
                        }
                        vessels = p_network->GetVessels();

                        // Set the h-solver
                        std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                        std::string solver_name;
                        if (h_solver==1)
                        {
                            solver_name = "ConstantHaematocrit/"; 
                            std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                            auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                            p_haematocrit_solver->SetVesselNetwork(p_network);
                            p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                            p_abstract_haematocrit_solver = p_haematocrit_solver;     
                        }
                        else if (h_solver==2)
                        {
                            solver_name = "PriesHaematocrit/"; 
                            std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                            auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                            // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                            p_haematocrit_solver->SetVesselNetwork(p_network);
                            p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                            p_abstract_haematocrit_solver = p_haematocrit_solver;                
                        }
                        else if (h_solver==3)
                        {
                            solver_name = "MemoryHaematocrit/"; 
                            std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                            auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                            // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                            p_haematocrit_solver->SetVesselNetwork(p_network);
                            p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                            p_abstract_haematocrit_solver = p_haematocrit_solver;     
                        }
                        else if (h_solver==4)
                        {
                            solver_name = "FungHaematocrit/"; 
                            std::cout << "Now using FungHaematocritSolver..." << std::endl;
                            auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                            // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                            p_haematocrit_solver->SetVesselNetwork(p_network);
                            p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                            p_abstract_haematocrit_solver = p_haematocrit_solver;     
                        }
                                                
                        // Store the lambda values for file name
                        std::stringstream lambda_stream;
                        lambda_stream << std::fixed << std::setprecision(0) << lambda;
                        std::string lambda_string = lambda_stream.str();
                        std::string lambda_directory = network_name + solver_name + "/Lambda" + lambda_string; 

                        // Store the n_mean values for file name
                        std::stringstream mean_stream;
                        mean_stream << std::fixed << std::setprecision(0) << n_mean;
                        std::string mean_string = mean_stream.str();
                        // std::string lambda_directory = network_name + solver_name + "/Lambda" + lambda_string; 

                        // Set up the grid for the finite difference solver
                        auto p_grid = RegularGrid<2>::Create();
                        // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
                        // QLength grid_spacing(1.0*1_um);  // the simulation time gets quite long if you reduce the resolution further
                        // QLength grid_spacing = 10_um;
                        p_grid->SetSpacing(grid_spacing);
                        c_vector<unsigned, 3> dimensions;
                        dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                        dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                        dimensions[2] = 1;
                        p_grid->SetDimensions(dimensions);

                        // Choose the PDE
                        auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                        
                        // Set the diffusivity and decay terms
                        p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                        p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                        // Set up the discrete source
                        auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                        QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                                GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                        QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                                Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                        p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                        p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                        p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                        // Set up the finite difference solver for oxygen (which handles everything)
                        auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                        p_oxygen_solver->SetPde(p_oxygen_pde);
                        p_oxygen_solver->SetLabel("oxygen");
                        p_oxygen_solver->SetGrid(p_grid);

                        // Set up the viscosity calculator
                        auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                        p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                        p_viscosity_calculator->SetVesselNetwork(p_network);
                        p_viscosity_calculator->Calculate();
                        
                        // Set up the impedance calculator
                        auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                        p_impedance_calculator->SetVesselNetwork(p_network);
                        p_impedance_calculator->Calculate();

                        // Set up the flow solver
                        FlowSolver<2> flow_solver;
                        flow_solver.SetVesselNetwork(p_network);
                        flow_solver.SetUp();

                        // Set up pruning parameters
                        QLength max_radius_to_kill = 25.02_um;  // set threshold up to which vessels should be pruned
                        QLength radius_step = 0.1_um;  // Use the same number of data points as individual pruning
                        // QLength radius_step = 1_um;
                        QLength current_radius_threshold = 0_um;
    
                        // Prune all vessels up to specified dose 
                        while (current_radius_threshold <= max_radius_to_kill)
                        { 
                            // Display status message
                            double current_radius_threshold_um = current_radius_threshold*1000000;  // in micrometres
                            std::cout << "Now pruning up to vessels with radius = " << current_radius_threshold_um << " um" << std::endl;
                            
                            // Remove vessels under radius threshold
                            p_network->RemoveThinVessels(current_radius_threshold, false);  

                            // Set filename
                            std::stringstream alpha_stream;
                            std::stringstream threshold_stream;
                            alpha_stream << std::fixed << std::setprecision(2) << alpha;
                            threshold_stream << std::fixed << std::setprecision(2) << current_radius_threshold_um;
                            std::string alpha_string = alpha_stream.str();
                            std::string threshold_string = threshold_stream.str();
                            std::string file_name = lambda_directory + "/Alpha" + alpha_string + "/Mean" + mean_string + "/RadiusThreshold" + threshold_string;
                            std::string str_directory_name = file_name;
                            auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                            // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance. The loop runs the solvers, calculates the max. difference in the network between the previous H of a vessel and the new H in the vessel, and repeats this until all vessels converge on a value within a certain tolerance.
                            unsigned max_iter = 100;  // simulation fails if it doesn't converge in these many iterations
                            double tolerance2 = 0.01; //1.e-10
                            std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                            std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                            for (unsigned idx=0;idx<max_iter;idx++)
                            {
                                // Run the solvers
                                p_impedance_calculator->Calculate();
                                flow_solver.SetUp();
                                flow_solver.Solve();
                                p_abstract_haematocrit_solver->Calculate();
                                p_viscosity_calculator->Calculate();

                                // Check for convergence 
                                double max_difference = 0.0;
                                double h_for_max = 0.0;
                                double prev_for_max = 0.0;
                                for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)  // for all the segments 
                                {
                                    // Set segments with no flow to have no haematocrit
                                    if (fabs(segments[segment_index]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                    {
                                        segments[segment_index]->GetFlowProperties()->SetHaematocrit(0.0);
                                    }                                    
                                    double current_haematocrit = segments[segment_index]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                                    double difference = std::abs(current_haematocrit - previous_haematocrit[segment_index]);  // difference in new H in segment and old H
                                    if(difference>max_difference)  // if the difference is greater than previous max.
                                    {
                                        max_difference = difference;
                                        h_for_max = current_haematocrit;
                                        prev_for_max = previous_haematocrit[segment_index];
                                    }
                                    previous_haematocrit[segment_index] = current_haematocrit;
                                }
                                std::cout << "Segment H. at max difference: " << h_for_max << ", Prev. segment H at max difference:" << prev_for_max << std::endl;
                                if(max_difference<=tolerance2)  
                                {
                                    std::cout << "Converged after: " << idx << " iterations." <<  std::endl;
                                    broken_solver = 0;
                                    break;
                                }
                                else
                                {
                                    if(idx%1==0)
                                    {
                                        std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                                        // std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                                        // std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                                        // p_network->Write(output_file);
                                    }
                                }
                                if(idx==max_iter-1)  // if there is no convergence after all the iterations, print the error message
                                {
                                    std::cout << "Problem encountered in Dichotomous Network with h_solver = " << h_solver << " using n_alpha = " << n_alpha << " and lambda = " << lambda << std::endl;
                                    // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                                    error_log << "\n Problem encountered in Dichotomous Network with h_solver = " << h_solver << " using n_alpha = " << n_alpha << " and lambda = " << lambda; 
                                    broken_solver = 1;
                                    break;
                                }
                            }

                            // If solver doesn't converge, move on to next one
                            if (broken_solver == 1)
                            {
                                continue;
                            }

                            // Run the simulation 
                            SimulationTime::Instance()->SetStartTime(0.0);
                            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                            auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                            p_microvessel_solver->SetVesselNetwork(p_network);
                            p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                            p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                            p_microvessel_solver->Run();

                            // Print the average oxygenation
                            std::vector<double> solution = p_oxygen_solver->GetSolution();
                            double average_oxygen = 0.0;
                            for(unsigned jdx=0;jdx<solution.size();jdx++)
                            {
                                average_oxygen += solution[jdx];
                            }
                            average_oxygen /= double(solution.size());
                            std::cout << "Average oxygen: " << average_oxygen << std::endl;
                            // }
                            // Write the PQs
                            outfile.open("/tmp/narain/testoutput/TestDichotomousNetwork/forking_nonrandom_threshold_pruning_perfusion_quotients.txt", std::ios_base::app);
                            outfile << network_name << " " << solver_name << " " << lambda_string << " " << alpha_string << " " << mean_string << " " << threshold_string << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                            outfile.close();
                            
                            // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                            std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                            p_network->Write(output_file);

                            // Increase radius threshold
                            current_radius_threshold = current_radius_threshold + radius_step; 

                            // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                            ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                            ParameterCollection::Instance()->Destroy();
                            BaseUnits::Instance()->Destroy();
                            SimulationTime::Instance()->Destroy();
                            // }
                        }
                    }
                }
            }
        }
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a multi-generation forking network on a PDE grid with different h-splitting rules and stochastic pruning
    void xTestDichotomousNetworkWithStochasticPruningAndFlow2DAndVaryingMeans()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";

        // Set network in filename
        std::string network_name = "TestDichotomousNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("/tmp/narain/testoutput/TestDichotomousNetwork/forking_stochastic_pruning_perfusion_quotients.txt", std::ios_base::app);
        outfile.close();

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=4; h_solver++)
        {   
            // Generate the network for various lambdas
            for (unsigned k_aux=1; k_aux<2; k_aux++)  // generate the network for lambda = 4
            {
                // Set key vessel parameters
                double dimless_length = 1.0;  

                // ??? Non-dimensionalising the length
                for(unsigned i_aux=1; i_aux<order+1; i_aux++)
                {
                    dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
                }
   
                // Set up the reference length for the simulation
                QLength reference_length(1.0_um);
                BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

                // Set input radius
                QLength input_radius(7.5*GenericParameters::mpCapillaryRadius->GetValue());  // = 50_um
                // VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);

                // Set the viscosity
                QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
                // QDynamicViscosity viscosity = Owen11Parameters::mpPlasmaViscosity->GetValue();  // = 0.0012

                // Set the inlet and initial haematocrit
                double initial_haematocrit = 0.45;

                // Set threshold for perfusion quotient
		        QFlowRate threshold = 3.e-12*unit::metre_cubed_per_second;  // set flow threshold for what counts as a perfused vessel

                // Generate the networks
                VesselNetworkGenerator<2> network_generator;

                // Set lambda
                double lambda;  // lambda = length/diameter
                double twicelambda;  // used as an input parameter
                lambda = 2.0+double(k_aux)*2.0;
                twicelambda = 2.0*lambda;

                // Set the number of different means
                unsigned max_means = 3;

                // Run the simulation with different heterogeneities
                for (unsigned n_alpha=1; n_alpha<=max_alpha; n_alpha++)
                { 

                    // Offset the radii depending on the mean and SD
                    for (unsigned n_mean=0; n_mean<max_means; n_mean++)
                    { 

                        // Set up flag for broken solver
                        unsigned broken_solver = 0;

                        // Set alpha
                        double alpha = 1.0+(double)n_alpha*0.1;  // alpha determines the relative radius of the left vessel 

                        // Height of of first-order vessels
                        QLength main_vert_length = 0.9*twicelambda*input_radius*pow(2.0,-1.0/3.0);
                        
                        // Height of the domain
                        QLength domain_side_length_y = 4.0*main_vert_length+10.0_um;  // 10 um allows for enough room on top and bottom of network (considering the grid size is 10 um)

                        // Length of the domain
                        QLength domain_side_length_x = dimless_length*2.0*twicelambda*input_radius;

                        // Set the number of trials and the maximum beta value
                        int max_trials = 100;
                        int max_beta = 35;
                        
                        // Iterate over beta values
                        int beta = 1;
                        while (beta <= max_beta)
                        {
                            // Run the trials                    
                            int n_trial = 1;
                            while (n_trial <= max_trials)
                            {  
                                // Generate the network
                                std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateDichotomousNetworkForRadiotherapyRandomHeterogeneity(order, main_vert_length, input_radius, alpha, twicelambda);

                                // Identify input and output nodes and assign them properties
                                VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(0.0_um,2.0*main_vert_length));
                                VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network, Vertex<2>(domain_side_length_x, 2.0*main_vert_length));
                                p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
                                // p_network->GetNode(0)->GetFlowProperties()->SetPressure(Owen11Parameters::mpInletPressure->GetValue("User"));  // = 3333.06
                                p_inlet_node->GetFlowProperties()->SetPressure(3333.0_Pa);
                                p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
                                // p_network->GetNode(p_network->GetNumberOfNodes()-1)->GetFlowProperties()->SetPressure(Owen11Parameters::mpOutletPressure->GetValue("User"));  // = 1999.84
                                p_outlet_node->GetFlowProperties()->SetPressure(2000.0_Pa);

                                // Specify the minimum offsets for all alphas (75 um, 6 gens)
                                QLength min_array[] = {-1.77_um, -1.52_um, -0.87_um, 0.02_um, 1.04_um};
                                QLength mean_array[] = {3.97_um, 4.22_um, 4.87_um, 5.76_um, 6.78_um};
                                QLength max_array[] = {9.11_um, 9.36_um, 10.01_um, 10.9_um, 11.92_um};
                        
                                // Initialise the offsets
                                std::vector<QLength> min_list(min_array, min_array+5);
                                std::vector<QLength> mean_list(mean_array, mean_array+5);
                                std::vector<QLength> max_list(max_array, max_array+5);

                                // Set the radius offset based on the alpha
                                std::vector<std::shared_ptr<Vessel<2> > > vessels;
                                vessels = p_network->GetVessels();
                                if (n_mean==0)
                                    {
                                        QLength radius_offset = min_list[n_alpha];
                                        for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)
                                        {                            
                                            // Get the current segment's radius
                                            QLength current_radius = vessels[vessel_index]->GetRadius();

                                            // Calculate the new radius
                                            QLength new_vessel_radius = current_radius + (radius_offset*0.5);

                                            // Set the new radius
                                            vessels[vessel_index]->SetRadius(new_vessel_radius);
                                        }
                                    }
                                else if (n_mean==1)
                                {
                                    QLength radius_offset = mean_list[n_alpha];
                                    for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)
                                    {                            
                                        // Get the current segment's radius
                                        QLength current_radius = vessels[vessel_index]->GetRadius();

                                        // Calculate the new radius
                                        QLength new_vessel_radius = current_radius + (radius_offset*0.5);
        
                                        // Set the new radius
                                        vessels[vessel_index]->SetRadius(new_vessel_radius);
                                    }
                                }
                                else if (n_mean==2)
                                {
                                    QLength radius_offset = max_list[n_alpha];
                                    for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)
                                    {                            
                                        // Get the current segment's radius
                                        QLength current_radius = vessels[vessel_index]->GetRadius();

                                        // Calculate the new radius
                                        QLength new_vessel_radius = current_radius + (radius_offset*0.5);

                                        // Set the new radius
                                        vessels[vessel_index]->SetRadius(new_vessel_radius);
                                    }
                                }
                                vessels = p_network->GetVessels();

                                // Set the h-solver
                                std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                                std::string solver_name;
                                if (h_solver==1)
                                {
                                    solver_name = "ConstantHaematocrit/"; 
                                    std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                                    auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                                    p_haematocrit_solver->SetVesselNetwork(p_network);
                                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                                    p_abstract_haematocrit_solver = p_haematocrit_solver;     
                                }
                                else if (h_solver==2)
                                {
                                    solver_name = "PriesHaematocrit/"; 
                                    std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                                    auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                                    p_haematocrit_solver->SetVesselNetwork(p_network);
                                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                                    p_abstract_haematocrit_solver = p_haematocrit_solver;                
                                }
                                else if (h_solver==3)
                                {
                                    solver_name = "MemoryHaematocrit/"; 
                                    std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                                    auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                                    p_haematocrit_solver->SetVesselNetwork(p_network);
                                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                                    p_abstract_haematocrit_solver = p_haematocrit_solver;     
                                }
                                else if (h_solver==4)
                                {
                                    solver_name = "FungHaematocrit/"; 
                                    std::cout << "Now using FungHaematocritSolver..." << std::endl;
                                    auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                                    // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                                    p_haematocrit_solver->SetVesselNetwork(p_network);
                                    p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                                    p_abstract_haematocrit_solver = p_haematocrit_solver;     
                                }
                                                
                                // Store the lambda values for file name
                                std::stringstream lambda_stream;
                                lambda_stream << std::fixed << std::setprecision(0) << lambda;
                                std::string lambda_string = lambda_stream.str();
                                std::string lambda_directory = network_name + solver_name + "/Lambda" + lambda_string; 

                                // Store the n_mean values for file name
                                std::stringstream mean_stream;
                                mean_stream << std::fixed << std::setprecision(0) << n_mean;
                                std::string mean_string = mean_stream.str();
                                // std::string lambda_directory = network_name + solver_name + "/Lambda" + lambda_string; 

                                // Set up the grid for the finite difference solver
                                auto p_grid = RegularGrid<2>::Create();
                                // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
                                // QLength grid_spacing(1.0*1_um);  // the simulation time gets quite long if you reduce the resolution further
                                // QLength grid_spacing = 10_um;
                                p_grid->SetSpacing(grid_spacing);
                                c_vector<unsigned, 3> dimensions;
                                dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                                dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                                dimensions[2] = 1;
                                p_grid->SetDimensions(dimensions);

                                // Choose the PDE
                                auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                                
                                // Set the diffusivity and decay terms
                                p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                                p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                                // Set up the discrete source
                                auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                                QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                                        GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                                QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                                        Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                                p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                                p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                                p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                                p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                                // Set up the finite difference solver for oxygen (which handles everything)
                                auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                                p_oxygen_solver->SetPde(p_oxygen_pde);
                                p_oxygen_solver->SetLabel("oxygen");
                                p_oxygen_solver->SetGrid(p_grid);

                                // Set up the viscosity calculator
                                auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                                p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                                p_viscosity_calculator->SetVesselNetwork(p_network);
                                p_viscosity_calculator->Calculate();
                                
                                // Set up the impedance calculator
                                auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                                p_impedance_calculator->SetVesselNetwork(p_network);
                                p_impedance_calculator->Calculate();

                                // Set up the flow solver
                                FlowSolver<2> flow_solver;
                                flow_solver.SetVesselNetwork(p_network);
                                flow_solver.SetUp();
    
                                // Conduct stochastic pruning (thinner vessels are more likely to be pruned)
                                p_network->StochasticPruning(beta, false);  

                                // Set filename
                                std::stringstream alpha_stream;
                                std::stringstream trial_stream;
                                std::stringstream beta_stream;
                                alpha_stream << std::fixed << std::setprecision(2) << alpha;
                                trial_stream << std::fixed << std::setprecision(0) << n_trial;
                                beta_stream << std::fixed << std::setprecision(0) << beta;
                                std::string alpha_string = alpha_stream.str();
                                std::string trial_string = trial_stream.str();
                                std::string beta_string = beta_stream.str();
                                std::string file_name = lambda_directory + "/Alpha" + alpha_string + "/Mean" + mean_string + "/Beta" + beta_string + "/Trial" + trial_string;
                                std::string str_directory_name = file_name;
                                auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                                // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance. The loop runs the solvers, calculates the max. difference in the network between the previous H of a vessel and the new H in the vessel, and repeats this until all vessels converge on a value within a certain tolerance.
                                unsigned max_iter = 100;  // simulation fails if it doesn't converge in these many iterations
                                double tolerance2 = 0.01; //1.e-10
                                std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                                std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                                for (unsigned idx=0;idx<max_iter;idx++)
                                {
                                    // Run the solvers
                                    p_impedance_calculator->Calculate();
                                    flow_solver.SetUp();
                                    flow_solver.Solve();
                                    p_abstract_haematocrit_solver->Calculate();
                                    p_viscosity_calculator->Calculate();

                                    // Check for convergence 
                                    double max_difference = 0.0;
                                    double h_for_max = 0.0;
                                    double prev_for_max = 0.0;
                                    for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)  // for all the segments 
                                    {
                                        // Set segments with no flow to have no haematocrit
                                        if (fabs(segments[segment_index]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                        {
                                            segments[segment_index]->GetFlowProperties()->SetHaematocrit(0.0);
                                        }                                    
                                        double current_haematocrit = segments[segment_index]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                                        double difference = std::abs(current_haematocrit - previous_haematocrit[segment_index]);  // difference in new H in segment and old H
                                        if(difference>max_difference)  // if the difference is greater than previous max.
                                        {
                                            max_difference = difference;
                                            h_for_max = current_haematocrit;
                                            prev_for_max = previous_haematocrit[segment_index];
                                        }
                                        previous_haematocrit[segment_index] = current_haematocrit;
                                    }
                                    std::cout << "Segment H. at max difference: " << h_for_max << ", Prev. segment H at max difference:" << prev_for_max << std::endl;
                                    if(max_difference<=tolerance2)  
                                    {
                                        std::cout << "Converged after: " << idx << " iterations." <<  std::endl;
                                        broken_solver = 0;
                                        break;
                                    }
                                    else
                                    {
                                        if(idx%1==0)
                                        {
                                            std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                                            // std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                                            // std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                                            // p_network->Write(output_file);
                                        }
                                    }
                                    if(idx==max_iter-1)  // if there is no convergence after all the iterations, print the error message
                                    {
                                        std::cout << "Problem encountered in Dichotomous Network with h_solver = " << h_solver << " using n_alpha = " << n_alpha << " and lambda = " << lambda << " in trial = " << n_trial << std::endl;
                                        // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                                        error_log << "\n Problem encountered in Dichotomous Network with h_solver = " << h_solver << " using n_alpha = " << n_alpha << " and lambda = " << lambda << " in trial = " << n_trial; 
                                        broken_solver = 1;
                                        break;
                                    }
                                }

                                // If solver doesn't converge, move on to next one
                                if (broken_solver == 1)
                                {
                                    continue;
                                }

                                // Run the simulation 
                                SimulationTime::Instance()->SetStartTime(0.0);
                                SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                                auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                                p_microvessel_solver->SetVesselNetwork(p_network);
                                p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                                p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                                p_microvessel_solver->Run();

                                // Print the average oxygenation
                                std::vector<double> solution = p_oxygen_solver->GetSolution();
                                double average_oxygen = 0.0;
                                for(unsigned jdx=0;jdx<solution.size();jdx++)
                                {
                                    average_oxygen += solution[jdx];
                                }
                                average_oxygen /= double(solution.size());
                                std::cout << "Average oxygen: " << average_oxygen << std::endl;

                                // Write the PQs
                                outfile.open("/tmp/narain/testoutput/TestDichotomousNetwork/forking_random_stochastic_pruning_perfusion_quotients.txt", std::ios_base::app);
                                outfile << network_name << " " << solver_name << " " << lambda_string << " " << alpha_string << " "  << mean_string  << " " << beta_string << " " << trial_string << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                                outfile.close();
                                
                                // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                                std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                                p_network->Write(output_file);

                                // // ??? Assuming radius doesn't exceed threshold, write the network to the output file
                                // if(current_radius_threshold>12.9_um && current_radius_threshold<13.1_um)
                                // {
                                //     p_network->Write(output_file_final_RT);
                                // }
                                n_trial = n_trial + 1;  // move onto next trial

                                // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                                ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                                ParameterCollection::Instance()->Destroy();
                                BaseUnits::Instance()->Destroy();
                                SimulationTime::Instance()->Destroy();
                            }
                            // Move on to next beta
                            beta = beta + 1;
                        }
                    }
                }
            }
        }
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Voronoi Network
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Choose key parameters
    unsigned NumberOfSeedPoints = 400;  // change this to select which Voronoi architecture to use: 25, 100, 400
    unsigned NumberOfLayouts = 2;  // number of different point layouts to run simulations with (max. = 100)
    // double dimless_domain_size_x = 2050.0;  // x-coordinate of output node

	// unsigned thin_selections = 100;  // number of layouts from which to choose thin vessels (upto 100)
    // unsigned NHet = 5;  // number of heterogeneity selections (upto 5)
    // QLength large_vessel_radius = 10_um;
    // QLength small_vessel_radius = 5_um;    
    // QLength inlet_vessel_radius = 75_um;  // the maximum radius of the distribution
    // double dimless_domain_size_x = 2000.0; 
    // double dimless_domain_size_y = 2000.0 - 86.6025; 
    // double dimless_domain_size_y = 1818.65 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
    // unsigned dimless_vessel_length = 100.0;

    // Make a Voronoi network on a PDE grid as a Dirichlet BC in 2D 
    void xTestVoronoiNetworkLineSource2D()
    {
        // Run simulations with different layouts (could be used to compute some average later)
        for (unsigned layout=1;layout < NumberOfLayouts+1; layout++)   
        {

            // Set up the reference length for the simulation
            QLength reference_length(1.0_um);
            BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

            // Set key parameters
            // unsigned NumberOfSeedPoints = 400;  // change this to select which Voronoi architecture to use: 25, 100, 400
            // QLength input_radius(10.0 *GenericParameters::mpCapillaryRadius->GetValue());  // = 50_um
            // VesselNetworkPropertyManager<2>::SetSegmentRadii(p_network, vessel_radius);

            // Initialise the simulation space
            VesselNetworkGenerator<2> network_generator;
            double dimless_domain_size_x = 2000.0; 
            QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
            
            // Initialise the .vtp outputs
            auto p_file_handler = std::make_shared<OutputFileHandler>("TestVoronoiNetwork/DirichletHaematocrit/NoPruning/SeedPoints"+to_string(NumberOfSeedPoints)+"/Layout"+to_string(layout), true);  // create a folder
            // std::string output_file_initial = p_file_handler->GetOutputDirectoryFullPath().append("InitialNetwork.vtp");
            std::string output_file_final = p_file_handler->GetOutputDirectoryFullPath().append("FinalNetwork.vtp");
            // std::string output_file_final_RT = p_file_handler->GetOutputDirectoryFullPath().append("PostRTNetwork.vtp");

            // Set the simulation paramaters
            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
            // QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;  // minimum flow
            std::ofstream outfileMean;
            // unsigned NumberOfCutOffs = 10;  // upper limit for varying pruning length threshold
            // std::vector<std::vector<double> > QuotientMatrixMean(1,std::vector<double>(NumberOfCutOffs+1));
            // std::vector<std::vector<double> > QuotientMatrixAggregate(1,std::vector<double>(NumberOfCutOffs+1,0.0));

            // Read the network layout from a file
            std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(NumberOfSeedPoints)+"SeedPoints/EdgesMatrixSampleNumber"+to_string(layout)+".txt");
            std::vector<std::vector<double> > rEdgesMatrix;
            string line;

            // Break down the rows in the file into column values
            while (std::getline(in, line)) 
            {
                rEdgesMatrix.push_back(std::vector<double>());
                std::stringstream split(line);
                double value;
                while (split >> value)
                {
                        rEdgesMatrix.back().push_back(value);
                }
            }

            // Generate the network
            std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateVoronoiNetwork(rEdgesMatrix);

            // Set up the grid for the finite difference solver
            auto p_grid = RegularGrid<2>::Create();
            // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
            // QLength grid_spacing(1.0*1_um);  // the simulation time gets quite long if you reduce the resolution further
            // QLength grid_spacing = 10_um;
            p_grid->SetSpacing(grid_spacing);
            c_vector<unsigned, 3> dimensions;
            dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
            dimensions[1] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num_y
            dimensions[2] = 1;
            p_grid->SetDimensions(dimensions);
            // std::shared_ptr<Part<2> > p_domain = Part<2>::Create();
            // p_domain->AddRectangle(domain_x, domain_y);
            // std::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
            // p_grid->GenerateFromPart(p_domain, grid_spacing);  // set domain and grid spacing

            // Choose the PDE
            auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
            
            // Set the diffusivity and decay terms
            p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
            p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

            // Set the O2 concentration value
            QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") * GenericParameters::mpGasConcentrationAtStp->GetValue("User");
            QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp * Owen11Parameters::mpReferencePartialPressure->GetValue("User");

            // Set the boundary condition to be the network acting as a Dirichlet source
            std::shared_ptr<DiscreteContinuumBoundaryCondition<2> > p_vessel_boundary_condition = DiscreteContinuumBoundaryCondition<2>::Create();
            p_vessel_boundary_condition->SetValue(vessel_oxygen_concentration);  // set the O2 concentration at the source
            p_vessel_boundary_condition->SetType(BoundaryConditionType::VESSEL_LINE);
            p_vessel_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);

            // Set up the finite difference solver for oxygen (which handles everything)
            SimpleLinearEllipticFiniteDifferenceSolver<2> solver;
            solver.SetGrid(p_grid);
            solver.SetPde(p_oxygen_pde);
            solver.AddBoundaryCondition(p_vessel_boundary_condition);
            solver.SetVesselNetwork(p_network);
            solver.SetLabel("oxygen");
            solver.SetFileName("oxygen_solution_0");

            // Run the solver and write the output file (visualise with Paraview; set Filters->Alphabetical->Tube)
            solver.SetFileHandler(p_file_handler);
            solver.SetWriteSolution(true);
            solver.Solve();
            p_network->Write(p_file_handler->GetOutputDirectoryFullPath() + "voronoi_network_results.vtp");
            
            // Dump our parameter collection to an xml file and, importantly, clear it for the next test
            ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
            ParameterCollection::Instance()->Destroy();
            BaseUnits::Instance()->Destroy();
        }
    }

    // Make a 2D Voronoi network on a PDE grid with flow and H-splitting
    void xTestVoronoiNetworkWithFlow2D() 
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";

        // Set network in filename
        std::string network_name = "TestVoronoiNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("voronoi_lognormal_individual_perfusion_quotients.txt");
        outfile.close();

        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Define the key parameters
        // double dimless_domain_size_x = 2050.0;  // x-coordinate of output node
        // double dimless_domain_size_y = 1818.65 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        unsigned max_n_mu = 3;
        // unsigned n_vessels = 382;  // number of non-inlet/outlet vessels from which to select ones to make thin
        // double percToKill = 0.2;  // percentage of vessels to kill
        // unsigned ToBeKilled = (unsigned)(percToKill*n_vessels);  // number to kill
        unsigned ToBeKilled = 100;  // number to kill

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=1; h_solver++)
        {     
            // Initialise the simulation
            std::shared_ptr<VesselNetwork<2> > p_network;
            std::vector<std::shared_ptr<Vessel<2> > > vessels;
            VesselPtr<2> p_current_vessel;
            QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
            QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
            std::vector<std::vector<double> > QuotientMatrixMean(max_n_mu,std::vector<double>(ToBeKilled+1));
            std::vector<std::vector<double> > QuotientMatrixAggregate(max_n_mu,std::vector<double>(ToBeKilled+1,0.0));
            std::vector<std::vector<unsigned> > Order;
            double PerfusionQuotient;

            // Run simulations with different layouts (could be used to compute some average later)
            for (unsigned layout=1;layout < NumberOfLayouts+1; layout++)   
            { 
                // Read the network layout from a file
                VesselNetworkGenerator<2> network_generator;
                std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(NumberOfSeedPoints)+"SeedPoints/EdgesMatrixSampleNumber"+to_string(layout)+".txt");
                std::vector<std::vector<double> > rEdgesMatrix;
                string line;
                while (std::getline(in, line)) 
                {
                    rEdgesMatrix.push_back(std::vector<double>());
                    std::stringstream split(line);
                    double value;
                    while (split >> value)
                    {
                            rEdgesMatrix.back().push_back(value);
                    }
                }

                // Loop over increasing level of heteregeneity 
                for(unsigned n_mu=0; n_mu<=max_n_mu; n_mu++)
                {
                    // Convert to actual mu
                    int mu = (3*n_mu)+10;  // alpha determines the mean
                    
                    // Choose a radii list file and read the file to an array
                    std::vector<std::vector<double> > radii_array;
                    radii_array.clear();
                    // int list_index=0;
                    string line_1;
                    std::ifstream radius_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/100SeedPoints/square_voronoi_radius_log_normal_distribution/mu_"+to_string(mu)+"/radii_list_"+to_string(layout)+".txt");
                    while (std::getline(radius_list_file, line_1)) 
                    {
                        radii_array.push_back(std::vector<double>());
                        std::stringstream split(line_1);
                        double value;
                        while (split >> value)
                        {
                            radii_array.back().push_back(value);
                        }
                    }

                    // Generate the network
                    std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateVoronoiNetwork(rEdgesMatrix);

                    // Assign the radii
                    vessels = p_network->GetVessels();
                    for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)  // for all the segments in the network
                    {
                        // Get segment radius from array
                        QLength new_vessel_radius(double(radii_array[vessel_index][0])*unit::microns);
                        vessels[vessel_index]->SetRadius(new_vessel_radius);
                    }

                    // // Set inlet and outlet nodes (if node (start/end) is at the edge of the domain space, make it input or output depending on which side), and assign each non-inlet/outlet vessel a radius from the list  
                    std::vector<std::shared_ptr<Vessel<2> > > vessels = p_network->GetVessels();
                    for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
                    {
                        if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
                        {
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] <  tolerance)
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
                            }
                            //if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
                            else
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
                            }
                            (*vessel_iterator)->SetRadius(inlet_vessel_radius);
                        }
                        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)
                        {
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] <  tolerance)
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
                            }
                            //if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
                            else
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
                            }
                            (*vessel_iterator)->SetRadius(inlet_vessel_radius);
                        }
                    }
                    vessels = p_network->GetVessels();
            
                    // Set the haematocrit solver
                    std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                    std::string solver_name;
                    if (h_solver==1)
                    {
                        solver_name = "ConstantHaematocrit"; 
                        std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==2)
                    {
                        solver_name = "PriesHaematocrit"; 
                        std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                
                    }
                    else if (h_solver==3)
                    {
                        solver_name = "MemoryHaematocrit"; 
                        std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==4)
                    {
                        solver_name = "FungHaematocrit"; 
                        std::cout << "Now using FungHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==5)
                    {
                        solver_name = "YangHaematocrit"; 
                        std::cout << "Now using YangHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();   
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }

                    // Set up the viscosity solver
                    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                    p_viscosity_calculator->SetVesselNetwork(p_network);
                    p_viscosity_calculator->Calculate();
                    
                    // Set up the impedance solver
                    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                    p_impedance_calculator->SetVesselNetwork(p_network);
                    p_impedance_calculator->Calculate();

                    // Set up the flow solver 
                    FlowSolver<2> flow_solver;
                    flow_solver.SetVesselNetwork(p_network);
                    // flow_solver.SetUp();
                    flow_solver.SetUseDirectSolver(true);
                    // flow_solver.Solve();

                    // Prune all vessels up to specified dose 
                    for(unsigned KilledVessels=0; KilledVessels < ToBeKilled+1; KilledVessels++)
                    { 
                        // Display status message
                        std::cout << "Now killed " << KilledVessels << " vessels." << std::endl;  

                        // Set filename
                        std::stringstream selection_stream;
                        std::stringstream sd_stream;
                        std::stringstream mean_stream;
                        std::stringstream kill_stream;
                        // selection_stream << std::fixed << std::setprecision(0) << to_string(list_number);                    
                        selection_stream << std::fixed << std::setprecision(0) << to_string(layout);                    
                        // sd_stream << std::fixed << std::setprecision(0) << sigma;                     
                        mean_stream << std::fixed << std::setprecision(0) << mu;                                       
                        kill_stream << std::fixed << std::setprecision(0) << KilledVessels;           
                        std::string selection_string = selection_stream.str();
                        // std::string sd_string = sd_stream.str();                    
                        std::string mean_string = mean_stream.str();                    
                        std::string kill_string = kill_stream.str();                    
                        // std::string file_name = network_name + solver_name + "/Sigma" + sd_string + "/Mu" + mean_string + "/Selection" + selection_string + "/Kills" + kill_string;
                        std::string file_name = network_name + solver_name + "/Mu" + mean_string + "/Selection" + selection_string + "/Kills" + kill_string;
                        std::string str_directory_name = file_name;
                        std::cout << str_directory_name << std::endl;
                        auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                        // Set up flag for broken solver
                        unsigned broken_solver = 0;

                        // // Set the simulation paramaters
                        // double tolerance = 0.001;  // tolerance for solution convergence ???
                        // // QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;  // minimum flow
                        // std::ofstream outfileMean;
                        // // unsigned NumberOfCutOffs = 10;  // upper limit for varying pruning length threshold
                        // // std::vector<std::vector<double> > QuotientMatrixMean(1,std::vector<double>(NumberOfCutOffs+1));
                        // // std::vector<std::vector<double> > QuotientMatrixAggregate(1,std::vector<double>(NumberOfCutOffs+1,0.0));

                        // // Set up the grid for the finite difference solver
                        // auto p_grid = RegularGrid<2>::Create();
                        // // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
                        // // QLength grid_spacing(1.0*1_um);  // the simulation time gets quite long if you reduce the resolution further
                        // // QLength grid_spacing = 10_um;
                        // p_grid->SetSpacing(grid_spacing);
                        // c_vector<unsigned, 3> dimensions;
                        // dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                        // dimensions[1] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num_y
                        // dimensions[2] = 1;
                        // p_grid->SetDimensions(dimensions);
                        // // std::shared_ptr<Part<2> > p_domain = Part<2>::Create();
                        // // p_domain->AddRectangle(domain_x, domain_y);
                        // // std::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
                        // // p_grid->GenerateFromPart(p_domain, grid_spacing);  // set domain and grid spacing

                        // // Choose the PDE
                        // auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                        
                        // // Set the diffusivity and decay terms
                        // p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                        // p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                        // // Set up the discrete source
                        // auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                        // QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                        //         GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                        // QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                        //         Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                        // p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                        // p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                        // p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        // p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                        // // Set up the finite difference solver for oxygen (which handles everything)
                        // auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                        // p_oxygen_solver->SetPde(p_oxygen_pde);
                        // p_oxygen_solver->SetLabel("oxygen");
                        // p_oxygen_solver->SetGrid(p_grid);
                        // // solver.SetFileName("oxygen_solution_0");

                        // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                        // unsigned max_iter = 1000;  // 1000 
                        // double tolerance2 = 1.e-10;
                        std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                        // std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                        // for(unsigned idx=0;idx<max_iter;idx++)
                        // {
                        // Run the solvers
                        p_impedance_calculator->Calculate();
                        p_abstract_haematocrit_solver->Calculate();
                        p_viscosity_calculator->Calculate();
                        flow_solver.SetUp();
                        flow_solver.Solve();

                        // Get the residual
                        // double max_difference = 0.0;
                        // double h_for_max = 0.0;
                        // double prev_for_max = 0.0;
                        for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                        {
                            // Set segments with no flow to have no haematocrit
                            if (fabs(segments[jdx]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                            {
                                segments[jdx]->GetFlowProperties()->SetHaematocrit(0.0);
                            }   
                            
                            // double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                            // double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                            // if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                            // {
                            //     max_difference = difference;
                            //     h_for_max = current_haematocrit;
                            //     prev_for_max = previous_haematocrit[jdx];
                            // }
                            // previous_haematocrit[jdx] = current_haematocrit;
                        }
                        //     std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                        //     // Print the final or intermediary convergence results
                        //     if(max_difference<=tolerance2)  
                        //     {
                        //         std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                        //         broken_solver = 0;
                        //         break;
                        //     }
                        //     else
                        //     {
                        //         if(idx%1==0)
                        //         {
                        //             std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                        //             // std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                        //             // std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                        //             // p_network->Write(output_file);
                        //         }
                        //     }

                        //     // If there is no convergence after all the iterations, print the error message.
                        //     if(idx==max_iter-1)
                        //     {
                        //         std::cout << "Problem encountered in Voronoi network with h_solver = " << h_solver << " using layout = " << layout << std::endl;
                        //         error_log << "\n Problem encountered in Voronoi network with h_solver = " << to_string(h_solver) << " using layout = " << to_string(layout); 
                        //         // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                        //         broken_solver = 1;
                        //         break;
                        //     }
                        // }

                        // If solver doesn't converge, move on to next one
                        if (broken_solver == 1)
                        {
                            continue;
                        }

                        // Run the simulation 
                        // SimulationTime::Instance()->SetStartTime(0.0);
                        // SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                        // auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                        // p_microvessel_solver->SetVesselNetwork(p_network);
                        // p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                        // p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                        // p_microvessel_solver->Run();
                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                        p_network->Write(output_file);

                        // Print the average oxygenation
                        // std::vector<double> solution = p_oxygen_solver->GetSolution();
                        // double average_oxygen = 0.0;
                        // for(unsigned jdx=0;jdx<solution.size();jdx++)
                        // {
                        //     average_oxygen += solution[jdx];
                        // }
                        // average_oxygen /= double(solution.size());
                        // std::cout << "Average oxygen: " << average_oxygen << std::endl;
                        
                        // // Write the PQs
                        // QFlowRate threshold = 5.e-14*unit::metre_cubed_per_second;
                        QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;
                        // QFlowRate threshold = 2.e-13*unit::metre_cubed_per_second; 
                        PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
                        // QuotientMatrixAggregate[i][KilledVessels] += PerfusionQuotient;  // store the PQ with the associated NHet and number of kills
                        outfile.open("voronoi_lognormal_individual_perfusion_quotients.txt", std::ios_base::app);
                        // outfile << network_name << " " << solver_name << " " << heterogeneity_string << " " << selection_string << " " << kill_string << " " << PerfusionQuotient << " \n"; 
                        outfile << network_name << " " << solver_name << " " << mean_string << " " << selection_string << " " << kill_string << " " << PerfusionQuotient << " \n"; 
                        outfile.close();

                        // Remove the smallest vessel
                        vessels = p_network->GetVessels();
                        QLength minimum_radius = inlet_vessel_radius;
                        unsigned int minimum_index = 0;
                        for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)  // for all the segments in the network
                        {
                            // Get the current segment's radius
                            QLength current_radius = vessels[vessel_index]->GetRadius();

                            // If the current radius is less than the minimum radius, record the new minimum
                            if (current_radius < minimum_radius)
                            {
                                minimum_radius = current_radius;
                                minimum_index = vessel_index;
                            }                  
                        }
                        p_network->RemoveVessel(vessels[minimum_index], true);  // remove the vessel

                        // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                        ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                        ParameterCollection::Instance()->Destroy();
                        BaseUnits::Instance()->Destroy();
                        SimulationTime::Instance()->Destroy();
                    }
            
                } 
            }
        }

        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hexagonal Network 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Choose heterogeneity parameters
	unsigned thin_selections = 100;  // number of layouts from which to choose thin vessels (upto 100)
    unsigned NHet = 5;  // number of heterogeneity selections (upto 5)
    QLength large_vessel_radius = 10_um;
    QLength small_vessel_radius = 5_um;    
    QLength inlet_vessel_radius = 37.5_um;  // the maximum radius of the distribution
    // double dimless_domain_size_x = 2000.0; 
    // double dimless_domain_size_y = 2000.0 - 86.6025; 
    double dimless_domain_size_x = 2050.0;  // x-coordinate of output node
    double dimless_domain_size_y = 1818.65 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
    unsigned dimless_vessel_length = 100.0;

    // Make a hexagonal network on a PDE grid as a Dirichlet BC in 2D 
    void xTestHexagonalNetworkLineSource2D()
    {
        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);
       
        // Define the key parameters
        std::vector<std::vector<unsigned> > Order;
        std::ofstream outfileMean;
        unsigned n_vessels = 386;  // number of non-inlet/outlet vessels from which to select ones to make thin
        std::shared_ptr<VesselNetwork<2> > p_network;
        std::vector<std::shared_ptr<Vessel<2> > > vessels;
        string line2;
        VesselPtr<2> p_current_vessel;
        unsigned ToBeThin;
        unsigned NThin;
        std::vector<std::vector<double> > rEdgesMatrix;
        string line;

        // Initialise the simulation space
        QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
        QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
        std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;

        // Read the network from a file
        VesselNetworkGenerator<2> network_generator;
        std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix.txt");

        // // Create a folder for the outputs
        // auto p_file_handler = std::make_shared<OutputFileHandler>("TestHexagonalNetwork/DirichletHaematocrit/NoPruning/", true);

        // Break down the rows in the file into column values
        while (std::getline(in, line)) 
        {
            rEdgesMatrix.push_back(std::vector<double>());
            std::stringstream split(line);
            double value;
            while (split >> value)
            {
                rEdgesMatrix.back().push_back(value);
            }
        }

        // Loop over different thin vessel layouts
        for(unsigned att=0; att < thin_selections;att++)
	    {
            // Choose a layout for the thin vessels
            outfileMean.open("HexagonalFromVoronoiMean.txt");
	        std::ifstream in2("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/Selection"+to_string(att+1)+".txt");
        	Order.clear();

            // Break down the row into column values
            while (std::getline(in2, line2)) 
	        {
				Order.push_back(std::vector<unsigned>());
			    std::stringstream split2(line2);
			    unsigned value2;
				while (split2 >> value2)
		        {
					Order.back().push_back(value2);
			    }
		    }

		    // Loop over increasing level of heteregeneity (% of thin vessels)
            for(unsigned i=0; i<NHet;i++)
		    {
                // Initialise outputs
                // std::string output_file_two_radii = p_file_handler->GetOutputDirectoryFullPath().append("TwoRadiiNetwork/NHet"+to_string(NHet)+"/Selection"+to_string(i)+"/hexagonal_network_results.vtp");

                std::ostringstream strs;
                strs << std::fixed << std::setprecision( 2 );
                strs << "TestHexagonalNetwork/NoPruning/DirichletHaematocrit/Selection" << att << "/NHet" << i;
                std::string str_directory_name = strs.str();
                auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);         

                // Generate the network
                p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);
                auto p_segment = p_network->GetVesselSegments()[0];
                p_segment->SetRadius(large_vessel_radius);
                // p_segment->GetFlowProperties()->SetHaematocrit(inlet_haematocrit);
                VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   
                double percOfThin = 0.05*i;
                ToBeThin = (unsigned)(percOfThin*n_vessels);
                NThin = 0;

                // Make a certain % of vessels thin
                while( NThin < ToBeThin)
                {
                    p_current_vessel = p_network->GetVessel(Order[NThin][0]-1);
                    p_current_vessel->SetRadius(small_vessel_radius);
                    NThin++;
                }

                // Set up the grid for the finite difference solver
                auto p_grid = RegularGrid<2>::Create();
                // QLength grid_spacing(Owen11Parameters::mpLatticeSpacing->GetValue("User"));
                // QLength grid_spacing(1.0*1_um);  // the simulation time gets quite long if you reduce the resolution further
                // QLength grid_spacing = 10_um;
                p_grid->SetSpacing(grid_spacing);
                c_vector<unsigned, 3> dimensions;
                dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                dimensions[2] = 1;
                p_grid->SetDimensions(dimensions);

                // Choose the PDE
                auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();

                // Set the diffusivity and decay terms
                p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                // Set the O2 concentration value
                QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") * GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp * Owen11Parameters::mpReferencePartialPressure->GetValue("User");

                // Set the boundary condition to be the network acting as a Dirichlet source
                std::shared_ptr<DiscreteContinuumBoundaryCondition<2> > p_vessel_boundary_condition = DiscreteContinuumBoundaryCondition<2>::Create();
                p_vessel_boundary_condition->SetValue(vessel_oxygen_concentration);  // set the O2 concentration at the source
                p_vessel_boundary_condition->SetType(BoundaryConditionType::VESSEL_LINE);
                p_vessel_boundary_condition->SetSource(BoundaryConditionSource::PRESCRIBED);

                // Set up the finite difference solver for oxygen (which handles everything)
                SimpleLinearEllipticFiniteDifferenceSolver<2> solver;
                solver.SetGrid(p_grid);
                solver.SetPde(p_oxygen_pde);
                solver.AddBoundaryCondition(p_vessel_boundary_condition);
                solver.SetVesselNetwork(p_network);
                solver.SetLabel("oxygen");
                solver.SetFileName("oxygen_solution_0");
                
                // Run the solver and write the output file (visualise with Paraview; set Filters->Alphabetical->Tube)
                solver.SetFileHandler(p_file_handler);
                solver.SetWriteSolution(true);
                solver.Solve();
                p_network->Write(p_file_handler->GetOutputDirectoryFullPath() + "hexagonal_network_results.vtp");
                
                // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                ParameterCollection::Instance()->Destroy();
                BaseUnits::Instance()->Destroy();
            }
        }
    }

    // Make a 2D hexagonal network on a PDE grid with flow and H-splitting
    void xTestHexagonalNetworkWithFlow2D()
    {
        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Define the key parameters
        std::vector<std::vector<unsigned> > Order;
        unsigned n_vessels = 386;  // number of non-inlet/outlet vessels from which to select ones to make thin
        std::shared_ptr<VesselNetwork<2> > p_network;
        std::vector<std::shared_ptr<Vessel<2> > > vessels;
        string line2;
        VesselPtr<2> p_current_vessel;
        unsigned ToBeThin;
        unsigned NThin;

        // Initialise the simulation space
        QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
        QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
        std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
        std::ofstream outfileMean;

        // Define the key parameters
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        double initial_haematocrit = 0.45;
        double tolerance = 0.001;
        
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";

        // Read the network from a file
        VesselNetworkGenerator<2> network_generator;
        std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix.txt");
        std::vector<std::vector<double> > rEdgesMatrix;
        string line;
        
        // Break down the rows in the file into column values
        while (std::getline(in, line)) 
        {
            rEdgesMatrix.push_back(std::vector<double>());
            std::stringstream split(line);
            double value;
            while (split >> value)
            {
                rEdgesMatrix.back().push_back(value);
            }
        }

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<6; h_solver++)
        {
            // Loop over different thin vessel layouts
            for(unsigned att=0; att < thin_selections;att++)
            {
                // Choose a layout for the thin vessels
                outfileMean.open("HexagonalFromVoronoiMean.txt");
                std::ifstream in2("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/Selection"+to_string(att+1)+".txt");
                Order.clear();

                // Break down the row into column values
                while (std::getline(in2, line2)) 
                {
                    Order.push_back(std::vector<unsigned>());
                    std::stringstream split2(line2);
                    unsigned value2;
                    while (split2 >> value2)
                    {
                        Order.back().push_back(value2);
                    }
                }

                // Loop over increasing level of heteregeneity (% of thin vessels)
                for(unsigned i=0; i<NHet;i++)
                {

                    // Set file name
                    std::ostringstream strs;
                    strs << std::fixed << std::setprecision( 2 );
                    strs << "TestHexagonalNetwork/NoPruning/";        

                    // std::ostringstream strs;
                    if (h_solver==1)
                    {
                        strs << "ConstantHaematocrit";
                    }
                    else if (h_solver==2)
                    {
                        strs << "PriesHaematocrit";
                    }
                    else if (h_solver==3)
                    {
                        strs << "MemoryHaematocrit";
                    }
                    else if (h_solver==4)
                    {
                        strs << "FungHaematocrit";
                    }
                    else if (h_solver==5)
                    {
                        strs << "YangHaematocrit";
                    }
                    strs << "/Selection" << att << "/NHet" << i;
                    std::string str_directory_name = strs.str();
                    auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                    // Generate the network
                    p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);
                    
                    // Set inlet and outlet nodes
                    auto p_segment = p_network->GetVesselSegments()[0];
                    vessels = p_network->GetVessels();
                    for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
                    {
                        if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
                        {
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < -0.5*dimless_vessel_length+tolerance)
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
                            }
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
                            }
                        }
                        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)   	
                        //std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
                        {
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] < -0.5*dimless_vessel_length+tolerance)
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
                            }
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
                            }
                        }
                    }
                    p_segment = p_network->GetVesselSegments()[0];
                    p_segment->SetRadius(large_vessel_radius);
                    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   

                    double percOfThin = 0.05*i;
                    ToBeThin = (unsigned)(percOfThin*n_vessels);
                    NThin = 0;
                    
                    // Make a certain % of vessels thin
                    while( NThin < ToBeThin)
                    {
                        p_current_vessel = p_network->GetVessel(Order[NThin][0]-1);
                        p_current_vessel->SetRadius(small_vessel_radius);
                        NThin++;
                    }

                    // Set the haematocrit solver
                    std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                    if (h_solver==1)
                    {
                        std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;       
                    }
                    else if (h_solver==2)
                    {
                        std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                    
                    }
                    else if (h_solver==3)
                    {
                        std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                        // auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;      
                    }
                    else if (h_solver==4)
                    {
                        std::cout << "Now using FungHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;      
                    }
                    else if (h_solver==5)
                    {
                        std::cout << "Now using YangHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;      
                    }

                    // Set up the viscosity solver
                    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                    p_viscosity_calculator->SetVesselNetwork(p_network);
                    p_viscosity_calculator->Calculate();
                    
                    // Set up the impedance solver
                    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                    p_impedance_calculator->SetVesselNetwork(p_network);
                    p_impedance_calculator->Calculate();

                    // Set up the flow solver 
                    FlowSolver<2> flow_solver;
                    flow_solver.SetVesselNetwork(p_network);
                    // flow_solver.SetUp();
                    flow_solver.SetUseDirectSolver(true);
                    // flow_solver.Solve();

                    // Set up the grid for the finite difference solver
                    auto p_grid = RegularGrid<2>::Create();
                    p_grid->SetSpacing(grid_spacing);
                    c_vector<unsigned, 3> dimensions;
                    dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                    dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                    dimensions[2] = 1;
                    p_grid->SetDimensions(dimensions);

                    // Choose the PDE
                    auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                    
                    // Set the diffusivity and decay terms
                    p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                    p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                    // Set up the discrete source
                    auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                    QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                            GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                    QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                            Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                    p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                    p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                    p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                    // Set up the finite difference solver for oxygen (which handles everything)
                    auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                    p_oxygen_solver->SetPde(p_oxygen_pde);
                    p_oxygen_solver->SetLabel("oxygen");
                    p_oxygen_solver->SetGrid(p_grid);
                    // solver.SetFileName("oxygen_solution_0");

                    // Set up flag for broken solver
                    unsigned broken_solver = 0;

                    // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                    unsigned max_iter = 1000;  // 1000 
                    double tolerance2 = 1.e-10;
                    std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                    std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                    for(unsigned idx=0;idx<max_iter;idx++)
                    {
                        // Run the solvers
                        p_impedance_calculator->Calculate();
                        flow_solver.SetUp();
                        flow_solver.Solve();
                        p_abstract_haematocrit_solver->Calculate();
                        p_viscosity_calculator->Calculate();

                        // Get the residual
                        double max_difference = 0.0;
                        double h_for_max = 0.0;
                        double prev_for_max = 0.0;
                        for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                        {
                            double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                            double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                            if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                            {
                                max_difference = difference;
                                h_for_max = current_haematocrit;
                                prev_for_max = previous_haematocrit[jdx];
                            }
                            previous_haematocrit[jdx] = current_haematocrit;
                        }
                        std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                        // Print the final or intermediary convergence results
                        if(max_difference<=tolerance2)  
                        {
                            std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                            broken_solver = 0;
                            break;
                        }
                        else
                        {
                            if(idx%1==0)
                            {
                                std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                                std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                                std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                                p_network->Write(output_file);
                            }
                        }

                        // If there is no convergence after all the iterations, print the error message.
                        if(idx==max_iter-1)
                        {
                            std::cout << "Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(att) << " and NHet = " << i << std::endl;
                            error_log << "\n Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(att) << " and NHet = " << i; 
                            broken_solver = 1;
                            break;
                            // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                        }
                    }

                    // If solver doesn't converge, move on to next one
                    if (broken_solver == 1)
                    {
                        continue;
                    }

                    // Run the simulation 
                    SimulationTime::Instance()->SetStartTime(0.0);
                    SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                    auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                    p_microvessel_solver->SetVesselNetwork(p_network);
                    p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                    p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                    p_microvessel_solver->Run();

                    // Print the average oxygenation
                    std::vector<double> solution = p_oxygen_solver->GetSolution();
                    double average_oxygen = 0.0;
                    for(unsigned jdx=0;jdx<solution.size();jdx++)
                    {
                        average_oxygen += solution[jdx];
                    }
                    average_oxygen /= double(solution.size());
                    std::cout << "Average oxygen: " << average_oxygen << std::endl;
                    
                    // Write the output file (visualise with Paraview: set Filters->Alphabetical->Tube)
                    std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                    p_network->Write(output_file);

                    // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                    ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                    ParameterCollection::Instance()->Destroy();
                    BaseUnits::Instance()->Destroy();
                    SimulationTime::Instance()->Destroy();
                }
            }
        }
         
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a full 2D hexagonal network on a PDE grid with flow and H-splitting (correcting for the unexplained -50 um offset in the other simulation)
    void xTestOffsetHexagonalNeighbourhoodWithFlow2D()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";
        
        // Set network in filename
        std::string network_name = "TestHexagonalNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("hexagonal_perfusion_quotients.txt");
        outfile.close();

        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Define the key parameters
        // double dimless_domain_size_x = 2050.0;  // x-coordinate of output node
        // double dimless_domain_size_y = 1818.65 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
        // unsigned dimless_vessel_length = 100.0;
        unsigned n_vessels = 386;  // number of non-inlet/outlet vessels from which to select ones to make thin
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        double initial_haematocrit = 0.45;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        double percToKill = 0.2;  // percentage of vessels to kill
        unsigned ToBeKilled = (unsigned)(percToKill*n_vessels);  // number to kill

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=1; h_solver++)
        {
            // Initialise the simulation
            std::shared_ptr<VesselNetwork<2> > p_network;
            std::vector<std::shared_ptr<Vessel<2> > > vessels;
            string line2;
            VesselPtr<2> p_current_vessel;
            unsigned ToBeThin;
            unsigned NThin;
            QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
            QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
            // std::ofstream outfileMean;
            std::vector<std::vector<double> > QuotientMatrixMean(NHet,std::vector<double>(ToBeKilled+1));
            std::vector<std::vector<double> > QuotientMatrixAggregate(NHet,std::vector<double>(ToBeKilled+1,0.0));
            std::vector<std::vector<unsigned> > Order;
            double PerfusionQuotient;

            // Read the network layout from a file
            VesselNetworkGenerator<2> network_generator;
            std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix_Offset.txt");
            std::vector<std::vector<double> > rEdgesMatrix;
            string line;
            
            // Break down the rows in the file into column values
            while (std::getline(in, line)) 
            {
                rEdgesMatrix.push_back(std::vector<double>());
                std::stringstream split(line);
                double value;
                while (split >> value)
                {
                    rEdgesMatrix.back().push_back(value);
                }
            }

            // Loop over different thin vessel layouts
            for(unsigned att=0; att<thin_selections; att++)
            {
                // Choose a layout for the thin vessels
                // outfileMean.open("HexagonalFromVoronoiMean.txt");
                std::ifstream in2("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/Selection"+to_string(att+1)+".txt");
                Order.clear();

                // Break down the row into column values
                while (std::getline(in2, line2)) 
                {
                    Order.push_back(std::vector<unsigned>());
                    std::stringstream split2(line2);
                    unsigned value2;
                    while (split2 >> value2)
                    {
                        Order.back().push_back(value2);
                    }
                }

                // Loop over increasing level of heteregeneity (% of thin vessels)
                for(unsigned i=0; i<NHet;i++)
                {
                    // Generate the network
                    p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);

                    // Set inlet and outlet nodes (if node (start/end) is at the edge of the domain space, make it input or output depending on which side)
                    auto p_segment = p_network->GetVesselSegments()[0];
                    vessels = p_network->GetVessels();
                    for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
                    {
                        if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)
                        {
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
                            }
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
                            }
                        }
                        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)   	
                        {
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);
                            }
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);
                            }
                        }
                    }
                    p_segment = p_network->GetVesselSegments()[0];
                    p_segment->SetRadius(large_vessel_radius);
                    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   

                    // Make a certain % of vessels thin
                    double percOfThin = 0.05*i;
                    ToBeThin = (unsigned)(percOfThin*n_vessels);
                    NThin = 0;                    
                    while( NThin < ToBeThin)
                    {
                        p_current_vessel = p_network->GetVessel(Order[NThin][0]-1);  // ??? '-1' accounts for vessel ID change
                        p_current_vessel->SetRadius(small_vessel_radius);
                        NThin++;
                    }
                    vessels = p_network->GetVessels();

                    // Set the haematocrit solver
                    std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                    std::string solver_name;
                    if (h_solver==1)
                    {
                        solver_name = "ConstantHaematocrit"; 
                        std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==2)
                    {
                        solver_name = "PriesHaematocrit"; 
                        std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                
                    }
                    else if (h_solver==3)
                    {
                        solver_name = "MemoryHaematocrit"; 
                        std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==4)
                    {
                        solver_name = "FungHaematocrit"; 
                        std::cout << "Now using FungHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==5)
                    {
                        solver_name = "YangHaematocrit"; 
                        std::cout << "Now using YangHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();   
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }

                    // Set up the viscosity solver
                    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                    p_viscosity_calculator->SetVesselNetwork(p_network);
                    p_viscosity_calculator->Calculate();
                    
                    // Set up the impedance solver
                    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                    p_impedance_calculator->SetVesselNetwork(p_network);
                    p_impedance_calculator->Calculate();

                    // Set up the flow solver 
                    FlowSolver<2> flow_solver;
                    flow_solver.SetVesselNetwork(p_network);
                    // flow_solver.SetUp();
                    flow_solver.SetUseDirectSolver(true);
                    // flow_solver.Solve();

                    // Set up the grid for the finite difference solver
                    auto p_grid = RegularGrid<2>::Create();
                    p_grid->SetSpacing(grid_spacing);
                    c_vector<unsigned, 3> dimensions;
                    dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                    dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                    dimensions[2] = 1;
                    p_grid->SetDimensions(dimensions);

                    // Choose the PDE
                    auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                    
                    // Set the diffusivity and decay terms
                    p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                    p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                    // Set up the discrete source
                    auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                    QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                            GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                    QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                            Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                    p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                    p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                    p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                    // Set up the finite difference solver for oxygen (which handles everything)
                    auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                    p_oxygen_solver->SetPde(p_oxygen_pde);
                    p_oxygen_solver->SetLabel("oxygen");
                    p_oxygen_solver->SetGrid(p_grid);
                    // solver.SetFileName("oxygen_solution_0");

                    // Prune all vessels up to specified dose          
                    for(unsigned KilledVessels=0; KilledVessels < ToBeKilled+1; KilledVessels++)
                    {
                        // Display status message
                        std::cout << "Now killed " << KilledVessels << " vessels." << std::endl;

                        // Write snapshot
                        // unsigned PictureNumber = (unsigned)(0.5*ToBeKilled);
                        // if(KilledVessels==PictureNumber)
                        // {
                        //     p_network->Write(output_file_RT);
                        // }
                        // else if(KilledVessels==0)
                        // {
                        //     p_network->Write(output_file_two_radii);
                        // }
                        // p_network->Write(output_file_two_radii);

                        //p_network->RemoveVessel(p_network->GetVessel(Order[KilledVessels][0]-1),true);
                    
                        // Set up flag for broken solver
                        unsigned broken_solver = 0;

                        // Set filename
                        std::stringstream selection_stream;
                        std::stringstream heterogeneity_stream;
                        std::stringstream kill_stream;
                        selection_stream << std::fixed << std::setprecision(0) << to_string(att+1);                    heterogeneity_stream << std::fixed << std::setprecision(0) << i;                              kill_stream << std::fixed << std::setprecision(0) << KilledVessels;           
                        std::string selection_string = selection_stream.str();
                        std::string heterogeneity_string = heterogeneity_stream.str();                    
                        std::string kill_string = kill_stream.str();                    
                        std::string file_name = network_name + solver_name + "/Selection" + selection_string + "/NHet" + heterogeneity_string + "/Kills" + kill_string;
                        std::string str_directory_name = file_name;
                         std::cout << str_directory_name << std::endl;
                        auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                        // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                        unsigned max_iter = 1000;  // 1000 
                        double tolerance2 = 1.e-10;
                        std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                        std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                        for(unsigned idx=0; idx<max_iter; idx++)
                        {
                            // Run the solvers
                            p_impedance_calculator->Calculate();
                            p_abstract_haematocrit_solver->Calculate();
                            p_viscosity_calculator->Calculate();
                            flow_solver.SetUp();
                            flow_solver.Solve();

                            // Check for convergence 
                            double max_difference = 0.0;
                            double h_for_max = 0.0;
                            double prev_for_max = 0.0;
                            for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                            {
                                double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                                double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                                if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                                {
                                    max_difference = difference;
                                    h_for_max = current_haematocrit;
                                    prev_for_max = previous_haematocrit[jdx];
                                }
                                previous_haematocrit[jdx] = current_haematocrit;
                            }
                            std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                            // Print the final or intermediary convergence results
                            if(max_difference<=tolerance2)  
                            {
                                std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                                broken_solver = 0;
                                break;
                            }
                            else
                            {
                                if(idx%1==0)
                                {
                                    std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                                    // std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                                    // std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                                    // p_network->Write(output_file);
                                }
                            }

                            // If there is no convergence after all the iterations, print the error message.
                            if(idx==max_iter-1)
                            {
                                std::cout << "Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(att) << " and NHet = " << i << std::endl;
                                error_log << "\n Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(att) << " and NHet = " << i; 
                                broken_solver = 1;
                                break;
                                // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                            }
                        }

                        // If solver doesn't converge, move on to next one
                        if (broken_solver == 1)
                        {
                            continue;
                        }

                        // Run the simulation 
                        SimulationTime::Instance()->SetStartTime(0.0);
                        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                        auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                        p_microvessel_solver->SetVesselNetwork(p_network);
                        p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                        p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                        p_microvessel_solver->Run();
                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                        p_network->Write(output_file);

                        // Print the average oxygenation
                        std::vector<double> solution = p_oxygen_solver->GetSolution();
                        double average_oxygen = 0.0;
                        for(unsigned jdx=0;jdx<solution.size();jdx++)
                        {
                            average_oxygen += solution[jdx];
                        }
                        average_oxygen /= double(solution.size());
                        std::cout << "Average oxygen: " << average_oxygen << std::endl;

                        // Write the PQs
                        QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;
                        // outfile.open("hexagonal_perfusion_quotients.txt", std::ios_base::app);
                        // outfile << network_name << " " << solver_name << " " << selection_string << " " << heterogeneity_string << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                        // outfile.close();  
                        PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
                        // QuotientMatrixAggregate[i][KilledVessels] += PerfusionQuotient;  // store the PQ with the associated NHet and number of kills
                        outfile.open("hexagonal_perfusion_quotients.txt", std::ios_base::app);
                        outfile << network_name << " " << solver_name << " " << heterogeneity_string << " " << selection_string << " " << kill_string << " " << PerfusionQuotient << " \n"; 
                        outfile.close();

                        // Remove the vessel
                        p_network->RemoveVessel(vessels[Order[KilledVessels][0]-1],true);  // remove the vessel ('-1' accounts for vessel ID change after pruning)

                        // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                        ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                        ParameterCollection::Instance()->Destroy();
                        BaseUnits::Instance()->Destroy();
                        SimulationTime::Instance()->Destroy();
                    }
                
                    // std::cout << "Attempts " << att+1 << " \n";
                    // outfileMean << "Attempts " << att+1 << " \n";
                    // Write the resulting perfusion quotients into an output file
                    // if (i==NHet-1)
                    // {   
                    //     outfile.open("hexagonal_perfusion_quotients.txt");
                    //     for (unsigned iii = 0; iii < NHet; iii++)
                    //     {
                    //         for(unsigned jjj = 0; jjj < ToBeKilled+1; jjj++)
                    //         {
                    //             QuotientMatrixMean[iii][jjj]= QuotientMatrixAggregate[iii][jjj]/double(att+1);
                    //             //std::cout << "For heterogeneity level " << iii << ", and radiotherapy dose " << jjj << ", the aggregate quotient is: " << QuotientMatrixAggregate[iii][jjj] << " , and the mean quotient is: " << QuotientMatrixMean[iii][jjj] << endl;
                    //             outfile << network_name << " " << solver_name << " " << iii << " " << jjj <<  " " << QuotientMatrixMean[iii][jjj] << " \n"; 
                    //             // outfile.close();
                    //             // outfileMean << iii << " " << jjj << " " << QuotientMatrixMean[iii][jjj] << " \n";
                    //         }
                    //     }    
                    //     outfile.close();
                    // }
                }
            }
        }
         
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a full 2D hexagonal network on a PDE grid with flow, H-splitting, and radii set according to a log normal distribution and pruning based on radius thresholds (correcting for the unexplained -50 um offset in the other simulation)
    void xTestOffsetHeterogeneousHexagonalNeighbourhoodWithRadiusThresholdPruning2D()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";
        
        // Set network in filename
        std::string network_name = "TestHexagonalNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("hex_lognormal_radius_threshold_perfusion_quotients.txt");
        outfile.close();

        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Define the key parameters
        unsigned n_vessels = 386;  // number of non-inlet/outlet vessels from which to select ones to make thin
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        double initial_haematocrit = 0.45;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        double percToKill = 0.2;  // percentage of vessels to kill
        unsigned ToBeKilled = (unsigned)(percToKill*n_vessels);  // number to kill
        unsigned max_n_alpha = 3;

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=1; h_solver++)
        {
            // Initialise the simulation
            std::shared_ptr<VesselNetwork<2> > p_network;
            std::vector<std::shared_ptr<Vessel<2> > > vessels;
            string line2;
            VesselPtr<2> p_current_vessel;
            QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
            QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
            std::vector<std::vector<double> > QuotientMatrixMean(max_n_alpha,std::vector<double>(ToBeKilled+1));
            std::vector<std::vector<double> > QuotientMatrixAggregate(max_n_alpha,std::vector<double>(ToBeKilled+1,0.0));
            std::vector<std::vector<unsigned> > Order;
            double PerfusionQuotient;

            // Read the network layout from a file and store it in an array
            VesselNetworkGenerator<2> network_generator;
            std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix_Offset.txt");
            std::vector<std::vector<double> > rEdgesMatrix;
            string line;
            while (std::getline(in, line)) 
            {
                rEdgesMatrix.push_back(std::vector<double>());
                std::stringstream split(line);
                double value;
                while (split >> value)
                {
                    rEdgesMatrix.back().push_back(value);
                }
            }

            // Loop over increasing level of heteregeneity
            for(unsigned n_alpha=0; n_alpha<=max_n_alpha; n_alpha++)
            {
                // Convert to actual mu
                int alpha = (3*n_alpha)+10;  // alpha determines the mean
                
                // Loop over different thin vessel layouts
                for(unsigned list_number=1; list_number<=thin_selections; list_number++)
                {                
                    // Choose a radii list file and read the file to an array
                    std::vector<std::vector<double> > radii_array;
                    radii_array.clear();
                    string line;
                    std::ifstream radius_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/hexagonal_radius_log_normal_distribution/mu_"+to_string(alpha)+"/radii_list_"+to_string(list_number)+".txt");
                    while (std::getline(radius_list_file, line)) 
                    {
                        radii_array.push_back(std::vector<double>());
                        std::stringstream split(line);
                        double value;
                        while (split >> value)
                        {
                            radii_array.back().push_back(value);
                        }
                    }

                    // Choose the corresponding vessel ID order file and read the file to an array
                    std::vector<std::vector<double> > id_array;
                    id_array.clear();
                    string line_2;
                    std::ifstream id_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/hexagonal_radius_log_normal_distribution/mu_"+to_string(alpha)+"/id_list_"+to_string(list_number)+".txt");
                    while (std::getline(id_list_file, line_2)) 
                    {
                        id_array.push_back(std::vector<double>());
                        std::stringstream split(line_2);
                        double value;
                        while (split >> value)
                        {
                            id_array.back().push_back(value);
                        }
                    }

                    // Generate the network
                    p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);

                    // Set inlet and outlet nodes (if node (start/end) is at the edge of the domain space, make it input or output depending on which side), and assign each non-inlet/outlet vessel a radius from the list
                    auto p_segment = p_network->GetVesselSegments()[0];
                    p_segment->SetRadius(inlet_vessel_radius);
                    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   
                    vessels = p_network->GetVessels();
                    for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
                    {
                        if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the start node equals one (no other segments at nodes)
                        {
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the start node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);  // set the node as the input node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the start node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the ouput node pressure
                            }
                        }
                        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the end node equals one (no other segments at nodes)	
                        {
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the end node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);;  // set the node as the input node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the end node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the output node pressure
                            }
                        }
                    }  

                    // Assign vessel radii based on imported arrays
                    unsigned int vessel_index = 0;
                    while(vessel_index < n_vessels)
                    {
                        p_current_vessel = p_network->GetVessel(id_array[vessel_index][0]);

                        // Get segment radius from array
                        QLength new_vessel_radius(double(radii_array[vessel_index][0])*unit::microns);

                        p_current_vessel->SetRadius(new_vessel_radius);
                        vessel_index++;
                    }
                    vessels = p_network->GetVessels();

                    // Set the haematocrit solver
                    std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                    std::string solver_name;
                    if (h_solver==1)
                    {
                        solver_name = "ConstantHaematocrit"; 
                        std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==2)
                    {
                        solver_name = "PriesHaematocrit"; 
                        std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                
                    }
                    else if (h_solver==3)
                    {
                        solver_name = "MemoryHaematocrit"; 
                        std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==4)
                    {
                        solver_name = "FungHaematocrit"; 
                        std::cout << "Now using FungHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==5)
                    {
                        solver_name = "YangHaematocrit"; 
                        std::cout << "Now using YangHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();   
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }

                    // Set up the viscosity solver
                    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                    p_viscosity_calculator->SetVesselNetwork(p_network);
                    p_viscosity_calculator->Calculate();
                    
                    // Set up the impedance solver
                    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                    p_impedance_calculator->SetVesselNetwork(p_network);
                    p_impedance_calculator->Calculate();

                    // Set up the flow solver 
                    FlowSolver<2> flow_solver;
                    flow_solver.SetVesselNetwork(p_network);
                    flow_solver.SetUseDirectSolver(true);

                    // Set up the grid for the finite difference solver
                    // auto p_grid = RegularGrid<2>::Create();
                    // p_grid->SetSpacing(grid_spacing);
                    // c_vector<unsigned, 3> dimensions;
                    // dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                    // dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                    // dimensions[2] = 1;
                    // p_grid->SetDimensions(dimensions);

                    // Choose the PDE
                    // auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                    
                    // Set the diffusivity and decay terms
                    // p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                    // p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                    // Set up the discrete source
                    // auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                    // QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                    //         GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                    // QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                    //         Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                    // p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                    // p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                    // p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    // p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                    // Set up the finite difference solver for oxygen (which handles everything)
                    // auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                    // p_oxygen_solver->SetPde(p_oxygen_pde);
                    // p_oxygen_solver->SetLabel("oxygen");
                    // p_oxygen_solver->SetGrid(p_grid);

                    // Set up pruning parameters
                    QLength radius_ceiling = 21_um;  // set threshold up to and excluding which vessels should be pruned
                    QLength radius_step = 1_um;
                    QLength current_radius_threshold = 0_um;

                    // Prunes vessels from thinnest up to simulate RT
                    while (current_radius_threshold <= radius_ceiling)
                    { 
                        // Display status message
                        double current_radius_threshold_um = current_radius_threshold*1000000;  // in micrometres
                        std::cout << "Now pruning up to vessels with radius = " << current_radius_threshold_um << " um" << std::endl;

                        // Remove vessels under radius threshold
                        p_network->RemoveThinVessels(current_radius_threshold, false);  

                        // Set filename
                        std::stringstream selection_stream;
                        std::stringstream heterogeneity_stream;
                        std::stringstream threshold_stream;
                        selection_stream << std::fixed << std::setprecision(0) << to_string(list_number);                    
                        heterogeneity_stream << std::fixed << std::setprecision(0) << alpha;                              
                        threshold_stream << std::fixed << std::setprecision(0) << current_radius_threshold_um;           
                        std::string selection_string = selection_stream.str();
                        std::string heterogeneity_string = heterogeneity_stream.str();                    
                        std::string threshold_string = threshold_stream.str();                    
                        std::string file_name = network_name + solver_name + "/Mu" + heterogeneity_string + "/Selection" + selection_string + "/RadiusThreshold" + threshold_string;
                        std::string str_directory_name = file_name;
                        std::cout << str_directory_name << std::endl;
                        auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                        // Set up flag for broken solver
                        unsigned broken_solver = 0;

                        // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                        unsigned max_iter = 5;  // 1000 
                        double tolerance2 = 1.e-10;
                        std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                        std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                        for(unsigned idx=0; idx<max_iter; idx++)
                        {
                            // Run the solvers
                            p_impedance_calculator->Calculate();
                            p_abstract_haematocrit_solver->Calculate();
                            p_viscosity_calculator->Calculate();
                            flow_solver.SetUp();
                            flow_solver.Solve();

                            // Get the residual
                            double max_difference = 0.0;
                            double h_for_max = 0.0;
                            double prev_for_max = 0.0;
                            for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                            {
                                // Set segments with no flow to have no haematocrit
                                if (fabs(segments[jdx]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                {
                                    segments[jdx]->GetFlowProperties()->SetHaematocrit(0.0);
                                }   
                                double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                                double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                                if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                                {
                                    max_difference = difference;
                                    h_for_max = current_haematocrit;
                                    prev_for_max = previous_haematocrit[jdx];
                                }
                                previous_haematocrit[jdx] = current_haematocrit;
                            }
                            std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                            // Print the final or intermediary convergence results
                            if(max_difference<=tolerance2)  
                            {
                                std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                                broken_solver = 0;
                                break;
                            }
                            else
                            {
                                if(idx%1==0)
                                {
                                    std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                                    // std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                                    // std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                                    // p_network->Write(output_file);
                                }
                            }

                            // If there is no convergence after all the iterations, print the error message.
                            if(idx==max_iter-1)
                            {
                                std::cout << "Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha << std::endl;
                                error_log << "\n Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha; 
                                broken_solver = 1;
                                break;
                                // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                            }
                        }

                        // If solver doesn't converge, move on to next one
                        if (broken_solver == 1)
                        {
                            // Move onto next radius threshold for next step
                            current_radius_threshold = current_radius_threshold + radius_step;  
                            continue;
                        }

                        // Run the simulation 
                        // SimulationTime::Instance()->SetStartTime(0.0);
                        // SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                        // auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                        // p_microvessel_solver->SetVesselNetwork(p_network);
                        // p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                        // p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                        // p_microvessel_solver->Run();
                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                        p_network->Write(output_file);

                        // // Print the average oxygenation
                        // std::vector<double> solution = p_oxygen_solver->GetSolution();
                        // double average_oxygen = 0.0;
                        // for(unsigned jdx=0;jdx<solution.size();jdx++)
                        // {
                        //     average_oxygen += solution[jdx];
                        // }
                        // average_oxygen /= double(solution.size());
                        // std::cout << "Average oxygen: " << average_oxygen << std::endl;

                        // Move onto next radius threshold for next step
                        current_radius_threshold = current_radius_threshold + radius_step;  

                        // // Write the PQs
                        QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;
                        PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
                        // QuotientMatrixAggregate[i][KilledVessels] += PerfusionQuotient;  // store the PQ with the associated NHet and number of kills
                        outfile.open("hex_lognormal_radius_threshold_perfusion_quotients.txt", std::ios_base::app);
                        outfile << network_name << " " << solver_name << " " << heterogeneity_string << " " << selection_string << " " << threshold_string << " " << PerfusionQuotient << " \n"; 
                        outfile.close();

                        // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                        ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                        ParameterCollection::Instance()->Destroy();
                        BaseUnits::Instance()->Destroy();
                        SimulationTime::Instance()->Destroy();
                    }
                
                    // std::cout << "Attempts " << att+1 << " \n";
                    // outfileMean << "Attempts " << att+1 << " \n";
                    // Write the resulting perfusion quotients into an output file
                    // if (alpha==max_alpha)
                    // {   
                    //     outfile.open("hexagonal_perfusion_quotients.txt");
                    //     for (unsigned iii = 0; iii < max_alpha; iii++)
                    //     {
                    //         for(unsigned jjj = 0; jjj < ToBeKilled+1; jjj++)
                    //         {
                    //             QuotientMatrixMean[iii][jjj]= QuotientMatrixAggregate[iii][jjj]/double(list_number+1);
                    //             //std::cout << "For heterogeneity level " << iii << ", and radiotherapy dose " << jjj << ", the aggregate quotient is: " << QuotientMatrixAggregate[iii][jjj] << " , and the mean quotient is: " << QuotientMatrixMean[iii][jjj] << endl;
                    //             outfile << network_name << " " << solver_name << " " << iii << " " << jjj <<  " " << QuotientMatrixMean[iii][jjj] << " \n"; 
                    //             // outfile.close();
                    //             // outfileMean << iii << " " << jjj << " " << QuotientMatrixMean[iii][jjj] << " \n";
                    //         }
                    //     }    
                    //     outfile.close();
                    // }
                }
            }
        }
         
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a full 2D hexagonal network on a PDE grid with flow, H-splitting, radii set according to a log normal distribution, and kills-based pruning (correcting for the unexplained -50 um offset in the other simulation)
    void xTestOffsetHeterogeneousHexagonalNeighbourhoodWithIndividualPruning2D()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";
        
        // Set network in filename
        std::string network_name = "TestHexagonalNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("hex_lognormal_individual_perfusion_quotients.txt");
        outfile.close();

        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Define the key parameters
        // double dimless_domain_size_x = 2050.0;  // x-coordinate of output node
        // double dimless_domain_size_y = 1818.65 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
        // unsigned dimless_vessel_length = 100.0;
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        double initial_haematocrit = 0.45;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        unsigned max_n_alpha = 3;
        unsigned n_vessels = 386;  // number of non-inlet/outlet vessels from which to select ones to make thin
        // double percToKill = 0.2;  // percentage of vessels to kill
        // unsigned ToBeKilled = (unsigned)(percToKill*n_vessels);  // number to kill
        unsigned ToBeKilled = 200;  // number to kill

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=1; h_solver++)
        {
            // Initialise the simulation
            std::shared_ptr<VesselNetwork<2> > p_network;
            std::vector<std::shared_ptr<Vessel<2> > > vessels;
            VesselPtr<2> p_current_vessel;
            QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
            QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
            std::vector<std::vector<double> > QuotientMatrixMean(max_n_alpha,std::vector<double>(ToBeKilled+1));
            std::vector<std::vector<double> > QuotientMatrixAggregate(max_n_alpha,std::vector<double>(ToBeKilled+1,0.0));
            std::vector<std::vector<unsigned> > Order;
            double PerfusionQuotient;

            // Read the network layout from a file and store it in an array
            VesselNetworkGenerator<2> network_generator;
            std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix_Offset.txt");
            std::vector<std::vector<double> > rEdgesMatrix;
            string line;
            while (std::getline(in, line)) 
            {
                rEdgesMatrix.push_back(std::vector<double>());
                std::stringstream split(line);
                double value;
                while (split >> value)
                {
                    rEdgesMatrix.back().push_back(value);
                }
            }

            // Loop over increasing level of heteregeneity (% of thin vessels)
            for(unsigned n_alpha=0; n_alpha<=max_n_alpha; n_alpha++)
            {
                // Convert to actual mu
                int alpha = (3*n_alpha)+10;  // alpha determines the mean

                // Loop over different thin vessel layouts
                for(unsigned list_number=1; list_number<=thin_selections; list_number++)
                {                
                    // Choose a radii list file and read the file to an array
                    std::vector<std::vector<double> > radii_array;
                    radii_array.clear();
                    // int list_index=0;
                    string line_1;
                    std::ifstream radius_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/hexagonal_radius_log_normal_distribution/mu_"+to_string(alpha)+"/radii_list_"+to_string(list_number)+".txt");
                    while (std::getline(radius_list_file, line_1)) 
                    {
                        radii_array.push_back(std::vector<double>());
                        std::stringstream split(line_1);
                        double value;
                        while (split >> value)
                        {
                            radii_array.back().push_back(value);
                        }
                    }

                    // Choose the corresponding vessel ID order file and read the file to an array
                    std::vector<std::vector<double> > id_array;
                    id_array.clear();
                    // int list_index=0;
                    string line_2;
                    std::ifstream id_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/hexagonal_radius_log_normal_distribution/mu_"+to_string(alpha)+"/id_list_"+to_string(list_number)+".txt");
                    while (std::getline(id_list_file, line_2)) 
                    {
                        id_array.push_back(std::vector<double>());
                        std::stringstream split(line_2);
                        double value;
                        while (split >> value)
                        {
                            id_array.back().push_back(value);
                        }
                    }

                    // Generate the network
                    p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);

                    // Set inlet and outlet nodes (if node (start/end) is at the edge of the domain space, make it input or output depending on which side), and assign each non-inlet/outlet vessel a radius from the list
                    auto p_segment = p_network->GetVesselSegments()[0];
                    p_segment->SetRadius(inlet_vessel_radius);
                    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   
                    vessels = p_network->GetVessels();
                    for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
                    {
                        if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the start node equals one (no other segments at nodes)
                        {
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the start node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);  // set the node as the input node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the start node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the ouput node pressure
                            }
                        }
                        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the end node equals one (no other segments at nodes)	
                        {
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the end node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);;  // set the node as the input node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the end node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the output node pressure
                            }
                        }
                    }
                
                    // Assign vessel radii based on imported arrays
                    unsigned int vessel_index = 0;
                    while(vessel_index < n_vessels)
                    {
                        p_current_vessel = p_network->GetVessel(id_array[vessel_index][0]);

                        // Get segment radius from array
                        QLength new_vessel_radius(double(radii_array[vessel_index][0])*unit::microns);

                        p_current_vessel->SetRadius(new_vessel_radius);
                        vessel_index++;
                    }
                    vessels = p_network->GetVessels();

                    // Set the haematocrit solver
                    std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                    std::string solver_name;
                    if (h_solver==1)
                    {
                        solver_name = "ConstantHaematocrit"; 
                        std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==2)
                    {
                        solver_name = "PriesHaematocrit"; 
                        std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                
                    }
                    else if (h_solver==3)
                    {
                        solver_name = "MemoryHaematocrit"; 
                        std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==4)
                    {
                        solver_name = "FungHaematocrit"; 
                        std::cout << "Now using FungHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==5)
                    {
                        solver_name = "YangHaematocrit"; 
                        std::cout << "Now using YangHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();   
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }

                    // Set up the viscosity solver
                    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                    p_viscosity_calculator->SetVesselNetwork(p_network);
                    p_viscosity_calculator->Calculate();
                    
                    // Set up the impedance solver
                    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                    p_impedance_calculator->SetVesselNetwork(p_network);
                    p_impedance_calculator->Calculate();

                    // Set up the flow solver 
                    FlowSolver<2> flow_solver;
                    flow_solver.SetVesselNetwork(p_network);
                    // flow_solver.SetUp();
                    flow_solver.SetUseDirectSolver(true);
                    // flow_solver.Solve();

                    // Set up the grid for the finite difference solver
                    // auto p_grid = RegularGrid<2>::Create();
                    // p_grid->SetSpacing(grid_spacing);
                    // c_vector<unsigned, 3> dimensions;
                    // dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                    // dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                    // dimensions[2] = 1;
                    // p_grid->SetDimensions(dimensions);

                    // // Choose the PDE
                    // auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                    
                    // // Set the diffusivity and decay terms
                    // p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                    // p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                    // // Set up the discrete source
                    // auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                    // QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                    //         GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                    // QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                    //         Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                    // p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                    // p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                    // p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    // p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                    // // Set up the finite difference solver for oxygen (which handles everything)
                    // auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                    // p_oxygen_solver->SetPde(p_oxygen_pde);
                    // p_oxygen_solver->SetLabel("oxygen");
                    // p_oxygen_solver->SetGrid(p_grid);
                    // solver.SetFileName("oxygen_solution_0");

                    // Prune all vessels up to specified dose 
                    for(unsigned KilledVessels=0; KilledVessels < ToBeKilled+1; KilledVessels++)
                    { 
                        // Display status message
                        std::cout << "Now killed " << KilledVessels << " vessels." << std::endl;  

                        // Set filename
                        std::stringstream selection_stream;
                        std::stringstream heterogeneity_stream;
                        std::stringstream kill_stream;
                        selection_stream << std::fixed << std::setprecision(0) << to_string(list_number);                    
                        heterogeneity_stream << std::fixed << std::setprecision(0) << alpha;                              
                        kill_stream << std::fixed << std::setprecision(0) << KilledVessels;           
                        std::string selection_string = selection_stream.str();
                        std::string heterogeneity_string = heterogeneity_stream.str();                    
                        std::string kill_string = kill_stream.str();                    
                        std::string file_name = network_name + solver_name + "/Mu" + heterogeneity_string + "/Selection" + selection_string + "/Kills" + kill_string;
                        std::string str_directory_name = file_name;
                        std::cout << str_directory_name << std::endl;
                        auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                        // Set up flag for broken solver
                        unsigned broken_solver = 0;

                        // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                        unsigned max_iter = 5;  // 1000 
                        double tolerance2 = 1.e-10;
                        std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                        std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                        for(unsigned idx=0; idx<max_iter; idx++)
                        {
                            // Run the solvers
                            p_impedance_calculator->Calculate();
                            p_abstract_haematocrit_solver->Calculate();
                            p_viscosity_calculator->Calculate();
                            flow_solver.SetUp();
                            flow_solver.Solve();

                            // Get the residual
                            double max_difference = 0.0;
                            double h_for_max = 0.0;
                            double prev_for_max = 0.0;
                            for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                            {
                                // Set segments with no flow to have no haematocrit
                                if (fabs(segments[jdx]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                {
                                    segments[jdx]->GetFlowProperties()->SetHaematocrit(0.0);
                                }   
                                double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                                double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                                if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                                {
                                    max_difference = difference;
                                    h_for_max = current_haematocrit;
                                    prev_for_max = previous_haematocrit[jdx];
                                }
                                previous_haematocrit[jdx] = current_haematocrit;
                            }
                            std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                            // Print the final or intermediary convergence results
                            if(max_difference<=tolerance2)  
                            {
                                std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                                broken_solver = 0;
                                break;
                            }
                            else
                            {
                                if(idx%1==0)
                                {
                                    std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                                    std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                                    std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                                    p_network->Write(output_file);
                                }
                            }

                            // If there is no convergence after all the iterations, print the error message.
                            if(idx==max_iter-1)
                            {
                                std::cout << "Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha << std::endl;
                                error_log << "\n Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha; 
                                broken_solver = 1;
                                break;
                                // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                            }
                        }

                        // If solver doesn't converge, move on to next one
                        if (broken_solver == 1)
                        {
                            continue;
                        }

                        // Run the simulation 
                        // SimulationTime::Instance()->SetStartTime(0.0);
                        // SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                        // auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                        // p_microvessel_solver->SetVesselNetwork(p_network);
                        // p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                        // p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                        // p_microvessel_solver->Run();
                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                        p_network->Write(output_file);

                        // Print the average oxygenation
                        // std::vector<double> solution = p_oxygen_solver->GetSolution();
                        // double average_oxygen = 0.0;
                        // for(unsigned jdx=0;jdx<solution.size();jdx++)
                        // {
                        //     average_oxygen += solution[jdx];
                        // }
                        // average_oxygen /= double(solution.size());
                        // std::cout << "Average oxygen: " << average_oxygen << std::endl;

                        // // Write the PQs
                        QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;
                        // // outfile.open("hexagonal_perfusion_quotients.txt", std::ios_base::app);
                        // // outfile << network_name << " " << solver_name << " " << selection_string << " " << heterogeneity_string << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                        // // outfile.close();  
                        PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
                        // QuotientMatrixAggregate[i][KilledVessels] += PerfusionQuotient;  // store the PQ with the associated NHet and number of kills
                        outfile.open("hex_lognormal_individual_perfusion_quotients.txt", std::ios_base::app);
                        outfile << network_name << " " << solver_name << " " << heterogeneity_string << " " << selection_string << " " << kill_string << " " << PerfusionQuotient << " \n"; 
                        outfile.close();

                        // Remove the smallest vessel
                        QLength minimum_radius = inlet_vessel_radius;
                        unsigned int minimum_index = 0;
                        for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)  // for all the segments in the network
                        {
                            // Get the current segment's radius
                            QLength current_radius = vessels[vessel_index]->GetRadius();

                            // If the current radius is less than the minimum radius, record the new minimum
                            if (current_radius < minimum_radius)
                            {
                                minimum_radius = current_radius;
                                minimum_index = vessel_index;
                            }
                        }
                        p_network->RemoveVessel(vessels[minimum_index], true);  // remove the vessel

                        // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                        ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                        ParameterCollection::Instance()->Destroy();
                        BaseUnits::Instance()->Destroy();
                        SimulationTime::Instance()->Destroy();
                    }
                }
            }
        }
         
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a full 2D hexagonal network on a PDE grid with flow, H-splitting, radii set according to a log normal distribution, and stochastic pruning (correcting for the unexplained -50 um offset in the other simulation)
    void xTestOffsetHeterogeneousHexagonalNeighbourhoodWithStochasticPruning2D()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";
        
        // Set network in filename
        std::string network_name = "TestHexagonalNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("hex_lognormal_stochastic_perfusion_quotients.txt");
        outfile.close();

        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Define the key parameters
        // double dimless_domain_size_x = 2050.0;  // x-coordinate of output node
        // double dimless_domain_size_y = 1818.65 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
        // unsigned dimless_vessel_length = 100.0;
        unsigned n_vessels = 386;  // number of non-inlet/outlet vessels from which to select ones to make thin
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        double initial_haematocrit = 0.45;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        double percToKill = 0.2;  // percentage of vessels to kill
        unsigned ToBeKilled = (unsigned)(percToKill*n_vessels);  // number to kill
        max_alpha = 4;

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=1; h_solver++)
        {
            // Initialise the simulation
            std::shared_ptr<VesselNetwork<2> > p_network;
            std::vector<std::shared_ptr<Vessel<2> > > vessels;
            string line2;
            VesselPtr<2> p_current_vessel;
            // unsigned ToBeThin;
            // unsigned NThin;
            QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
            QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
            // std::ofstream outfileMean;
            std::vector<std::vector<double> > QuotientMatrixMean(max_alpha,std::vector<double>(ToBeKilled+1));
            std::vector<std::vector<double> > QuotientMatrixAggregate(max_alpha,std::vector<double>(ToBeKilled+1,0.0));
            std::vector<std::vector<unsigned> > Order;
            double PerfusionQuotient;

            // Read the network layout from a file and store it in an array
            VesselNetworkGenerator<2> network_generator;
            std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix_Offset.txt");
            std::vector<std::vector<double> > rEdgesMatrix;
            string line;
            while (std::getline(in, line)) 
            {
                rEdgesMatrix.push_back(std::vector<double>());
                std::stringstream split(line);
                double value;
                while (split >> value)
                {
                    rEdgesMatrix.back().push_back(value);
                }
            }

            // Loop over increasing level of heteregeneity (% of thin vessels)
            for(unsigned alpha=6; alpha<=8; alpha++)
            {
                // Loop over different thin vessel layouts
                for(unsigned list_number=1; list_number<=thin_selections; list_number++)
                {   
                    // Choose a radii list file and read the file to an array
                    std::vector<std::vector<double> > radii_array;
                    radii_array.clear();
                    // int list_index=0;
                    string line;
                    std::ifstream radius_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/hexagonal_radius_log_normal_distribution/sigma_"+to_string(alpha)+"/radii_list_"+to_string(list_number)+".txt");
                    while (std::getline(radius_list_file, line)) 
                    {
                        radii_array.push_back(std::vector<double>());
                        std::stringstream split(line);
                        double value;
                        while (split >> value)
                        {
                            radii_array.back().push_back(value);
                        }
                    }

                    // Choose the corresponding vessel ID order file and read the file to an array
                    std::vector<std::vector<double> > id_array;
                    id_array.clear();
                    // int list_index=0;
                    string line_2;
                    std::ifstream id_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/hexagonal_radius_log_normal_distribution/sigma_"+to_string(alpha)+"/id_list_"+to_string(list_number)+".txt");
                    while (std::getline(id_list_file, line_2)) 
                    {
                        id_array.push_back(std::vector<double>());
                        std::stringstream split(line_2);
                        double value;
                        while (split >> value)
                        {
                            id_array.back().push_back(value);
                        }
                    }

                    // Set the maximum beta value
                    int max_beta = 35;

                    // Iterate over beta values
                    int beta = 1;
                    while (beta <= max_beta)
                    {
                        // Generate the network
                        p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);

                        // Set inlet and outlet nodes (if node (start/end) is at the edge of the domain space, make it input or output depending on which side), and assign each non-inlet/outlet vessel a radius from the list
                        auto p_segment = p_network->GetVesselSegments()[0];
                        p_segment->SetRadius(inlet_vessel_radius);
                        VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   
                        vessels = p_network->GetVessels();
                        for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
                        {
                            if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the start node equals one (no other segments at nodes)
                            {
                                if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the start node is located approx. at the start of the domain 
                                {
                                    (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);  // set the node as the input node
                                    (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                                }
                                if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the start node is located approx. at the end of the domain 
                                {
                                    (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                    (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the ouput node pressure
                                }
                            }
                            if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the end node equals one (no other segments at nodes)	
                            {
                                if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the end node is located approx. at the start of the domain 
                                {
                                    (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);;  // set the node as the input node
                                    (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                                }
                                if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the end node is located approx. at the end of the domain 
                                {
                                    (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                    (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the output node pressure
                                }
                            }
                        }  

                        // Assign vessel radii based on imported arrays
                        unsigned int vessel_index = 0;
                        while(vessel_index < n_vessels)
                        {
                            p_current_vessel = p_network->GetVessel(id_array[vessel_index][0]);

                            // Get segment radius from array
                            QLength new_vessel_radius(double(radii_array[vessel_index][0])*unit::microns);

                            p_current_vessel->SetRadius(new_vessel_radius);
                            vessel_index++;
                        }
                        vessels = p_network->GetVessels();
                        
                        // Set the haematocrit solver
                        std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                        std::string solver_name;
                        if (h_solver==1)
                        {
                            solver_name = "ConstantHaematocrit"; 
                            std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                            auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                            p_haematocrit_solver->SetVesselNetwork(p_network);
                            p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                            p_abstract_haematocrit_solver = p_haematocrit_solver;     
                        }
                        else if (h_solver==2)
                        {
                            solver_name = "PriesHaematocrit"; 
                            std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                            auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                            // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                            p_haematocrit_solver->SetVesselNetwork(p_network);
                            p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                            p_abstract_haematocrit_solver = p_haematocrit_solver;                
                        }
                        else if (h_solver==3)
                        {
                            solver_name = "MemoryHaematocrit"; 
                            std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                            auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                            // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                            p_haematocrit_solver->SetVesselNetwork(p_network);
                            p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                            p_abstract_haematocrit_solver = p_haematocrit_solver;     
                        }
                        else if (h_solver==4)
                        {
                            solver_name = "FungHaematocrit"; 
                            std::cout << "Now using FungHaematocritSolver..." << std::endl;
                            auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                            // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                            p_haematocrit_solver->SetVesselNetwork(p_network);
                            p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                            p_abstract_haematocrit_solver = p_haematocrit_solver;     
                        }
                        else if (h_solver==5)
                        {
                            solver_name = "YangHaematocrit"; 
                            std::cout << "Now using YangHaematocritSolver..." << std::endl;
                            auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                            // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();   
                            p_haematocrit_solver->SetVesselNetwork(p_network);
                            p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                            p_abstract_haematocrit_solver = p_haematocrit_solver;     
                        }

                        // Set up the viscosity solver
                        auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                        p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                        p_viscosity_calculator->SetVesselNetwork(p_network);
                        p_viscosity_calculator->Calculate();
                        
                        // Set up the impedance solver
                        auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                        p_impedance_calculator->SetVesselNetwork(p_network);
                        p_impedance_calculator->Calculate();

                        // Set up the flow solver 
                        FlowSolver<2> flow_solver;
                        flow_solver.SetVesselNetwork(p_network);
                        // flow_solver.SetUp();
                        flow_solver.SetUseDirectSolver(true);
                        // flow_solver.Solve();

                        // Set up the grid for the finite difference solver
                        auto p_grid = RegularGrid<2>::Create();
                        p_grid->SetSpacing(grid_spacing);
                        c_vector<unsigned, 3> dimensions;
                        dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                        dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                        dimensions[2] = 1;
                        p_grid->SetDimensions(dimensions);

                        // Choose the PDE
                        auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                        
                        // Set the diffusivity and decay terms
                        p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                        p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                        // Set up the discrete source
                        auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                        QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                                GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                        QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                                Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                        p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                        p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                        p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                        // Set up the finite difference solver for oxygen (which handles everything)
                        auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                        p_oxygen_solver->SetPde(p_oxygen_pde);
                        p_oxygen_solver->SetLabel("oxygen");
                        p_oxygen_solver->SetGrid(p_grid);
                        // solver.SetFileName("oxygen_solution_0");

                        // Display status message
                        std::cout << "Now pruning with beta = " << beta << std::endl;

                        // Conduct stochastic pruning (thinner vessels are more likely to be pruned)
                        p_network->StochasticPruning(beta, false);  

                        // Set filename
                        std::stringstream selection_stream;
                        std::stringstream heterogeneity_stream;
                        std::stringstream beta_stream;
                        selection_stream << std::fixed << std::setprecision(0) << to_string(list_number);                    
                        heterogeneity_stream << std::fixed << std::setprecision(0) << alpha;                              
                        beta_stream << std::fixed << std::setprecision(0) << beta;
                        std::string selection_string = selection_stream.str();
                        std::string heterogeneity_string = heterogeneity_stream.str();                    
                        std::string beta_string = beta_stream.str();
                        std::string file_name = network_name + solver_name + "/Sigma" + heterogeneity_string + "/Selection" + selection_string + "/Beta" + beta_string;
                        std::string str_directory_name = file_name;
                        std::cout << str_directory_name << std::endl;
                        auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                        // Set up flag for broken solver
                        unsigned broken_solver = 0;

                        // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                        unsigned max_iter = 1000;  // 1000 
                        double tolerance2 = 1.e-10;
                        std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                        std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                        for(unsigned idx=0; idx<max_iter; idx++)
                        {
                            // Run the solvers
                            p_impedance_calculator->Calculate();
                            p_abstract_haematocrit_solver->Calculate();
                            p_viscosity_calculator->Calculate();
                            flow_solver.SetUp();
                            flow_solver.Solve();

                            // Get the residual
                            double max_difference = 0.0;
                            double h_for_max = 0.0;
                            double prev_for_max = 0.0;
                            for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                            {
                                // Set segments with no flow to have no haematocrit
                                if (fabs(segments[jdx]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                {
                                    segments[jdx]->GetFlowProperties()->SetHaematocrit(0.0);
                                }   
                                double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                                double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                                if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                                {
                                    max_difference = difference;
                                    h_for_max = current_haematocrit;
                                    prev_for_max = previous_haematocrit[jdx];
                                }
                                previous_haematocrit[jdx] = current_haematocrit;
                            }
                            std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                            // Print the final or intermediary convergence results
                            if(max_difference<=tolerance2)  
                            {
                                std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                                broken_solver = 0;
                                break;
                            }
                            else
                            {
                                if(idx%1==0)
                                {
                                    std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                                    // std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                                    // std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                                    // p_network->Write(output_file);
                                }
                            }

                            // If there is no convergence after all the iterations, print the error message.
                            if(idx==max_iter-1)
                            {
                                std::cout << "Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha << std::endl;
                                error_log << "\n Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha; 
                                broken_solver = 1;
                                break;
                                // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                            }
                        }

                        // If solver doesn't converge, move on to next one
                        if (broken_solver == 1)
                        {
                            continue;
                        }

                        // Run the simulation 
                        SimulationTime::Instance()->SetStartTime(0.0);
                        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                        auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                        p_microvessel_solver->SetVesselNetwork(p_network);
                        p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                        p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                        p_microvessel_solver->Run();
                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                        p_network->Write(output_file);

                        // Print the average oxygenation
                        std::vector<double> solution = p_oxygen_solver->GetSolution();
                        double average_oxygen = 0.0;
                        for(unsigned jdx=0;jdx<solution.size();jdx++)
                        {
                            average_oxygen += solution[jdx];
                        }
                        average_oxygen /= double(solution.size());
                        std::cout << "Average oxygen: " << average_oxygen << std::endl;

                        // // Write the PQs
                        QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;
                        // // outfile.open("hexagonal_perfusion_quotients.txt", std::ios_base::app);
                        // // outfile << network_name << " " << solver_name << " " << selection_string << " " << heterogeneity_string << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                        // // outfile.close();  
                        PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
                        // QuotientMatrixAggregate[i][KilledVessels] += PerfusionQuotient;  // store the PQ with the associated NHet and number of kills
                        outfile.open("hex_lognormal_stochastic_perfusion_quotients.txt", std::ios_base::app);
                        outfile << network_name << " " << solver_name << " " << heterogeneity_string << " " << selection_string << " " << beta_string << " " << PerfusionQuotient << " \n"; 
                        outfile.close();

                        // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                        ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                        ParameterCollection::Instance()->Destroy();
                        BaseUnits::Instance()->Destroy();
                        SimulationTime::Instance()->Destroy();

                        // Move on to next beta
                        beta = beta + 1;
                    }
                
                    // std::cout << "Attempts " << att+1 << " \n";
                    // outfileMean << "Attempts " << att+1 << " \n";
                    // Write the resulting perfusion quotients into an output file
                    // if (alpha==max_alpha)
                    // {   
                    //     outfile.open("hexagonal_perfusion_quotients.txt");
                    //     for (unsigned iii = 0; iii < max_alpha; iii++)
                    //     {
                    //         for(unsigned jjj = 0; jjj < ToBeKilled+1; jjj++)
                    //         {
                    //             QuotientMatrixMean[iii][jjj]= QuotientMatrixAggregate[iii][jjj]/double(list_number+1);
                    //             //std::cout << "For heterogeneity level " << iii << ", and radiotherapy dose " << jjj << ", the aggregate quotient is: " << QuotientMatrixAggregate[iii][jjj] << " , and the mean quotient is: " << QuotientMatrixMean[iii][jjj] << endl;
                    //             outfile << network_name << " " << solver_name << " " << iii << " " << jjj <<  " " << QuotientMatrixMean[iii][jjj] << " \n"; 
                    //             // outfile.close();
                    //             // outfileMean << iii << " " << jjj << " " << QuotientMatrixMean[iii][jjj] << " \n";
                    //         }
                    //     }    
                    //     outfile.close();
                    // }
                }
            }
        }
         
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a full 2D hexagonal network on a PDE grid with flow, H-splitting, and radii set according to a log normal distribution and pruning based on radius thresholds (correcting for the unexplained -50 um offset in the other simulation with the Constant haematocrit solver
    void xTestConstantOffsetHeterogeneousHexagonalNeighbourhoodWithRadiusThresholdPruning2D()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";
        
        // Set network in filename
        std::string network_name = "TestHexagonalNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("hex_lognormal_radius_threshold_perfusion_quotients.txt");
        outfile.close();

        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Define the key parameters
        unsigned n_vessels = 386;  // number of non-inlet/outlet vessels from which to select ones to make thin
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        // double initial_haematocrit = 0.45;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        double percToKill = 0.2;  // percentage of vessels to kill
        unsigned ToBeKilled = (unsigned)(percToKill*n_vessels);  // number to kill
        unsigned max_n_alpha = 3;

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=1; h_solver++)
        {
            // Initialise the simulation
            std::shared_ptr<VesselNetwork<2> > p_network;
            std::vector<std::shared_ptr<Vessel<2> > > vessels;
            string line2;
            VesselPtr<2> p_current_vessel;
            QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
            QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
            std::vector<std::vector<double> > QuotientMatrixMean(max_n_alpha,std::vector<double>(ToBeKilled+1));
            std::vector<std::vector<double> > QuotientMatrixAggregate(max_n_alpha,std::vector<double>(ToBeKilled+1,0.0));
            std::vector<std::vector<unsigned> > Order;
            double PerfusionQuotient;

            // Read the network layout from a file and store it in an array
            VesselNetworkGenerator<2> network_generator;
            std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix_Offset.txt");
            std::vector<std::vector<double> > rEdgesMatrix;
            string line;
            while (std::getline(in, line)) 
            {
                rEdgesMatrix.push_back(std::vector<double>());
                std::stringstream split(line);
                double value;
                while (split >> value)
                {
                    rEdgesMatrix.back().push_back(value);
                }
            }

            // Loop over increasing level of heteregeneity
            for(unsigned n_alpha=0; n_alpha<=max_n_alpha; n_alpha++)
            {
                // Convert to actual mu
                int alpha = (3*n_alpha)+10;  // alpha determines the mean
                
                // Loop over different thin vessel layouts
                for(unsigned list_number=1; list_number<=thin_selections; list_number++)
                {                
                    // Choose a radii list file and read the file to an array
                    std::vector<std::vector<double> > radii_array;
                    radii_array.clear();
                    string line;
                    std::ifstream radius_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/hexagonal_radius_log_normal_distribution/mu_"+to_string(alpha)+"/radii_list_"+to_string(list_number)+".txt");
                    while (std::getline(radius_list_file, line)) 
                    {
                        radii_array.push_back(std::vector<double>());
                        std::stringstream split(line);
                        double value;
                        while (split >> value)
                        {
                            radii_array.back().push_back(value);
                        }
                    }

                    // Choose the corresponding vessel ID order file and read the file to an array
                    std::vector<std::vector<double> > id_array;
                    id_array.clear();
                    string line_2;
                    std::ifstream id_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/hexagonal_radius_log_normal_distribution/mu_"+to_string(alpha)+"/id_list_"+to_string(list_number)+".txt");
                    while (std::getline(id_list_file, line_2)) 
                    {
                        id_array.push_back(std::vector<double>());
                        std::stringstream split(line_2);
                        double value;
                        while (split >> value)
                        {
                            id_array.back().push_back(value);
                        }
                    }

                    // Generate the network
                    p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);

                    // Set inlet and outlet nodes (if node (start/end) is at the edge of the domain space, make it input or output depending on which side), and assign each non-inlet/outlet vessel a radius from the list
                    auto p_segment = p_network->GetVesselSegments()[0];
                    p_segment->SetRadius(inlet_vessel_radius);
                    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   
                    vessels = p_network->GetVessels();
                    for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
                    {
                        if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the start node equals one (no other segments at nodes)
                        {
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the start node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);  // set the node as the input node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the start node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the ouput node pressure
                            }
                        }
                        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the end node equals one (no other segments at nodes)	
                        {
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the end node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);;  // set the node as the input node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the end node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the output node pressure
                            }
                        }
                    }  

                    // Assign vessel radii based on imported arrays
                    unsigned int vessel_index = 0;
                    while(vessel_index < n_vessels)
                    {
                        p_current_vessel = p_network->GetVessel(id_array[vessel_index][0]);

                        // Get segment radius from array
                        QLength new_vessel_radius(double(radii_array[vessel_index][0])*unit::microns);

                        p_current_vessel->SetRadius(new_vessel_radius);
                        vessel_index++;
                    }
                    vessels = p_network->GetVessels();

                    // Set the haematocrit solver
                    std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                    std::string solver_name;
                    if (h_solver==1)
                    {
                        solver_name = "ConstantHaematocrit"; 
                        std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==2)
                    {
                        solver_name = "PriesHaematocrit"; 
                        std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                
                    }
                    else if (h_solver==3)
                    {
                        solver_name = "MemoryHaematocrit"; 
                        std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==4)
                    {
                        solver_name = "FungHaematocrit"; 
                        std::cout << "Now using FungHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==5)
                    {
                        solver_name = "YangHaematocrit"; 
                        std::cout << "Now using YangHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();   
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }

                    // Set up the viscosity solver
                    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                    p_viscosity_calculator->SetVesselNetwork(p_network);
                    p_viscosity_calculator->Calculate();
                    
                    // Set up the impedance solver
                    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                    p_impedance_calculator->SetVesselNetwork(p_network);
                    p_impedance_calculator->Calculate();

                    // Set up the flow solver 
                    FlowSolver<2> flow_solver;
                    flow_solver.SetVesselNetwork(p_network);
                    flow_solver.SetUseDirectSolver(true);

                    // Set up the grid for the finite difference solver
                    // auto p_grid = RegularGrid<2>::Create();
                    // p_grid->SetSpacing(grid_spacing);
                    // c_vector<unsigned, 3> dimensions;
                    // dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                    // dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                    // dimensions[2] = 1;
                    // p_grid->SetDimensions(dimensions);

                    // Choose the PDE
                    // auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                    
                    // Set the diffusivity and decay terms
                    // p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                    // p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                    // Set up the discrete source
                    // auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                    // QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                    //         GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                    // QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                    //         Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                    // p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                    // p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                    // p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    // p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                    // Set up the finite difference solver for oxygen (which handles everything)
                    // auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                    // p_oxygen_solver->SetPde(p_oxygen_pde);
                    // p_oxygen_solver->SetLabel("oxygen");
                    // p_oxygen_solver->SetGrid(p_grid);

                    // Set up pruning parameters
                    QLength radius_ceiling = 21_um;  // set threshold up to and excluding which vessels should be pruned
                    QLength radius_step = 1_um;
                    QLength current_radius_threshold = 0_um;

                    // Prunes vessels from thinnest up to simulate RT
                    while (current_radius_threshold <= radius_ceiling)
                    { 
                        // Display status message
                        double current_radius_threshold_um = current_radius_threshold*1000000;  // in micrometres
                        std::cout << "Now pruning up to vessels with radius = " << current_radius_threshold_um << " um" << std::endl;

                        // Remove vessels under radius threshold
                        p_network->RemoveThinVessels(current_radius_threshold, false);  

                        // Set filename
                        std::stringstream selection_stream;
                        std::stringstream heterogeneity_stream;
                        std::stringstream threshold_stream;
                        selection_stream << std::fixed << std::setprecision(0) << to_string(list_number);                    
                        heterogeneity_stream << std::fixed << std::setprecision(0) << alpha;                              
                        threshold_stream << std::fixed << std::setprecision(0) << current_radius_threshold_um;           
                        std::string selection_string = selection_stream.str();
                        std::string heterogeneity_string = heterogeneity_stream.str();                    
                        std::string threshold_string = threshold_stream.str();                    
                        std::string file_name = network_name + solver_name + "/Mu" + heterogeneity_string + "/Selection" + selection_string + "/RadiusThreshold" + threshold_string;
                        std::string str_directory_name = file_name;
                        std::cout << str_directory_name << std::endl;
                        auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                        // Set up flag for broken solver
                        unsigned broken_solver = 0;

                        // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                        // unsigned max_iter = 5;  // 1000 
                        // double tolerance2 = 1.e-10;
                        std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                        // std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                        // for(unsigned idx=0; idx<max_iter; idx++)
                        // {
                            // Run the solvers
                            p_impedance_calculator->Calculate();
                            p_abstract_haematocrit_solver->Calculate();
                            p_viscosity_calculator->Calculate();
                            flow_solver.SetUp();
                            flow_solver.Solve();

                            // // Get the residual
                            // double max_difference = 0.0;
                            // double h_for_max = 0.0;
                            // double prev_for_max = 0.0;
                            for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                            {
                                // Set segments with no flow to have no haematocrit
                                if (fabs(segments[jdx]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                {
                                    segments[jdx]->GetFlowProperties()->SetHaematocrit(0.0);
                                }   
                                // double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                                // double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                                // if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                                // {
                                //     max_difference = difference;
                                //     h_for_max = current_haematocrit;
                                //     prev_for_max = previous_haematocrit[jdx];
                                // }
                                // previous_haematocrit[jdx] = current_haematocrit;
                            }
                            // std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                            // // Print the final or intermediary convergence results
                            // if(max_difference<=tolerance2)  
                            // {
                            //     std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                            //     broken_solver = 0;
                            //     break;
                            // }
                            // else
                            // {
                            //     if(idx%1==0)
                            //     {
                            //         std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                            //         // std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                            //         // std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                            //         // p_network->Write(output_file);
                            //     }
                            // }

                            // If there is no convergence after all the iterations, print the error message.
                            // if(idx==max_iter-1)
                            // {
                            //     std::cout << "Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha << std::endl;
                            //     error_log << "\n Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha; 
                            //     broken_solver = 1;
                            //     break;
                            //     // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                            // }
                        // }

                        // If solver doesn't converge, move on to next one
                        if (broken_solver == 1)
                        {
                            // Move onto next radius threshold for next step
                            current_radius_threshold = current_radius_threshold + radius_step;  
                            continue;
                        }

                        // Run the simulation 
                        // SimulationTime::Instance()->SetStartTime(0.0);
                        // SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                        // auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                        // p_microvessel_solver->SetVesselNetwork(p_network);
                        // p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                        // p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                        // p_microvessel_solver->Run();
                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                        p_network->Write(output_file);

                        // // Print the average oxygenation
                        // std::vector<double> solution = p_oxygen_solver->GetSolution();
                        // double average_oxygen = 0.0;
                        // for(unsigned jdx=0;jdx<solution.size();jdx++)
                        // {
                        //     average_oxygen += solution[jdx];
                        // }
                        // average_oxygen /= double(solution.size());
                        // std::cout << "Average oxygen: " << average_oxygen << std::endl;

                        // Move onto next radius threshold for next step
                        current_radius_threshold = current_radius_threshold + radius_step;  

                        // // Write the PQs
                        QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;
                        PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
                        // QuotientMatrixAggregate[i][KilledVessels] += PerfusionQuotient;  // store the PQ with the associated NHet and number of kills
                        outfile.open("hex_lognormal_radius_threshold_perfusion_quotients.txt", std::ios_base::app);
                        outfile << network_name << " " << solver_name << " " << heterogeneity_string << " " << selection_string << " " << threshold_string << " " << PerfusionQuotient << " \n"; 
                        outfile.close();

                        // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                        ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                        ParameterCollection::Instance()->Destroy();
                        BaseUnits::Instance()->Destroy();
                        SimulationTime::Instance()->Destroy();
                    }
                
                    // std::cout << "Attempts " << att+1 << " \n";
                    // outfileMean << "Attempts " << att+1 << " \n";
                    // Write the resulting perfusion quotients into an output file
                    // if (alpha==max_alpha)
                    // {   
                    //     outfile.open("hexagonal_perfusion_quotients.txt");
                    //     for (unsigned iii = 0; iii < max_alpha; iii++)
                    //     {
                    //         for(unsigned jjj = 0; jjj < ToBeKilled+1; jjj++)
                    //         {
                    //             QuotientMatrixMean[iii][jjj]= QuotientMatrixAggregate[iii][jjj]/double(list_number+1);
                    //             //std::cout << "For heterogeneity level " << iii << ", and radiotherapy dose " << jjj << ", the aggregate quotient is: " << QuotientMatrixAggregate[iii][jjj] << " , and the mean quotient is: " << QuotientMatrixMean[iii][jjj] << endl;
                    //             outfile << network_name << " " << solver_name << " " << iii << " " << jjj <<  " " << QuotientMatrixMean[iii][jjj] << " \n"; 
                    //             // outfile.close();
                    //             // outfileMean << iii << " " << jjj << " " << QuotientMatrixMean[iii][jjj] << " \n";
                    //         }
                    //     }    
                    //     outfile.close();
                    // }
                }
            }
        }
         
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a full 2D hexagonal network on a PDE grid with flow, H-splitting, radii set according to a log normal distribution, and kills-based pruning (correcting for the unexplained -50 um offset in the other simulation) with the Constant haematocrit solver
    void xTestConstantOffsetHeterogeneousHexagonalNeighbourhoodWithIndividualPruning2D()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";
        
        // Set network in filename
        std::string network_name = "TestHexagonalNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("/tmp/narain/testoutput/TestHexagonalNetwork/hex_lognormal_individual_perfusion_quotients.txt");
        outfile.close();

        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Define the key parameters
        // double dimless_domain_size_x = 2050.0;  // x-coordinate of output node
        // double dimless_domain_size_y = 1818.65 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
        // unsigned dimless_vessel_length = 100.0;
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        // double initial_haematocrit = 0.45;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        unsigned max_n_alpha = 3;
        unsigned n_vessels = 386;  // number of non-inlet/outlet vessels from which to select ones to make thin
        // double percToKill = 0.2;  // percentage of vessels to kill
        // unsigned ToBeKilled = (unsigned)(percToKill*n_vessels);  // number to kill
        unsigned ToBeKilled = 200;  // number to kill

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=1; h_solver++)
        {
            // Initialise the simulation
            std::shared_ptr<VesselNetwork<2> > p_network;
            std::vector<std::shared_ptr<Vessel<2> > > vessels;
            VesselPtr<2> p_current_vessel;
            QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
            QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
            std::vector<std::vector<double> > QuotientMatrixMean(max_n_alpha,std::vector<double>(ToBeKilled+1));
            std::vector<std::vector<double> > QuotientMatrixAggregate(max_n_alpha,std::vector<double>(ToBeKilled+1,0.0));
            std::vector<std::vector<unsigned> > Order;
            double PerfusionQuotient;

            // Read the network layout from a file and store it in an array
            VesselNetworkGenerator<2> network_generator;
            std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix_Offset.txt");
            std::vector<std::vector<double> > rEdgesMatrix;
            string line;
            while (std::getline(in, line)) 
            {
                rEdgesMatrix.push_back(std::vector<double>());
                std::stringstream split(line);
                double value;
                while (split >> value)
                {
                    rEdgesMatrix.back().push_back(value);
                }
            }

            // Loop over increasing level of heteregeneity (% of thin vessels)
            for(unsigned n_alpha=0; n_alpha<=max_n_alpha; n_alpha++)
            {
                // Convert to actual mu
                int alpha = (3*n_alpha)+10;  // alpha determines the mean

                // Loop over different thin vessel layouts
                for(unsigned list_number=1; list_number<=thin_selections; list_number++)
                {                
                    // Choose a radii list file and read the file to an array
                    std::vector<std::vector<double> > radii_array;
                    radii_array.clear();
                    // int list_index=0;
                    string line_1;
                    std::ifstream radius_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/hexagonal_radius_log_normal_distribution/mu_"+to_string(alpha)+"/radii_list_"+to_string(list_number)+".txt");
                    while (std::getline(radius_list_file, line_1)) 
                    {
                        radii_array.push_back(std::vector<double>());
                        std::stringstream split(line_1);
                        double value;
                        while (split >> value)
                        {
                            radii_array.back().push_back(value);
                        }
                    }

                    // Choose the corresponding vessel ID order file and read the file to an array
                    std::vector<std::vector<double> > id_array;
                    id_array.clear();
                    // int list_index=0;
                    string line_2;
                    std::ifstream id_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/hexagonal_radius_log_normal_distribution/mu_"+to_string(alpha)+"/id_list_"+to_string(list_number)+".txt");
                    while (std::getline(id_list_file, line_2)) 
                    {
                        id_array.push_back(std::vector<double>());
                        std::stringstream split(line_2);
                        double value;
                        while (split >> value)
                        {
                            id_array.back().push_back(value);
                        }
                    }

                    // Generate the network
                    p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);

                    // Set inlet and outlet nodes (if node (start/end) is at the edge of the domain space, make it input or output depending on which side), and assign each non-inlet/outlet vessel a radius from the list
                    auto p_segment = p_network->GetVesselSegments()[0];
                    p_segment->SetRadius(inlet_vessel_radius);
                    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   
                    vessels = p_network->GetVessels();
                    for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
                    {
                        if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the start node equals one (no other segments at nodes)
                        {
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the start node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);  // set the node as the input node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the start node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the ouput node pressure
                            }
                        }
                        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the end node equals one (no other segments at nodes)	
                        {
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the end node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);;  // set the node as the input node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the end node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the output node pressure
                            }
                        }
                    }
                
                    // Assign vessel radii based on imported arrays
                    unsigned int vessel_index = 0;
                    while(vessel_index < n_vessels)
                    {
                        p_current_vessel = p_network->GetVessel(id_array[vessel_index][0]);

                        // Get segment radius from array
                        QLength new_vessel_radius(double(radii_array[vessel_index][0])*unit::microns);

                        p_current_vessel->SetRadius(new_vessel_radius);
                        vessel_index++;
                    }
                    vessels = p_network->GetVessels();

                    // Set the haematocrit solver
                    std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                    std::string solver_name;
                    if (h_solver==1)
                    {
                        solver_name = "ConstantHaematocrit"; 
                        std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==2)
                    {
                        solver_name = "PriesHaematocrit"; 
                        std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                
                    }
                    else if (h_solver==3)
                    {
                        solver_name = "MemoryHaematocrit"; 
                        std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==4)
                    {
                        solver_name = "FungHaematocrit"; 
                        std::cout << "Now using FungHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==5)
                    {
                        solver_name = "YangHaematocrit"; 
                        std::cout << "Now using YangHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();   
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }

                    // Set up the viscosity solver
                    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                    p_viscosity_calculator->SetVesselNetwork(p_network);
                    p_viscosity_calculator->Calculate();
                    
                    // Set up the impedance solver
                    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                    p_impedance_calculator->SetVesselNetwork(p_network);
                    p_impedance_calculator->Calculate();

                    // Set up the flow solver 
                    FlowSolver<2> flow_solver;
                    flow_solver.SetVesselNetwork(p_network);
                    // flow_solver.SetUp();
                    flow_solver.SetUseDirectSolver(true);
                    // flow_solver.Solve();

                    // Set up the grid for the finite difference solver
                    // auto p_grid = RegularGrid<2>::Create();
                    // p_grid->SetSpacing(grid_spacing);
                    // c_vector<unsigned, 3> dimensions;
                    // dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                    // dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                    // dimensions[2] = 1;
                    // p_grid->SetDimensions(dimensions);

                    // Choose the PDE
                    // auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                    
                    // Set the diffusivity and decay terms
                    // p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                    // p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                    // Set up the discrete source
                    // auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                    // QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                    //         GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                    // QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                    //         Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                    // p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                    // p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                    // p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    // p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                    // // Set up the finite difference solver for oxygen (which handles everything)
                    // auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                    // p_oxygen_solver->SetPde(p_oxygen_pde);
                    // p_oxygen_solver->SetLabel("oxygen");
                    // p_oxygen_solver->SetGrid(p_grid);
                    // solver.SetFileName("oxygen_solution_0");

                    // Prune all vessels up to specified dose 
                    for(unsigned KilledVessels=0; KilledVessels < ToBeKilled+1; KilledVessels++)
                    { 
                        // Display status message
                        std::cout << "Now killed " << KilledVessels << " vessels." << std::endl;  

                        // Set filename
                        std::stringstream selection_stream;
                        std::stringstream heterogeneity_stream;
                        std::stringstream kill_stream;
                        selection_stream << std::fixed << std::setprecision(0) << to_string(list_number);                    
                        heterogeneity_stream << std::fixed << std::setprecision(0) << alpha;                              
                        kill_stream << std::fixed << std::setprecision(0) << KilledVessels;           
                        std::string selection_string = selection_stream.str();
                        std::string heterogeneity_string = heterogeneity_stream.str();                    
                        std::string kill_string = kill_stream.str();                    
                        std::string file_name = network_name + solver_name + "/Mu" + heterogeneity_string + "/Selection" + selection_string + "/Kills" + kill_string;
                        std::string str_directory_name = file_name;
                        std::cout << str_directory_name << std::endl;
                        auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                        // Set up flag for broken solver
                        unsigned broken_solver = 0;

                        // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                        // unsigned max_iter = 5;  // 1000 
                        // double tolerance2 = 1.e-10;
                        std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                        // std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                        // for(unsigned idx=0; idx<max_iter; idx++)
                        // {
                            // Run the solvers
                            p_impedance_calculator->Calculate();
                            p_abstract_haematocrit_solver->Calculate();
                            p_viscosity_calculator->Calculate();
                            flow_solver.SetUp();
                            flow_solver.Solve();

                            // Get the residual
                            // double max_difference = 0.0;
                            // double h_for_max = 0.0;
                            // double prev_for_max = 0.0;
                            for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                            {
                                // Set segments with no flow to have no haematocrit
                                if (fabs(segments[jdx]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                {
                                    segments[jdx]->GetFlowProperties()->SetHaematocrit(0.0);
                                }   
                                // double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                                // double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                                // if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                                // {
                                //     max_difference = difference;
                                //     h_for_max = current_haematocrit;
                                //     prev_for_max = previous_haematocrit[jdx];
                                // }
                                // previous_haematocrit[jdx] = current_haematocrit;
                            }
                        //     std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                        //     // Print the final or intermediary convergence results
                        //     if(max_difference<=tolerance2)  
                        //     {
                        //         std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                        //         broken_solver = 0;
                        //         break;
                        //     }
                        //     else
                        //     {
                        //         if(idx%1==0)
                        //         {
                        //             std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                        //             std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                        //             std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                        //             p_network->Write(output_file);
                        //         }
                        //     }

                        //     // If there is no convergence after all the iterations, print the error message.
                        //     if(idx==max_iter-1)
                        //     {
                        //         std::cout << "Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha << std::endl;
                        //         error_log << "\n Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha; 
                        //         broken_solver = 1;
                        //         break;
                        //         // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                        //     }
                        // }

                        // If solver doesn't converge, move on to next one
                        if (broken_solver == 1)
                        {
                            continue;
                        }

                        // Run the simulation 
                        // SimulationTime::Instance()->SetStartTime(0.0);
                        // SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                        // auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                        // p_microvessel_solver->SetVesselNetwork(p_network);
                        // p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                        // p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                        // p_microvessel_solver->Run();
                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                        p_network->Write(output_file);

                        // Print the average oxygenation
                        // std::vector<double> solution = p_oxygen_solver->GetSolution();
                        // double average_oxygen = 0.0;
                        // for(unsigned jdx=0;jdx<solution.size();jdx++)
                        // {
                        //     average_oxygen += solution[jdx];
                        // }
                        // average_oxygen /= double(solution.size());
                        // std::cout << "Average oxygen: " << average_oxygen << std::endl;

                        // // Write the PQs
                        // QFlowRate threshold = 5.e-14*unit::metre_cubed_per_second;
                        QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;
                        // QFlowRate threshold = 2.e-13*unit::metre_cubed_per_second;
                        // // outfile.open("hexagonal_perfusion_quotients.txt", std::ios_base::app);
                        // // outfile << network_name << " " << solver_name << " " << selection_string << " " << heterogeneity_string << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                        // // outfile.close();  
                        PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
                        // QuotientMatrixAggregate[i][KilledVessels] += PerfusionQuotient;  // store the PQ with the associated NHet and number of kills
                        outfile.open("/tmp/narain/testoutput/TestHexagonalNetwork/hex_lognormal_individual_perfusion_quotients.txt", std::ios_base::app);
                        outfile << network_name << " " << solver_name << " " << heterogeneity_string << " " << selection_string << " " << kill_string << " " << PerfusionQuotient << " \n"; 
                        outfile.close();

                        // Remove the smallest vessel
                        QLength minimum_radius = inlet_vessel_radius;
                        unsigned int minimum_index = 0;
                        for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)  // for all the segments in the network
                        {
                            // Get the current segment's radius
                            QLength current_radius = vessels[vessel_index]->GetRadius();

                            // If the current radius is less than the minimum radius, record the new minimum
                            if (current_radius < minimum_radius)
                            {
                                minimum_radius = current_radius;
                                minimum_index = vessel_index;
                            }
                        }
                        p_network->RemoveVessel(vessels[minimum_index], true);  // remove the vessel

                        // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                        ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                        ParameterCollection::Instance()->Destroy();
                        BaseUnits::Instance()->Destroy();
                        SimulationTime::Instance()->Destroy();
                    }
                }
            }
        }
         
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a full 2D hexagonal network on a PDE grid with flow, H-splitting, radii set according to a log normal distribution, and kills-based pruning in order of flow rate (correcting for the unexplained -50 um offset in the other simulation)
    void xTestConstantOffsetHeterogeneousHexagonalNeighbourhoodWithIndividualFlowRatePruning2D()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";
        
        // Set network in filename
        std::string network_name = "TestHexagonalNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("hex_lognormal_individual_perfusion_quotients.txt");
        outfile.close();

        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Define the key parameters
        // double dimless_domain_size_x = 2050.0;  // x-coordinate of output node
        // double dimless_domain_size_y = 1818.65 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
        // unsigned dimless_vessel_length = 100.0;
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        // double initial_haematocrit = 0.45;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        unsigned max_n_alpha = 3;
        unsigned n_vessels = 386;  // number of non-inlet/outlet vessels from which to select ones to make thin
        // double percToKill = 0.2;  // percentage of vessels to kill
        // unsigned ToBeKilled = (unsigned)(percToKill*n_vessels);  // number to kill
        unsigned ToBeKilled = 200;  // number to kill

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=1; h_solver++)
        {
            // Initialise the simulation
            std::shared_ptr<VesselNetwork<2> > p_network;
            std::vector<std::shared_ptr<Vessel<2> > > vessels;
            VesselPtr<2> p_current_vessel;
            QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
            QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
            std::vector<std::vector<double> > QuotientMatrixMean(max_n_alpha,std::vector<double>(ToBeKilled+1));
            std::vector<std::vector<double> > QuotientMatrixAggregate(max_n_alpha,std::vector<double>(ToBeKilled+1,0.0));
            std::vector<std::vector<unsigned> > Order;
            double PerfusionQuotient;

            // Read the network layout from a file and store it in an array
            VesselNetworkGenerator<2> network_generator;
            std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix_Offset.txt");
            std::vector<std::vector<double> > rEdgesMatrix;
            string line;
            while (std::getline(in, line)) 
            {
                rEdgesMatrix.push_back(std::vector<double>());
                std::stringstream split(line);
                double value;
                while (split >> value)
                {
                    rEdgesMatrix.back().push_back(value);
                }
            }

            // Loop over increasing level of heteregeneity (% of thin vessels)
            for(unsigned n_alpha=0; n_alpha<=max_n_alpha; n_alpha++)
            {
                // Convert to actual mu
                int alpha = (3*n_alpha)+10;  // alpha determines the mean

                // Loop over different thin vessel layouts
                for(unsigned list_number=1; list_number<=thin_selections; list_number++)
                {                
                    // Choose a radii list file and read the file to an array
                    std::vector<std::vector<double> > radii_array;
                    radii_array.clear();
                    // int list_index=0;
                    string line_1;
                    std::ifstream radius_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/hexagonal_radius_log_normal_distribution/mu_"+to_string(alpha)+"/radii_list_"+to_string(list_number)+".txt");
                    while (std::getline(radius_list_file, line_1)) 
                    {
                        radii_array.push_back(std::vector<double>());
                        std::stringstream split(line_1);
                        double value;
                        while (split >> value)
                        {
                            radii_array.back().push_back(value);
                        }
                    }

                    // Choose the corresponding vessel ID order file and read the file to an array
                    std::vector<std::vector<double> > id_array;
                    id_array.clear();
                    // int list_index=0;
                    string line_2;
                    std::ifstream id_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/hexagonal_radius_log_normal_distribution/mu_"+to_string(alpha)+"/id_list_"+to_string(list_number)+".txt");
                    while (std::getline(id_list_file, line_2)) 
                    {
                        id_array.push_back(std::vector<double>());
                        std::stringstream split(line_2);
                        double value;
                        while (split >> value)
                        {
                            id_array.back().push_back(value);
                        }
                    }

                    // Generate the network
                    p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);

                    // Set inlet and outlet nodes (if node (start/end) is at the edge of the domain space, make it input or output depending on which side), and assign each non-inlet/outlet vessel a radius from the list
                    auto p_segment = p_network->GetVesselSegments()[0];
                    p_segment->SetRadius(inlet_vessel_radius);
                    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   
                    vessels = p_network->GetVessels();
                    for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
                    {
                        if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the start node equals one (no other segments at nodes)
                        {
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the start node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);  // set the node as the input node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the start node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the ouput node pressure
                            }
                        }
                        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the end node equals one (no other segments at nodes)	
                        {
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the end node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);;  // set the node as the input node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the end node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the output node pressure
                            }
                        }
                    }
                
                    // Assign vessel radii based on imported arrays
                    unsigned int vessel_index = 0;
                    while(vessel_index < n_vessels)
                    {
                        p_current_vessel = p_network->GetVessel(id_array[vessel_index][0]);

                        // Get segment radius from array
                        QLength new_vessel_radius(double(radii_array[vessel_index][0])*unit::microns);

                        p_current_vessel->SetRadius(new_vessel_radius);
                        vessel_index++;
                    }
                    vessels = p_network->GetVessels();

                    // Set the haematocrit solver
                    std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                    std::string solver_name;
                    if (h_solver==1)
                    {
                        solver_name = "ConstantHaematocrit"; 
                        std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==2)
                    {
                        solver_name = "PriesHaematocrit"; 
                        std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                
                    }
                    else if (h_solver==3)
                    {
                        solver_name = "MemoryHaematocrit"; 
                        std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==4)
                    {
                        solver_name = "FungHaematocrit"; 
                        std::cout << "Now using FungHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==5)
                    {
                        solver_name = "YangHaematocrit"; 
                        std::cout << "Now using YangHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();   
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }

                    // Set up the viscosity solver
                    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                    p_viscosity_calculator->SetVesselNetwork(p_network);
                    p_viscosity_calculator->Calculate();
                    
                    // Set up the impedance solver
                    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                    p_impedance_calculator->SetVesselNetwork(p_network);
                    p_impedance_calculator->Calculate();

                    // Set up the flow solver 
                    FlowSolver<2> flow_solver;
                    flow_solver.SetVesselNetwork(p_network);
                    // flow_solver.SetUp();
                    flow_solver.SetUseDirectSolver(true);
                    // flow_solver.Solve();

                    // Set up the grid for the finite difference solver
                    // auto p_grid = RegularGrid<2>::Create();
                    // p_grid->SetSpacing(grid_spacing);
                    // c_vector<unsigned, 3> dimensions;
                    // dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                    // dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                    // dimensions[2] = 1;
                    // p_grid->SetDimensions(dimensions);

                    // Choose the PDE
                    // auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                    
                    // Set the diffusivity and decay terms
                    // p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                    // p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                    // Set up the discrete source
                    // auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                    // QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                    //         GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                    // QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                    //         Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                    // p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                    // p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                    // p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    // p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                    // // Set up the finite difference solver for oxygen (which handles everything)
                    // auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                    // p_oxygen_solver->SetPde(p_oxygen_pde);
                    // p_oxygen_solver->SetLabel("oxygen");
                    // p_oxygen_solver->SetGrid(p_grid);
                    // solver.SetFileName("oxygen_solution_0");

                    // Prune all vessels up to specified dose 
                    for(unsigned KilledVessels=0; KilledVessels < ToBeKilled+1; KilledVessels++)
                    { 
                        // Display status message
                        std::cout << "Now killed " << KilledVessels << " vessels." << std::endl;  

                        // Set filename
                        std::stringstream selection_stream;
                        std::stringstream heterogeneity_stream;
                        std::stringstream kill_stream;
                        selection_stream << std::fixed << std::setprecision(0) << to_string(list_number);                    
                        heterogeneity_stream << std::fixed << std::setprecision(0) << alpha;                              
                        kill_stream << std::fixed << std::setprecision(0) << KilledVessels;           
                        std::string selection_string = selection_stream.str();
                        std::string heterogeneity_string = heterogeneity_stream.str();                    
                        std::string kill_string = kill_stream.str();                    
                        std::string file_name = network_name + solver_name + "/Mu" + heterogeneity_string + "/Selection" + selection_string + "/Kills" + kill_string;
                        std::string str_directory_name = file_name;
                        std::cout << str_directory_name << std::endl;
                        auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                        // Set up flag for broken solver
                        // unsigned broken_solver = 0;

                        // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                        // unsigned max_iter = 5;  // 1000 
                        // double tolerance2 = 1.e-10;
                        std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                        // std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                        // for(unsigned idx=0; idx<max_iter; idx++)
                        // {
                            // Run the solvers
                            p_impedance_calculator->Calculate();
                            p_abstract_haematocrit_solver->Calculate();
                            p_viscosity_calculator->Calculate();
                            flow_solver.SetUp();
                            flow_solver.Solve();

                            // Get the residual
                            // double max_difference = 0.0;
                            // double h_for_max = 0.0;
                            // double prev_for_max = 0.0;
                            for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                            {
                                // Set segments with no flow to have no haematocrit
                                if (fabs(segments[jdx]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                {
                                    segments[jdx]->GetFlowProperties()->SetHaematocrit(0.0);
                                }   
                                // double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                                // double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                                // if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                                // {
                                //     max_difference = difference;
                                //     h_for_max = current_haematocrit;
                                //     prev_for_max = previous_haematocrit[jdx];
                                // }
                                // previous_haematocrit[jdx] = current_haematocrit;
                            }
                        //     std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                        //     // Print the final or intermediary convergence results
                        //     if(max_difference<=tolerance2)  
                        //     {
                        //         std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                        //         broken_solver = 0;
                        //         break;
                        //     }
                        //     else
                        //     {
                        //         if(idx%1==0)
                        //         {
                        //             std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                        //             std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                        //             std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                        //             p_network->Write(output_file);
                        //         }
                        //     }

                        //     // If there is no convergence after all the iterations, print the error message.
                        //     if(idx==max_iter-1)
                        //     {
                        //         std::cout << "Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha << std::endl;
                        //         error_log << "\n Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha; 
                        //         broken_solver = 1;
                        //         break;
                        //         // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                        //     }
                        // }

                        // // If solver doesn't converge, move on to next one
                        // if (broken_solver == 1)
                        // {
                        //     continue;
                        // }

                        // Run the simulation 
                        // SimulationTime::Instance()->SetStartTime(0.0);
                        // SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                        // auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                        // p_microvessel_solver->SetVesselNetwork(p_network);
                        // p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                        // p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                        // p_microvessel_solver->Run();
                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                        p_network->Write(output_file);

                        // Print the average oxygenation
                        // std::vector<double> solution = p_oxygen_solver->GetSolution();
                        // double average_oxygen = 0.0;
                        // for(unsigned jdx=0;jdx<solution.size();jdx++)
                        // {
                        //     average_oxygen += solution[jdx];
                        // }
                        // average_oxygen /= double(solution.size());
                        // std::cout << "Average oxygen: " << average_oxygen << std::endl;

                        // // Write the PQs
                        // QFlowRate threshold = 5.e-14*unit::metre_cubed_per_second;
                        QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;
                        // QFlowRate threshold = 2.e-13*unit::metre_cubed_per_second;
                        // // outfile.open("hexagonal_perfusion_quotients.txt", std::ios_base::app);
                        // // outfile << network_name << " " << solver_name << " " << selection_string << " " << heterogeneity_string << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                        // // outfile.close();  
                        PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
                        // QuotientMatrixAggregate[i][KilledVessels] += PerfusionQuotient;  // store the PQ with the associated NHet and number of kills
                        outfile.open("hex_lognormal_individual_perfusion_quotients.txt", std::ios_base::app);
                        outfile << network_name << " " << solver_name << " " << heterogeneity_string << " " << selection_string << " " << kill_string << " " << PerfusionQuotient << " \n"; 
                        outfile.close();
                        
                        // Remove no-flow vessels
                        vessels = p_network->GetVessels();
                        for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)  // for all the segments in the network
                        {
                            // Get the current segment's flow rate
                            QFlowRate current_flow = vessels[vessel_index]->GetFlowProperties()->GetFlowRate();

                            // If the current flow rate is zero, remove the vessel
                            if (current_flow == 0)
                            {
                                p_network->RemoveVessel(vessels[vessel_index], true);  // remove the vessel
                            }
                        }
                        
                        // Remove the vessels without the lowest flows
                        vessels = p_network->GetVessels();
                        if (vessels.size() == 0)  // if there are no vessels left, stop pruning
                        {
                            break;
                        }
                        QFlowRate minimum_flow = 1.e-0*unit::metre_cubed_per_second;
                        unsigned int minimum_index = 0;
                        for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)  // for all the segments in the network
                        {
                            // Get the current segment's flow rate
                            QFlowRate current_flow = vessels[vessel_index]->GetFlowProperties()->GetFlowRate();

                            // If the current flow rate is non-zero and less than the minimum flow rate, record the new minimum
                            if (current_flow < minimum_flow && current_flow > 0)
                            {
                                minimum_flow = current_flow;
                                minimum_index = vessel_index;
                            }
                        }
                        p_network->RemoveVessel(vessels[minimum_index], true);  // remove the vessel with the lowest flow

                        // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                        ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                        ParameterCollection::Instance()->Destroy();
                        BaseUnits::Instance()->Destroy();
                        SimulationTime::Instance()->Destroy();
                    }
                }
            }
        }
         
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a full single-feed 2D hexagonal network on a PDE grid with flow, H-splitting, radii set according to a log normal distribution, and kills-based pruning (correcting for the unexplained -50 um offset in the other simulation)
    void xTestConstantOffsetSingleFeedHeterogeneousHexagonalNeighbourhoodWithIndividualPruning2D()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";
        
        // Set network in filename
        std::string network_name = "TestHexagonalNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("hex_lognormal_individual_perfusion_quotients.txt");
        outfile.close();

        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Define the key parameters
        // double dimless_domain_size_x = 2050.0;  // x-coordinate of output node
        // double dimless_domain_size_y = 1818.65 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
        // unsigned dimless_vessel_length = 100.0;
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        // double initial_haematocrit = 0.45;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        unsigned max_n_alpha = 3;
        unsigned n_vessels = 386;  // number of non-inlet/outlet vessels from which to select ones to make thin
        // double percToKill = 0.2;  // percentage of vessels to kill
        // unsigned ToBeKilled = (unsigned)(percToKill*n_vessels);  // number to kill
        unsigned ToBeKilled = 200;  // number to kill

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=1; h_solver++)
        {
            // Initialise the simulation
            std::shared_ptr<VesselNetwork<2> > p_network;
            std::vector<std::shared_ptr<Vessel<2> > > vessels;
            VesselPtr<2> p_current_vessel;
            QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
            QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
            std::vector<std::vector<double> > QuotientMatrixMean(max_n_alpha,std::vector<double>(ToBeKilled+1));
            std::vector<std::vector<double> > QuotientMatrixAggregate(max_n_alpha,std::vector<double>(ToBeKilled+1,0.0));
            std::vector<std::vector<unsigned> > Order;
            double PerfusionQuotient;

            // Read the network layout from a file and store it in an array
            VesselNetworkGenerator<2> network_generator;
            std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix_Offset.txt");
            std::vector<std::vector<double> > rEdgesMatrix;
            string line;
            while (std::getline(in, line)) 
            {
                rEdgesMatrix.push_back(std::vector<double>());
                std::stringstream split(line);
                double value;
                while (split >> value)
                {
                    rEdgesMatrix.back().push_back(value);
                }
            }

            // Loop over increasing level of heteregeneity (% of thin vessels)
            for(unsigned n_alpha=0; n_alpha<=max_n_alpha; n_alpha++)
            {
                // Convert to actual mu
                int alpha = (3*n_alpha)+10;  // alpha determines the mean

                // Loop over different thin vessel layouts
                for(unsigned list_number=1; list_number<=thin_selections; list_number++)
                {                
                    // Choose a radii list file and read the file to an array
                    std::vector<std::vector<double> > radii_array;
                    radii_array.clear();
                    // int list_index=0;
                    string line_1;
                    std::ifstream radius_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/hexagonal_radius_log_normal_distribution/mu_"+to_string(alpha)+"/radii_list_"+to_string(list_number)+".txt");
                    while (std::getline(radius_list_file, line_1)) 
                    {
                        radii_array.push_back(std::vector<double>());
                        std::stringstream split(line_1);
                        double value;
                        while (split >> value)
                        {
                            radii_array.back().push_back(value);
                        }
                    }

                    // Choose the corresponding vessel ID order file and read the file to an array
                    std::vector<std::vector<double> > id_array;
                    id_array.clear();
                    // int list_index=0;
                    string line_2;
                    std::ifstream id_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/hexagonal_radius_log_normal_distribution/mu_"+to_string(alpha)+"/id_list_"+to_string(list_number)+".txt");
                    while (std::getline(id_list_file, line_2)) 
                    {
                        id_array.push_back(std::vector<double>());
                        std::stringstream split(line_2);
                        double value;
                        while (split >> value)
                        {
                            id_array.back().push_back(value);
                        }
                    }

                    // Generate the network
                    p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);

                    // Set inlet and outlet nodes (if node (start/end) is at the edge of the domain space, make it input or output depending on which side), and assign each non-inlet/outlet vessel a radius from the list
                    auto p_segment = p_network->GetVesselSegments()[0];
                    p_segment->SetRadius(inlet_vessel_radius);
                    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   
                    vessels = p_network->GetVessels();
                    for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
                    {
                        if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the start node equals one (no other segments at nodes)
                        {
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the start node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);  // set the node as the input node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the start node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the ouput node pressure
                            }
                        }
                        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the end node equals one (no other segments at nodes)	
                        {
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the end node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);;  // set the node as the input node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the end node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the output node pressure
                            }
                        }
                    }
                
                    // Assign vessel radii based on imported arrays
                    unsigned int vessel_index = 0;
                    while(vessel_index < n_vessels)
                    {
                        p_current_vessel = p_network->GetVessel(id_array[vessel_index][0]);

                        // Get segment radius from array
                        QLength new_vessel_radius(double(radii_array[vessel_index][0])*unit::microns);

                        p_current_vessel->SetRadius(new_vessel_radius);
                        vessel_index++;
                    }
                    vessels = p_network->GetVessels();

                    // // Remove all but one inlet and one outlet
                    // unsigned int extra_vessel_index = 0;
                    // p_network->RemoveVessel(vessels[30], true);  // remove the extra connecting vessel
                    // while(extra_vessel_index < n_vessels)
                    // {
                    //     // Remove extra inlets
                    //     if (extra_vessel_index<11 && extra_vessel_index>0) 
                    //     {
                    //         // Remove the vessel
                    //         p_network->RemoveVessel(vessels[extra_vessel_index], true);  // remove the vessel
                    //     }
                    //     // Remove extra outlets
                    //     if (extra_vessel_index>396 && extra_vessel_index<406)
                    //     {
                    //         // Remove the vessel
                    //         p_network->RemoveVessel(vessels[extra_vessel_index], true);  // remove the vessel
                    //     } 
                    //     extra_vessel_index++;
                    // }                    
                    // vessels = p_network->GetVessels();


                    // Remove all but one inlet and one outlet
                    p_network->TransformHexagonalNetwork(false);  
                    vessels = p_network->GetVessels();

                    // Set the haematocrit solver
                    std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                    std::string solver_name;
                    if (h_solver==1)
                    {
                        solver_name = "ConstantHaematocrit"; 
                        std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==2)
                    {
                        solver_name = "PriesHaematocrit"; 
                        std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                
                    }
                    else if (h_solver==3)
                    {
                        solver_name = "MemoryHaematocrit"; 
                        std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==4)
                    {
                        solver_name = "FungHaematocrit"; 
                        std::cout << "Now using FungHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }

                    // Set up the viscosity solver
                    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                    p_viscosity_calculator->SetVesselNetwork(p_network);
                    p_viscosity_calculator->Calculate();
                    
                    // Set up the impedance solver
                    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                    p_impedance_calculator->SetVesselNetwork(p_network);
                    p_impedance_calculator->Calculate();

                    // Set up the flow solver 
                    FlowSolver<2> flow_solver;
                    flow_solver.SetVesselNetwork(p_network);
                    // flow_solver.SetUp();
                    flow_solver.SetUseDirectSolver(true);
                    // flow_solver.Solve();

                    // Set up the grid for the finite difference solver
                    // auto p_grid = RegularGrid<2>::Create();
                    // p_grid->SetSpacing(grid_spacing);
                    // c_vector<unsigned, 3> dimensions;
                    // dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                    // dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                    // dimensions[2] = 1;
                    // p_grid->SetDimensions(dimensions);

                    // Choose the PDE
                    // auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                    
                    // Set the diffusivity and decay terms
                    // p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                    // p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                    // Set up the discrete source
                    // auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                    // QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                    //         GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                    // QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                    //         Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                    // p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                    // p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                    // p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    // p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                    // // Set up the finite difference solver for oxygen (which handles everything)
                    // auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                    // p_oxygen_solver->SetPde(p_oxygen_pde);
                    // p_oxygen_solver->SetLabel("oxygen");
                    // p_oxygen_solver->SetGrid(p_grid);
                    // solver.SetFileName("oxygen_solution_0");

                    // Prune all vessels up to specified dose 
                    for(unsigned KilledVessels=0; KilledVessels < ToBeKilled+1; KilledVessels++)
                    { 
                        // Display status message
                        std::cout << "Now killed " << KilledVessels << " vessels." << std::endl;  

                        // Set filename
                        std::stringstream selection_stream;
                        std::stringstream heterogeneity_stream;
                        std::stringstream kill_stream;
                        selection_stream << std::fixed << std::setprecision(0) << to_string(list_number);                    
                        heterogeneity_stream << std::fixed << std::setprecision(0) << alpha;                              
                        kill_stream << std::fixed << std::setprecision(0) << KilledVessels;           
                        std::string selection_string = selection_stream.str();
                        std::string heterogeneity_string = heterogeneity_stream.str();                    
                        std::string kill_string = kill_stream.str();                    
                        std::string file_name = network_name + solver_name + "/Mu" + heterogeneity_string + "/Selection" + selection_string + "/Kills" + kill_string;
                        std::string str_directory_name = file_name;
                        std::cout << str_directory_name << std::endl;
                        auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                        // Set up flag for broken solver
                        unsigned broken_solver = 0;

                        // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                        // unsigned max_iter = 5;  // 1000 
                        // double tolerance2 = 1.e-10;
                        std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                        // std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                        // for(unsigned idx=0; idx<max_iter; idx++)
                        // {
                            // Run the solvers
                            p_impedance_calculator->Calculate();
                            p_abstract_haematocrit_solver->Calculate();
                            p_viscosity_calculator->Calculate();
                            flow_solver.SetUp();
                            flow_solver.Solve();

                            // Get the residual
                            // double max_difference = 0.0;
                            // double h_for_max = 0.0;
                            // double prev_for_max = 0.0;
                            for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                            {
                                // Set segments with no flow to have no haematocrit
                                if (fabs(segments[jdx]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                {
                                    segments[jdx]->GetFlowProperties()->SetHaematocrit(0.0);
                                }   
                                // double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                                // double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                                // if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                                // {
                                //     max_difference = difference;
                                //     h_for_max = current_haematocrit;
                                //     prev_for_max = previous_haematocrit[jdx];
                                // }
                                // previous_haematocrit[jdx] = current_haematocrit;
                            }
                        //     std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                        //     // Print the final or intermediary convergence results
                        //     if(max_difference<=tolerance2)  
                        //     {
                        //         std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                        //         broken_solver = 0;
                        //         break;
                        //     }
                        //     else
                        //     {
                        //         if(idx%1==0)
                        //         {
                        //             std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                        //             std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                        //             std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                        //             p_network->Write(output_file);
                        //         }
                        //     }

                        //     // If there is no convergence after all the iterations, print the error message.
                        //     if(idx==max_iter-1)
                        //     {
                        //         std::cout << "Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha << std::endl;
                        //         error_log << "\n Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha; 
                        //         broken_solver = 1;
                        //         break;
                        //         // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                        //     }
                        // }

                        // If solver doesn't converge, move on to next one
                        if (broken_solver == 1)
                        {
                            continue;
                        }

                        // Run the simulation 
                        // SimulationTime::Instance()->SetStartTime(0.0);
                        // SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                        // auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                        // p_microvessel_solver->SetVesselNetwork(p_network);
                        // p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                        // p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                        // p_microvessel_solver->Run();
                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                        p_network->Write(output_file);

                        // Print the average oxygenation
                        // std::vector<double> solution = p_oxygen_solver->GetSolution();
                        // double average_oxygen = 0.0;
                        // for(unsigned jdx=0;jdx<solution.size();jdx++)
                        // {
                        //     average_oxygen += solution[jdx];
                        // }
                        // average_oxygen /= double(solution.size());
                        // std::cout << "Average oxygen: " << average_oxygen << std::endl;

                        // // Write the PQs
                        // QFlowRate threshold = 5.e-14*unit::metre_cubed_per_second;
                        QFlowRate threshold = 1.e-14*unit::metre_cubed_per_second;
                        // QFlowRate threshold = 2.e-13*unit::metre_cubed_per_second;
                        // // outfile.open("hexagonal_perfusion_quotients.txt", std::ios_base::app);
                        // // outfile << network_name << " " << solver_name << " " << selection_string << " " << heterogeneity_string << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                        // // outfile.close();  
                        PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
                        // QuotientMatrixAggregate[i][KilledVessels] += PerfusionQuotient;  // store the PQ with the associated NHet and number of kills
                        outfile.open("hex_lognormal_individual_perfusion_quotients.txt", std::ios_base::app);
                        outfile << network_name << " " << solver_name << " " << heterogeneity_string << " " << selection_string << " " << kill_string << " " << PerfusionQuotient << " \n"; 
                        outfile.close();

                        // Remove the smallest vessel
                        QLength minimum_radius = inlet_vessel_radius;
                        unsigned int minimum_index = 0;
                        for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)  // for all the segments in the network
                        {
                            // Get the current segment's radius
                            QLength current_radius = vessels[vessel_index]->GetRadius();

                            // If the current radius is less than the minimum radius, record the new minimum
                            if (current_radius < minimum_radius)
                            {
                                minimum_radius = current_radius;
                                minimum_index = vessel_index;
                            }
                        }
                        p_network->RemoveVessel(vessels[minimum_index], true);  // remove the vessel

                        // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                        ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                        ParameterCollection::Instance()->Destroy();
                        BaseUnits::Instance()->Destroy();
                        SimulationTime::Instance()->Destroy();
                    }
                }
            }
        }
         
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a full tall 2D hexagonal network on a PDE grid with flow, H-splitting, radii set according to a log normal distribution, and kills-based pruning (correcting for the unexplained -50 um offset in the other simulation)
    void xTestConstantTallOffsetHeterogeneousHexagonalNeighbourhoodWithIndividualPruning2D()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";
        
        // Set network in filename
        std::string network_name = "TestHexagonalNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("hex_lognormal_individual_perfusion_quotients.txt");
        outfile.close();

        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Define the key parameters
        double dimless_domain_size_x = 1150.0;  // x-coordinate of output node
        double dimless_domain_size_y = 3550.7 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
        // unsigned dimless_vessel_length = 100.0;
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        // double initial_haematocrit = 0.45;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        unsigned max_n_alpha = 3;
        unsigned n_vessels = 403;  // number of non-inlet/outlet vessels from which to select ones to make thin
        // double percToKill = 0.2;  // percentage of vessels to kill
        // unsigned ToBeKilled = (unsigned)(percToKill*n_vessels);  // number to kill
        unsigned ToBeKilled = 200;  // number to kill

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=1; h_solver++)
        {
            // Initialise the simulation
            std::shared_ptr<VesselNetwork<2> > p_network;
            std::vector<std::shared_ptr<Vessel<2> > > vessels;
            VesselPtr<2> p_current_vessel;
            QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
            QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
            std::vector<std::vector<double> > QuotientMatrixMean(max_n_alpha,std::vector<double>(ToBeKilled+1));
            std::vector<std::vector<double> > QuotientMatrixAggregate(max_n_alpha,std::vector<double>(ToBeKilled+1,0.0));
            std::vector<std::vector<unsigned> > Order;
            double PerfusionQuotient;

            // Read the network layout from a file and store it in an array
            VesselNetworkGenerator<2> network_generator;
            std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/TallEdgesMatrix_Offset.txt");
            std::vector<std::vector<double> > rEdgesMatrix;
            string line;
            while (std::getline(in, line)) 
            {
                rEdgesMatrix.push_back(std::vector<double>());
                std::stringstream split(line);
                double value;
                while (split >> value)
                {
                    rEdgesMatrix.back().push_back(value);
                }
            }

            // Loop over increasing level of heteregeneity (% of thin vessels)
            for(unsigned n_alpha=0; n_alpha<=max_n_alpha; n_alpha++)
            {
                // Convert to actual mu
                int alpha = (3*n_alpha)+10;  // alpha determines the mean

                // Loop over different thin vessel layouts
                for(unsigned list_number=1; list_number<=thin_selections; list_number++)
                {                
                    // Choose a radii list file and read the file to an array
                    std::vector<std::vector<double> > radii_array;
                    radii_array.clear();
                    // int list_index=0;
                    string line_1;
                    std::ifstream radius_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/tall_hexagonal_radius_log_normal_distribution/mu_"+to_string(alpha)+"/radii_list_"+to_string(list_number)+".txt");
                    while (std::getline(radius_list_file, line_1)) 
                    {
                        radii_array.push_back(std::vector<double>());
                        std::stringstream split(line_1);
                        double value;
                        while (split >> value)
                        {
                            radii_array.back().push_back(value);
                        }
                    }

                    // Choose the corresponding vessel ID order file and read the file to an array
                    std::vector<std::vector<double> > id_array;
                    id_array.clear();
                    // int list_index=0;
                    string line_2;
                    std::ifstream id_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/tall_hexagonal_radius_log_normal_distribution/mu_"+to_string(alpha)+"/id_list_"+to_string(list_number)+".txt");
                    while (std::getline(id_list_file, line_2)) 
                    {
                        id_array.push_back(std::vector<double>());
                        std::stringstream split(line_2);
                        double value;
                        while (split >> value)
                        {
                            id_array.back().push_back(value);
                        }
                    }

                    // Generate the network
                    p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);

                    // Set inlet and outlet nodes (if node (start/end) is at the edge of the domain space, make it input or output depending on which side), and assign each non-inlet/outlet vessel a radius from the list
                    auto p_segment = p_network->GetVesselSegments()[0];
                    p_segment->SetRadius(inlet_vessel_radius);
                    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   
                    vessels = p_network->GetVessels();
                    for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
                    {
                        if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the start node equals one (no other segments at nodes)
                        {
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the start node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);  // set the node as the input node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the start node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the ouput node pressure
                            }
                        }
                        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the end node equals one (no other segments at nodes)	
                        {
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the end node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);;  // set the node as the input node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the end node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the output node pressure
                            }
                        }
                    }
                
                    // Assign vessel radii based on imported arrays
                    unsigned int vessel_index = 0;
                    while(vessel_index < n_vessels)
                    {
                        p_current_vessel = p_network->GetVessel(id_array[vessel_index][0]);

                        // Get segment radius from array
                        QLength new_vessel_radius(double(radii_array[vessel_index][0])*unit::microns);

                        p_current_vessel->SetRadius(new_vessel_radius);
                        vessel_index++;
                    }
                    vessels = p_network->GetVessels();

                    // Set the haematocrit solver
                    std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                    std::string solver_name;
                    if (h_solver==1)
                    {
                        solver_name = "ConstantHaematocrit"; 
                        std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==2)
                    {
                        solver_name = "PriesHaematocrit"; 
                        std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                
                    }
                    else if (h_solver==3)
                    {
                        solver_name = "MemoryHaematocrit"; 
                        std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==4)
                    {
                        solver_name = "FungHaematocrit"; 
                        std::cout << "Now using FungHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==5)
                    {
                        solver_name = "YangHaematocrit"; 
                        std::cout << "Now using YangHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();   
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }

                    // Set up the viscosity solver
                    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                    p_viscosity_calculator->SetVesselNetwork(p_network);
                    p_viscosity_calculator->Calculate();
                    
                    // Set up the impedance solver
                    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                    p_impedance_calculator->SetVesselNetwork(p_network);
                    p_impedance_calculator->Calculate();

                    // Set up the flow solver 
                    FlowSolver<2> flow_solver;
                    flow_solver.SetVesselNetwork(p_network);
                    // flow_solver.SetUp();
                    flow_solver.SetUseDirectSolver(true);
                    // flow_solver.Solve();

                    // Set up the grid for the finite difference solver
                    // auto p_grid = RegularGrid<2>::Create();
                    // p_grid->SetSpacing(grid_spacing);
                    // c_vector<unsigned, 3> dimensions;
                    // dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                    // dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                    // dimensions[2] = 1;
                    // p_grid->SetDimensions(dimensions);

                    // Choose the PDE
                    // auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                    
                    // Set the diffusivity and decay terms
                    // p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                    // p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                    // Set up the discrete source
                    // auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                    // QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                    //         GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                    // QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                    //         Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                    // p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                    // p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                    // p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    // p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                    // // Set up the finite difference solver for oxygen (which handles everything)
                    // auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                    // p_oxygen_solver->SetPde(p_oxygen_pde);
                    // p_oxygen_solver->SetLabel("oxygen");
                    // p_oxygen_solver->SetGrid(p_grid);
                    // solver.SetFileName("oxygen_solution_0");

                    // Prune all vessels up to specified dose 
                    for(unsigned KilledVessels=0; KilledVessels < ToBeKilled+1; KilledVessels++)
                    { 
                        // Display status message
                        std::cout << "Now killed " << KilledVessels << " vessels." << std::endl;  

                        // Set filename
                        std::stringstream selection_stream;
                        std::stringstream heterogeneity_stream;
                        std::stringstream kill_stream;
                        selection_stream << std::fixed << std::setprecision(0) << to_string(list_number);                    
                        heterogeneity_stream << std::fixed << std::setprecision(0) << alpha;                              
                        kill_stream << std::fixed << std::setprecision(0) << KilledVessels;           
                        std::string selection_string = selection_stream.str();
                        std::string heterogeneity_string = heterogeneity_stream.str();                    
                        std::string kill_string = kill_stream.str();                    
                        std::string file_name = network_name + solver_name + "/Mu" + heterogeneity_string + "/Selection" + selection_string + "/Kills" + kill_string;
                        std::string str_directory_name = file_name;
                        std::cout << str_directory_name << std::endl;
                        auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                        // Set up flag for broken solver
                        unsigned broken_solver = 0;

                        // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                        // unsigned max_iter = 5;  // 1000 
                        // double tolerance2 = 1.e-10;
                        std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                        // std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                        // for(unsigned idx=0; idx<max_iter; idx++)
                        // {
                            // Run the solvers
                            p_impedance_calculator->Calculate();
                            p_abstract_haematocrit_solver->Calculate();
                            p_viscosity_calculator->Calculate();
                            flow_solver.SetUp();
                            flow_solver.Solve();

                            // Get the residual
                            // double max_difference = 0.0;
                            // double h_for_max = 0.0;
                            // double prev_for_max = 0.0;
                            for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                            {
                                // Set segments with no flow to have no haematocrit
                                if (fabs(segments[jdx]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                {
                                    segments[jdx]->GetFlowProperties()->SetHaematocrit(0.0);
                                }   
                                // double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                                // double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                                // if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                                // {
                                //     max_difference = difference;
                                //     h_for_max = current_haematocrit;
                                //     prev_for_max = previous_haematocrit[jdx];
                                // }
                                // previous_haematocrit[jdx] = current_haematocrit;
                            }
                        //     std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                        //     // Print the final or intermediary convergence results
                        //     if(max_difference<=tolerance2)  
                        //     {
                        //         std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                        //         broken_solver = 0;
                        //         break;
                        //     }
                        //     else
                        //     {
                        //         if(idx%1==0)
                        //         {
                        //             std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                        //             std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                        //             std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                        //             p_network->Write(output_file);
                        //         }
                        //     }

                        //     // If there is no convergence after all the iterations, print the error message.
                        //     if(idx==max_iter-1)
                        //     {
                        //         std::cout << "Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha << std::endl;
                        //         error_log << "\n Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha; 
                        //         broken_solver = 1;
                        //         break;
                        //         // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                        //     }
                        // }

                        // If solver doesn't converge, move on to next one
                        if (broken_solver == 1)
                        {
                            continue;
                        }

                        // Run the simulation 
                        // SimulationTime::Instance()->SetStartTime(0.0);
                        // SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                        // auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                        // p_microvessel_solver->SetVesselNetwork(p_network);
                        // p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                        // p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                        // p_microvessel_solver->Run();
                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                        p_network->Write(output_file);

                        // Print the average oxygenation
                        // std::vector<double> solution = p_oxygen_solver->GetSolution();
                        // double average_oxygen = 0.0;
                        // for(unsigned jdx=0;jdx<solution.size();jdx++)
                        // {
                        //     average_oxygen += solution[jdx];
                        // }
                        // average_oxygen /= double(solution.size());
                        // std::cout << "Average oxygen: " << average_oxygen << std::endl;

                        // // Write the PQs
                        // QFlowRate threshold = 5.e-14*unit::metre_cubed_per_second;
                        QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;
                        // QFlowRate threshold = 2.e-13*unit::metre_cubed_per_second;
                        // // outfile.open("hexagonal_perfusion_quotients.txt", std::ios_base::app);
                        // // outfile << network_name << " " << solver_name << " " << selection_string << " " << heterogeneity_string << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                        // // outfile.close();  
                        PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
                        // QuotientMatrixAggregate[i][KilledVessels] += PerfusionQuotient;  // store the PQ with the associated NHet and number of kills
                        outfile.open("hex_lognormal_individual_perfusion_quotients.txt", std::ios_base::app);
                        outfile << network_name << " " << solver_name << " " << heterogeneity_string << " " << selection_string << " " << kill_string << " " << PerfusionQuotient << " \n"; 
                        outfile.close();

                        // Remove the smallest vessel
                        QLength minimum_radius = inlet_vessel_radius;
                        unsigned int minimum_index = 0;
                        for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)  // for all the segments in the network
                        {
                            // Get the current segment's radius
                            QLength current_radius = vessels[vessel_index]->GetRadius();

                            // If the current radius is less than the minimum radius, record the new minimum
                            if (current_radius < minimum_radius)
                            {
                                minimum_radius = current_radius;
                                minimum_index = vessel_index;
                            }
                        }
                        p_network->RemoveVessel(vessels[minimum_index], true);  // remove the vessel

                        // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                        ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                        ParameterCollection::Instance()->Destroy();
                        BaseUnits::Instance()->Destroy();
                        SimulationTime::Instance()->Destroy();
                    }
                }
            }
        }
         
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a full long 2D hexagonal network on a PDE grid with flow, H-splitting, radii set according to a log normal distribution, and kills-based pruning (correcting for the unexplained -50 um offset in the other simulation)
    void xTestConstantLongOffsetHeterogeneousHexagonalNeighbourhoodWithIndividualPruning2D()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";
        
        // Set network in filename
        std::string network_name = "TestHexagonalNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("hex_lognormal_individual_perfusion_quotients.txt");
        outfile.close();

        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Define the key parameters
        // double dimless_domain_size_x = 2050.0;  // x-coordinate of output node
        // double dimless_domain_size_y = 1818.65 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
        // unsigned dimless_vessel_length = 100.0;
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        // double initial_haematocrit = 0.45;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        unsigned max_n_alpha = 3;
        unsigned n_vessels = 382;  // number of non-inlet/outlet vessels from which to select ones to make thin
        // double percToKill = 0.2;  // percentage of vessels to kill
        // unsigned ToBeKilled = (unsigned)(percToKill*n_vessels);  // number to kill
        unsigned ToBeKilled = 200;  // number to kill

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=1; h_solver++)
        {
            // Initialise the simulation
            std::shared_ptr<VesselNetwork<2> > p_network;
            std::vector<std::shared_ptr<Vessel<2> > > vessels;
            VesselPtr<2> p_current_vessel;
            QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
            QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
            std::vector<std::vector<double> > QuotientMatrixMean(max_n_alpha,std::vector<double>(ToBeKilled+1));
            std::vector<std::vector<double> > QuotientMatrixAggregate(max_n_alpha,std::vector<double>(ToBeKilled+1,0.0));
            std::vector<std::vector<unsigned> > Order;
            double PerfusionQuotient;

            // Read the network layout from a file and store it in an array
            VesselNetworkGenerator<2> network_generator;
            std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/LongEdgesMatrix_Offset.txt");
            std::vector<std::vector<double> > rEdgesMatrix;
            string line;
            while (std::getline(in, line)) 
            {
                rEdgesMatrix.push_back(std::vector<double>());
                std::stringstream split(line);
                double value;
                while (split >> value)
                {
                    rEdgesMatrix.back().push_back(value);
                }
            }

            // Loop over increasing level of heteregeneity
            for(unsigned n_alpha=0; n_alpha<=max_n_alpha; n_alpha++)
            {
                // Convert to actual mu
                int alpha = (3*n_alpha)+10;  // alpha determines the mean

                // Loop over different thin vessel layouts
                for(unsigned list_number=1; list_number<=thin_selections; list_number++)
                {                
                    // Choose a radii list file and read the file to an array
                    std::vector<std::vector<double> > radii_array;
                    radii_array.clear();
                    // int list_index=0;
                    string line_1;
                    std::ifstream radius_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/long_hexagonal_radius_log_normal_distribution/mu_"+to_string(alpha)+"/radii_list_"+to_string(list_number)+".txt");
                    while (std::getline(radius_list_file, line_1)) 
                    {
                        radii_array.push_back(std::vector<double>());
                        std::stringstream split(line_1);
                        double value;
                        while (split >> value)
                        {
                            radii_array.back().push_back(value);
                        }
                    }

                    // Choose the corresponding vessel ID order file and read the file to an array
                    std::vector<std::vector<double> > id_array;
                    id_array.clear();
                    // int list_index=0;
                    string line_2;
                    std::ifstream id_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/long_hexagonal_radius_log_normal_distribution/mu_"+to_string(alpha)+"/id_list_"+to_string(list_number)+".txt");
                    while (std::getline(id_list_file, line_2)) 
                    {
                        id_array.push_back(std::vector<double>());
                        std::stringstream split(line_2);
                        double value;
                        while (split >> value)
                        {
                            id_array.back().push_back(value);
                        }
                    }

                    // Generate the network
                    p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);

                    // Set inlet and outlet nodes (if node (start/end) is at the edge of the domain space, make it input or output depending on which side), and assign each non-inlet/outlet vessel a radius from the list
                    auto p_segment = p_network->GetVesselSegments()[0];
                    p_segment->SetRadius(inlet_vessel_radius);
                    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   
                    vessels = p_network->GetVessels();
                    for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
                    {
                        if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the start node equals one (no other segments at nodes)
                        {
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the start node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);  // set the node as the input node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the start node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the ouput node pressure
                            }
                        }
                        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the end node equals one (no other segments at nodes)	
                        {
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the end node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);;  // set the node as the input node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the end node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the output node pressure
                            }
                        }
                    }
                
                    // Assign vessel radii based on imported arrays
                    unsigned int vessel_index = 0;
                    while(vessel_index < n_vessels)
                    {
                        p_current_vessel = p_network->GetVessel(id_array[vessel_index][0]);

                        // Get segment radius from array
                        QLength new_vessel_radius(double(radii_array[vessel_index][0])*unit::microns);

                        p_current_vessel->SetRadius(new_vessel_radius);

                        vessel_index++;
                    }
                    vessels = p_network->GetVessels();

                    // Set the haematocrit solver
                    std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                    std::string solver_name;
                    if (h_solver==1)
                    {
                        solver_name = "ConstantHaematocrit"; 
                        std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==2)
                    {
                        solver_name = "PriesHaematocrit"; 
                        std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                
                    }
                    else if (h_solver==3)
                    {
                        solver_name = "MemoryHaematocrit"; 
                        std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==4)
                    {
                        solver_name = "FungHaematocrit"; 
                        std::cout << "Now using FungHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==5)
                    {
                        solver_name = "YangHaematocrit"; 
                        std::cout << "Now using YangHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();   
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }

                    // Set up the viscosity solver
                    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                    p_viscosity_calculator->SetVesselNetwork(p_network);
                    p_viscosity_calculator->Calculate();
                    
                    // Set up the impedance solver
                    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                    p_impedance_calculator->SetVesselNetwork(p_network);
                    p_impedance_calculator->Calculate();

                    // Set up the flow solver 
                    FlowSolver<2> flow_solver;
                    flow_solver.SetVesselNetwork(p_network);
                    // flow_solver.SetUp();
                    flow_solver.SetUseDirectSolver(true);
                    // flow_solver.Solve();

                    // Set up the grid for the finite difference solver
                    // auto p_grid = RegularGrid<2>::Create();
                    // p_grid->SetSpacing(grid_spacing);
                    // c_vector<unsigned, 3> dimensions;
                    // dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                    // dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                    // dimensions[2] = 1;
                    // p_grid->SetDimensions(dimensions);

                    // Choose the PDE
                    // auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                    
                    // Set the diffusivity and decay terms
                    // p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                    // p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                    // Set up the discrete source
                    // auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                    // QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                    //         GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                    // QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                    //         Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                    // p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                    // p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                    // p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    // p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                    // // Set up the finite difference solver for oxygen (which handles everything)
                    // auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                    // p_oxygen_solver->SetPde(p_oxygen_pde);
                    // p_oxygen_solver->SetLabel("oxygen");
                    // p_oxygen_solver->SetGrid(p_grid);
                    // solver.SetFileName("oxygen_solution_0");

                    // Prune all vessels up to specified dose 
                    for(unsigned KilledVessels=0; KilledVessels < ToBeKilled+1; KilledVessels++)
                    { 
                        // Display status message
                        std::cout << "Now killed " << KilledVessels << " vessels." << std::endl;  

                        // Set filename
                        std::stringstream selection_stream;
                        std::stringstream heterogeneity_stream;
                        std::stringstream kill_stream;
                        selection_stream << std::fixed << std::setprecision(0) << to_string(list_number);                    
                        heterogeneity_stream << std::fixed << std::setprecision(0) << alpha;                              
                        kill_stream << std::fixed << std::setprecision(0) << KilledVessels;           
                        std::string selection_string = selection_stream.str();
                        std::string heterogeneity_string = heterogeneity_stream.str();                    
                        std::string kill_string = kill_stream.str();                    
                        std::string file_name = network_name + solver_name + "/Mu" + heterogeneity_string + "/Selection" + selection_string + "/Kills" + kill_string;
                        std::string str_directory_name = file_name;
                        std::cout << str_directory_name << std::endl;
                        auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                        // Set up flag for broken solver
                        unsigned broken_solver = 0;

                        // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                        // unsigned max_iter = 5;  // 1000 
                        // double tolerance2 = 1.e-10;
                        std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                        // std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                        // for(unsigned idx=0; idx<max_iter; idx++)
                        // {
                            // Run the solvers
                            p_impedance_calculator->Calculate();
                            p_abstract_haematocrit_solver->Calculate();
                            p_viscosity_calculator->Calculate();
                            flow_solver.SetUp();
                            flow_solver.Solve();

                            // Get the residual
                            // double max_difference = 0.0;
                            // double h_for_max = 0.0;
                            // double prev_for_max = 0.0;
                            for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                            {
                                // Set segments with no flow to have no haematocrit
                                if (fabs(segments[jdx]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                {
                                    segments[jdx]->GetFlowProperties()->SetHaematocrit(0.0);
                                }   
                                // double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                                // double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                                // if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                                // {
                                //     max_difference = difference;
                                //     h_for_max = current_haematocrit;
                                //     prev_for_max = previous_haematocrit[jdx];
                                // }
                                // previous_haematocrit[jdx] = current_haematocrit;
                            }
                        //     std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                        //     // Print the final or intermediary convergence results
                        //     if(max_difference<=tolerance2)  
                        //     {
                        //         std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                        //         broken_solver = 0;
                        //         break;
                        //     }
                        //     else
                        //     {
                        //         if(idx%1==0)
                        //         {
                        //             std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                        //             std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                        //             std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                        //             p_network->Write(output_file);
                        //         }
                        //     }

                        //     // If there is no convergence after all the iterations, print the error message.
                        //     if(idx==max_iter-1)
                        //     {
                        //         std::cout << "Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha << std::endl;
                        //         error_log << "\n Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha; 
                        //         broken_solver = 1;
                        //         break;
                        //         // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                        //     }
                        // }

                        // If solver doesn't converge, move on to next one
                        if (broken_solver == 1)
                        {
                            continue;
                        }

                        // Run the simulation 
                        // SimulationTime::Instance()->SetStartTime(0.0);
                        // SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                        // auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                        // p_microvessel_solver->SetVesselNetwork(p_network);
                        // p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                        // p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                        // p_microvessel_solver->Run();
                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                        p_network->Write(output_file);

                        // Print the average oxygenation
                        // std::vector<double> solution = p_oxygen_solver->GetSolution();
                        // double average_oxygen = 0.0;
                        // for(unsigned jdx=0;jdx<solution.size();jdx++)
                        // {
                        //     average_oxygen += solution[jdx];
                        // }
                        // average_oxygen /= double(solution.size());
                        // std::cout << "Average oxygen: " << average_oxygen << std::endl;

                        // // Write the PQs
                        // QFlowRate threshold = 5.e-14*unit::metre_cubed_per_second;
                        QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;
                        // QFlowRate threshold = 2.e-13*unit::metre_cubed_per_second;
                        // // outfile.open("hexagonal_perfusion_quotients.txt", std::ios_base::app);
                        // // outfile << network_name << " " << solver_name << " " << selection_string << " " << heterogeneity_string << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                        // // outfile.close();  
                        PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
                        // QuotientMatrixAggregate[i][KilledVessels] += PerfusionQuotient;  // store the PQ with the associated NHet and number of kills
                        outfile.open("hex_lognormal_individual_perfusion_quotients.txt", std::ios_base::app);
                        outfile << network_name << " " << solver_name << " " << heterogeneity_string << " " << selection_string << " " << kill_string << " " << PerfusionQuotient << " \n"; 
                        outfile.close();

                        // Remove the smallest vessel
                        QLength minimum_radius = inlet_vessel_radius;
                        unsigned int minimum_index = 0;
                        for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)  // for all the segments in the network
                        {
                            // Get the current segment's radius
                            QLength current_radius = vessels[vessel_index]->GetRadius();

                            // If the current radius is less than the minimum radius, record the new minimum
                            if (current_radius < minimum_radius)
                            {
                                minimum_radius = current_radius;
                                minimum_index = vessel_index;
                            }
                        }
                        p_network->RemoveVessel(vessels[minimum_index], true);  // remove the vessel

                        // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                        ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                        ParameterCollection::Instance()->Destroy();
                        BaseUnits::Instance()->Destroy();
                        SimulationTime::Instance()->Destroy();
                    }
                }
            }
        }
         
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a full 2D hexagonal network on a PDE grid with flow, H-splitting, all (not just non-inlet/outlet) radii set according to a log normal distribution, and kills-based pruning on all the vessels (correcting for the unexplained -50 um offset in the other simulation)
    void xTestConstantOffsetHeterogeneousHexagonalNeighbourhoodWithIndividualPruningForAllVesselsAndCompleteDiameterAssignment2D()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";
        
        // Set network in filename
        std::string network_name = "TestHexagonalNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("hex_lognormal_individual_perfusion_quotients.txt");
        outfile.close();

        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Define the key parameters
        // double dimless_domain_size_x = 2050.0;  // x-coordinate of output node
        // double dimless_domain_size_y = 1818.65 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
        // unsigned dimless_vessel_length = 100.0;
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        // double initial_haematocrit = 0.45;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        unsigned max_n_alpha = 3;
        unsigned n_vessels = 407;  // number of vessels from which to select ones to make thin
        // double percToKill = 0.2;  // percentage of vessels to kill
        // unsigned ToBeKilled = (unsigned)(percToKill*n_vessels);  // number to kill
        unsigned ToBeKilled = 200;  // number to kill

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=1; h_solver++)
        {
            // Initialise the simulation
            std::shared_ptr<VesselNetwork<2> > p_network;
            std::vector<std::shared_ptr<Vessel<2> > > vessels;
            VesselPtr<2> p_current_vessel;
            QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
            QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
            std::vector<std::vector<double> > QuotientMatrixMean(max_n_alpha,std::vector<double>(ToBeKilled+1));
            std::vector<std::vector<double> > QuotientMatrixAggregate(max_n_alpha,std::vector<double>(ToBeKilled+1,0.0));
            std::vector<std::vector<unsigned> > Order;
            double PerfusionQuotient;

            // Read the network layout from a file and store it in an array
            VesselNetworkGenerator<2> network_generator;
            std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix_Offset.txt");
            std::vector<std::vector<double> > rEdgesMatrix;
            string line;
            while (std::getline(in, line)) 
            {
                rEdgesMatrix.push_back(std::vector<double>());
                std::stringstream split(line);
                double value;
                while (split >> value)
                {
                    rEdgesMatrix.back().push_back(value);
                }
            }

            // Loop over increasing level of heteregeneity (% of thin vessels)
            for(unsigned n_alpha=0; n_alpha<=max_n_alpha; n_alpha++)
            {
                // Convert to actual mu
                int alpha = (3*n_alpha)+10;  // alpha determines the mean

                // Loop over different thin vessel layouts
                for(unsigned list_number=1; list_number<=thin_selections; list_number++)
                {                
                    // Choose a radii list file and read the file to an array
                    std::vector<std::vector<double> > radii_array;
                    radii_array.clear();
                    // int list_index=0;
                    string line_1;
                    std::ifstream radius_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/square_hexagonal_radius_log_normal_distribution/mu_"+to_string(alpha)+"/radii_list_"+to_string(list_number)+".txt");
                    while (std::getline(radius_list_file, line_1)) 
                    {
                        radii_array.push_back(std::vector<double>());
                        std::stringstream split(line_1);
                        double value;
                        while (split >> value)
                        {
                            radii_array.back().push_back(value);
                        }
                    }

                    // Choose the corresponding vessel ID order file and read the file to an array
                    std::vector<std::vector<double> > id_array;
                    id_array.clear();
                    // int list_index=0;
                    string line_2;
                    std::ifstream id_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/square_hexagonal_radius_log_normal_distribution/mu_"+to_string(alpha)+"/id_list_"+to_string(list_number)+".txt");
                    while (std::getline(id_list_file, line_2)) 
                    {
                        id_array.push_back(std::vector<double>());
                        std::stringstream split(line_2);
                        double value;
                        while (split >> value)
                        {
                            id_array.back().push_back(value);
                        }
                    }

                    // Generate the network
                    p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);

                    // Set inlet and outlet nodes (if node (start/end) is at the edge of the domain space, make it input or output depending on which side), and assign each non-inlet/outlet vessel a radius from the list
                    auto p_segment = p_network->GetVesselSegments()[0];
                    p_segment->SetRadius(inlet_vessel_radius);
                    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   
                    vessels = p_network->GetVessels();
                    for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
                    {
                        if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the start node equals one (no other segments at nodes)
                        {
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the start node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);  // set the node as the input node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the start node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the ouput node pressure
                            }
                        }
                        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the end node equals one (no other segments at nodes)	
                        {
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the end node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);;  // set the node as the input node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the end node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the output node pressure
                            }
                        }
                    }
                
                    // Assign vessel radii based on imported arrays
                    unsigned int vessel_index = 0;
                    while(vessel_index < n_vessels)
                    {
                        p_current_vessel = p_network->GetVessel(id_array[vessel_index][0]);

                        // Get segment radius from array
                        QLength new_vessel_radius(double(radii_array[vessel_index][0])*unit::microns);

                        p_current_vessel->SetRadius(new_vessel_radius);
                        vessel_index++;
                    }
                    vessels = p_network->GetVessels();

                    // Set the haematocrit solver
                    std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                    std::string solver_name;
                    if (h_solver==1)
                    {
                        solver_name = "ConstantHaematocrit"; 
                        std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==2)
                    {
                        solver_name = "PriesHaematocrit"; 
                        std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                
                    }
                    else if (h_solver==3)
                    {
                        solver_name = "MemoryHaematocrit"; 
                        std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==4)
                    {
                        solver_name = "FungHaematocrit"; 
                        std::cout << "Now using FungHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==5)
                    {
                        solver_name = "YangHaematocrit"; 
                        std::cout << "Now using YangHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();   
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }

                    // Set up the viscosity solver
                    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                    p_viscosity_calculator->SetVesselNetwork(p_network);
                    p_viscosity_calculator->Calculate();
                    
                    // Set up the impedance solver
                    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                    p_impedance_calculator->SetVesselNetwork(p_network);
                    p_impedance_calculator->Calculate();

                    // Set up the flow solver 
                    FlowSolver<2> flow_solver;
                    flow_solver.SetVesselNetwork(p_network);
                    // flow_solver.SetUp();
                    flow_solver.SetUseDirectSolver(true);
                    // flow_solver.Solve();

                    // Set up the grid for the finite difference solver
                    // auto p_grid = RegularGrid<2>::Create();
                    // p_grid->SetSpacing(grid_spacing);
                    // c_vector<unsigned, 3> dimensions;
                    // dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                    // dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                    // dimensions[2] = 1;
                    // p_grid->SetDimensions(dimensions);

                    // Choose the PDE
                    // auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                    
                    // Set the diffusivity and decay terms
                    // p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                    // p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                    // Set up the discrete source
                    // auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                    // QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                    //         GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                    // QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                    //         Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                    // p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                    // p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                    // p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    // p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                    // // Set up the finite difference solver for oxygen (which handles everything)
                    // auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                    // p_oxygen_solver->SetPde(p_oxygen_pde);
                    // p_oxygen_solver->SetLabel("oxygen");
                    // p_oxygen_solver->SetGrid(p_grid);
                    // solver.SetFileName("oxygen_solution_0");

                    // Prune all vessels up to specified dose 
                    for(unsigned KilledVessels=0; KilledVessels < ToBeKilled+1; KilledVessels++)
                    { 
                        // Display status message
                        std::cout << "Now killed " << KilledVessels << " vessels." << std::endl;  

                        // Set filename
                        std::stringstream selection_stream;
                        std::stringstream heterogeneity_stream;
                        std::stringstream kill_stream;
                        selection_stream << std::fixed << std::setprecision(0) << to_string(list_number);                    
                        heterogeneity_stream << std::fixed << std::setprecision(0) << alpha;                              
                        kill_stream << std::fixed << std::setprecision(0) << KilledVessels;           
                        std::string selection_string = selection_stream.str();
                        std::string heterogeneity_string = heterogeneity_stream.str();                    
                        std::string kill_string = kill_stream.str();                    
                        std::string file_name = network_name + solver_name + "/Mu" + heterogeneity_string + "/Selection" + selection_string + "/Kills" + kill_string;
                        std::string str_directory_name = file_name;
                        std::cout << str_directory_name << std::endl;
                        auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                        // Set up flag for broken solver
                        unsigned broken_solver = 0;

                        // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                        // unsigned max_iter = 5;  // 1000 
                        // double tolerance2 = 1.e-10;
                        std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                        // std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                        // for(unsigned idx=0; idx<max_iter; idx++)
                        // {
                            // Run the solvers
                            p_impedance_calculator->Calculate();
                            p_abstract_haematocrit_solver->Calculate();
                            p_viscosity_calculator->Calculate();
                            flow_solver.SetUp();
                            flow_solver.Solve();

                            // Get the residual
                            // double max_difference = 0.0;
                            // double h_for_max = 0.0;
                            // double prev_for_max = 0.0;
                            for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                            {
                                // Set segments with no flow to have no haematocrit
                                if (fabs(segments[jdx]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                {
                                    segments[jdx]->GetFlowProperties()->SetHaematocrit(0.0);
                                }   
                                // double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                                // double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                                // if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                                // {
                                //     max_difference = difference;
                                //     h_for_max = current_haematocrit;
                                //     prev_for_max = previous_haematocrit[jdx];
                                // }
                                // previous_haematocrit[jdx] = current_haematocrit;
                            }
                        //     std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                        //     // Print the final or intermediary convergence results
                        //     if(max_difference<=tolerance2)  
                        //     {
                        //         std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                        //         broken_solver = 0;
                        //         break;
                        //     }
                        //     else
                        //     {
                        //         if(idx%1==0)
                        //         {
                        //             std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                        //             std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                        //             std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                        //             p_network->Write(output_file);
                        //         }
                        //     }

                        //     // If there is no convergence after all the iterations, print the error message.
                        //     if(idx==max_iter-1)
                        //     {
                        //         std::cout << "Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha << std::endl;
                        //         error_log << "\n Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha; 
                        //         broken_solver = 1;
                        //         break;
                        //         // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                        //     }
                        // }

                        // If solver doesn't converge, move on to next one
                        if (broken_solver == 1)
                        {
                            continue;
                        }

                        // Run the simulation 
                        // SimulationTime::Instance()->SetStartTime(0.0);
                        // SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                        // auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                        // p_microvessel_solver->SetVesselNetwork(p_network);
                        // p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                        // p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                        // p_microvessel_solver->Run();
                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                        p_network->Write(output_file);

                        // Print the average oxygenation
                        // std::vector<double> solution = p_oxygen_solver->GetSolution();
                        // double average_oxygen = 0.0;
                        // for(unsigned jdx=0;jdx<solution.size();jdx++)
                        // {
                        //     average_oxygen += solution[jdx];
                        // }
                        // average_oxygen /= double(solution.size());
                        // std::cout << "Average oxygen: " << average_oxygen << std::endl;

                        // // Write the PQs
                        // QFlowRate threshold = 5.e-14*unit::metre_cubed_per_second;
                        QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;
                        // QFlowRate threshold = 2.e-13*unit::metre_cubed_per_second;
                        // // outfile.open("hexagonal_perfusion_quotients.txt", std::ios_base::app);
                        // // outfile << network_name << " " << solver_name << " " << selection_string << " " << heterogeneity_string << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                        // // outfile.close();  
                        PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
                        // QuotientMatrixAggregate[i][KilledVessels] += PerfusionQuotient;  // store the PQ with the associated NHet and number of kills
                        outfile.open("hex_lognormal_individual_perfusion_quotients.txt", std::ios_base::app);
                        outfile << network_name << " " << solver_name << " " << heterogeneity_string << " " << selection_string << " " << kill_string << " " << PerfusionQuotient << " \n"; 
                        outfile.close();

                        // Remove the smallest vessel
                        QLength minimum_radius = inlet_vessel_radius;
                        unsigned int minimum_index = 0;
                        for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)  // for all the segments in the network
                        {
                            // Get the current segment's radius
                            QLength current_radius = vessels[vessel_index]->GetRadius();

                            // If the current radius is less than the minimum radius, record the new minimum
                            if (current_radius < minimum_radius)
                            {
                                minimum_radius = current_radius;
                                minimum_index = vessel_index;
                            }
                        }
                        p_network->RemoveVessel(vessels[minimum_index], true);  // remove the vessel

                        // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                        ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                        ParameterCollection::Instance()->Destroy();
                        BaseUnits::Instance()->Destroy();
                        SimulationTime::Instance()->Destroy();
                    }
                }
            }
        }
         
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a full 2D hexagonal network on a PDE grid with flow, H-splitting, all (not just non-inlet/outlet) radii set according to a log normal distribution, and kills-based pruning on the non-inlet/outlet vessels (correcting for the unexplained -50 um offset in the other simulation)
    void xTestConstantOffsetHeterogeneousHexagonalNeighbourhoodWithIndividualPruningForNonInletOutletVesselsAndCompleteDiameterAssignment2D()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";
        
        // Set network in filename
        std::string network_name = "TestHexagonalNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("hex_lognormal_individual_perfusion_quotients.txt");
        outfile.close();

        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Define the key parameters
        // double dimless_domain_size_x = 2050.0;  // x-coordinate of output node
        // double dimless_domain_size_y = 1818.65 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
        // unsigned dimless_vessel_length = 100.0;
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        // double initial_haematocrit = 0.45;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        unsigned max_n_alpha = 3;
        unsigned n_vessels = 407;  // number of vessels from which to select ones to make thin
        // double percToKill = 0.2;  // percentage of vessels to kill
        // unsigned ToBeKilled = (unsigned)(percToKill*n_vessels);  // number to kill
        unsigned ToBeKilled = 200;  // number to kill

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=1; h_solver++)
        {
            // Initialise the simulation
            std::shared_ptr<VesselNetwork<2> > p_network;
            std::vector<std::shared_ptr<Vessel<2> > > vessels;
            VesselPtr<2> p_current_vessel;
            QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
            QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
            std::vector<std::vector<double> > QuotientMatrixMean(max_n_alpha,std::vector<double>(ToBeKilled+1));
            std::vector<std::vector<double> > QuotientMatrixAggregate(max_n_alpha,std::vector<double>(ToBeKilled+1,0.0));
            std::vector<std::vector<unsigned> > Order;
            double PerfusionQuotient;

            // Read the network layout from a file and store it in an array
            VesselNetworkGenerator<2> network_generator;
            std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix_Offset.txt");
            std::vector<std::vector<double> > rEdgesMatrix;
            string line;
            while (std::getline(in, line)) 
            {
                rEdgesMatrix.push_back(std::vector<double>());
                std::stringstream split(line);
                double value;
                while (split >> value)
                {
                    rEdgesMatrix.back().push_back(value);
                }
            }

            // Loop over increasing level of heteregeneity (% of thin vessels)
            for(unsigned n_alpha=0; n_alpha<=max_n_alpha; n_alpha++)
            {
                // Convert to actual mu
                int alpha = (3*n_alpha)+10;  // alpha determines the mean

                // Loop over different thin vessel layouts
                for(unsigned list_number=1; list_number<=thin_selections; list_number++)
                {                
                    // Choose a radii list file and read the file to an array
                    std::vector<std::vector<double> > radii_array;
                    radii_array.clear();
                    // int list_index=0;
                    string line_1;
                    std::ifstream radius_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/square_hexagonal_radius_log_normal_distribution/mu_"+to_string(alpha)+"/radii_list_"+to_string(list_number)+".txt");
                    while (std::getline(radius_list_file, line_1)) 
                    {
                        radii_array.push_back(std::vector<double>());
                        std::stringstream split(line_1);
                        double value;
                        while (split >> value)
                        {
                            radii_array.back().push_back(value);
                        }
                    }

                    // Choose the corresponding vessel ID order file and read the file to an array
                    std::vector<std::vector<double> > id_array;
                    id_array.clear();
                    // int list_index=0;
                    string line_2;
                    std::ifstream id_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/square_hexagonal_radius_log_normal_distribution/mu_"+to_string(alpha)+"/id_list_"+to_string(list_number)+".txt");
                    while (std::getline(id_list_file, line_2)) 
                    {
                        id_array.push_back(std::vector<double>());
                        std::stringstream split(line_2);
                        double value;
                        while (split >> value)
                        {
                            id_array.back().push_back(value);
                        }
                    }

                    // Generate the network
                    p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);

                    // Set inlet and outlet nodes (if node (start/end) is at the edge of the domain space, make it input or output depending on which side), and assign each non-inlet/outlet vessel a radius from the list
                    auto p_segment = p_network->GetVesselSegments()[0];
                    p_segment->SetRadius(inlet_vessel_radius);
                    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   
                    vessels = p_network->GetVessels();
                    for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
                    {
                        if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the start node equals one (no other segments at nodes)
                        {
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the start node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);  // set the node as the input node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the start node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the ouput node pressure
                            }
                        }
                        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the end node equals one (no other segments at nodes)	
                        {
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the end node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);;  // set the node as the input node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the end node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the output node pressure
                            }
                        }
                    }
                
                    // Assign vessel radii based on imported arrays
                    unsigned int vessel_index = 0;
                    while(vessel_index < n_vessels)
                    {
                        p_current_vessel = p_network->GetVessel(id_array[vessel_index][0]);

                        // Get segment radius from array
                        QLength new_vessel_radius(double(radii_array[vessel_index][0])*unit::microns);

                        p_current_vessel->SetRadius(new_vessel_radius);
                        vessel_index++;
                    }
                    vessels = p_network->GetVessels();

                    // Set the haematocrit solver
                    std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                    std::string solver_name;
                    if (h_solver==1)
                    {
                        solver_name = "ConstantHaematocrit"; 
                        std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==2)
                    {
                        solver_name = "PriesHaematocrit"; 
                        std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                
                    }
                    else if (h_solver==3)
                    {
                        solver_name = "MemoryHaematocrit"; 
                        std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==4)
                    {
                        solver_name = "FungHaematocrit"; 
                        std::cout << "Now using FungHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==5)
                    {
                        solver_name = "YangHaematocrit"; 
                        std::cout << "Now using YangHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();   
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }

                    // Set up the viscosity solver
                    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                    p_viscosity_calculator->SetVesselNetwork(p_network);
                    p_viscosity_calculator->Calculate();
                    
                    // Set up the impedance solver
                    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                    p_impedance_calculator->SetVesselNetwork(p_network);
                    p_impedance_calculator->Calculate();

                    // Set up the flow solver 
                    FlowSolver<2> flow_solver;
                    flow_solver.SetVesselNetwork(p_network);
                    // flow_solver.SetUp();
                    flow_solver.SetUseDirectSolver(true);
                    // flow_solver.Solve();

                    // Set up the grid for the finite difference solver
                    // auto p_grid = RegularGrid<2>::Create();
                    // p_grid->SetSpacing(grid_spacing);
                    // c_vector<unsigned, 3> dimensions;
                    // dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                    // dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                    // dimensions[2] = 1;
                    // p_grid->SetDimensions(dimensions);

                    // Choose the PDE
                    // auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                    
                    // Set the diffusivity and decay terms
                    // p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                    // p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                    // Set up the discrete source
                    // auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                    // QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                    //         GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                    // QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                    //         Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                    // p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                    // p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                    // p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    // p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                    // // Set up the finite difference solver for oxygen (which handles everything)
                    // auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                    // p_oxygen_solver->SetPde(p_oxygen_pde);
                    // p_oxygen_solver->SetLabel("oxygen");
                    // p_oxygen_solver->SetGrid(p_grid);
                    // solver.SetFileName("oxygen_solution_0");

                    // Prune all vessels up to specified dose 
                    for(unsigned KilledVessels=0; KilledVessels < ToBeKilled+1; KilledVessels++)
                    { 
                        // Display status message
                        std::cout << "Now killed " << KilledVessels << " vessels." << std::endl;  

                        // Set filename
                        std::stringstream selection_stream;
                        std::stringstream heterogeneity_stream;
                        std::stringstream kill_stream;
                        selection_stream << std::fixed << std::setprecision(0) << to_string(list_number);                    
                        heterogeneity_stream << std::fixed << std::setprecision(0) << alpha;                              
                        kill_stream << std::fixed << std::setprecision(0) << KilledVessels;           
                        std::string selection_string = selection_stream.str();
                        std::string heterogeneity_string = heterogeneity_stream.str();                    
                        std::string kill_string = kill_stream.str();                    
                        std::string file_name = network_name + solver_name + "/Mu" + heterogeneity_string + "/Selection" + selection_string + "/Kills" + kill_string;
                        std::string str_directory_name = file_name;
                        std::cout << str_directory_name << std::endl;
                        auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                        // Set up flag for broken solver
                        unsigned broken_solver = 0;

                        // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                        // unsigned max_iter = 5;  // 1000 
                        // double tolerance2 = 1.e-10;
                        std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                        // std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                        // for(unsigned idx=0; idx<max_iter; idx++)
                        // {
                            // Run the solvers
                            p_impedance_calculator->Calculate();
                            p_abstract_haematocrit_solver->Calculate();
                            p_viscosity_calculator->Calculate();
                            flow_solver.SetUp();
                            flow_solver.Solve();

                            // Get the residual
                            // double max_difference = 0.0;
                            // double h_for_max = 0.0;
                            // double prev_for_max = 0.0;
                            for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                            {
                                // Set segments with no flow to have no haematocrit
                                if (fabs(segments[jdx]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                {
                                    segments[jdx]->GetFlowProperties()->SetHaematocrit(0.0);
                                }   
                                // double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                                // double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                                // if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                                // {
                                //     max_difference = difference;
                                //     h_for_max = current_haematocrit;
                                //     prev_for_max = previous_haematocrit[jdx];
                                // }
                                // previous_haematocrit[jdx] = current_haematocrit;
                            }
                        //     std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                        //     // Print the final or intermediary convergence results
                        //     if(max_difference<=tolerance2)  
                        //     {
                        //         std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                        //         broken_solver = 0;
                        //         break;
                        //     }
                        //     else
                        //     {
                        //         if(idx%1==0)
                        //         {
                        //             std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                        //             std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                        //             std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                        //             p_network->Write(output_file);
                        //         }
                        //     }

                        //     // If there is no convergence after all the iterations, print the error message.
                        //     if(idx==max_iter-1)
                        //     {
                        //         std::cout << "Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha << std::endl;
                        //         error_log << "\n Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha; 
                        //         broken_solver = 1;
                        //         break;
                        //         // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                        //     }
                        // }

                        // If solver doesn't converge, move on to next one
                        if (broken_solver == 1)
                        {
                            continue;
                        }

                        // Run the simulation 
                        // SimulationTime::Instance()->SetStartTime(0.0);
                        // SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                        // auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                        // p_microvessel_solver->SetVesselNetwork(p_network);
                        // p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                        // p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                        // p_microvessel_solver->Run();
                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                        p_network->Write(output_file);

                        // Print the average oxygenation
                        // std::vector<double> solution = p_oxygen_solver->GetSolution();
                        // double average_oxygen = 0.0;
                        // for(unsigned jdx=0;jdx<solution.size();jdx++)
                        // {
                        //     average_oxygen += solution[jdx];
                        // }
                        // average_oxygen /= double(solution.size());
                        // std::cout << "Average oxygen: " << average_oxygen << std::endl;

                        // // Write the PQs
                        // QFlowRate threshold = 5.e-14*unit::metre_cubed_per_second;
                        QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;
                        // QFlowRate threshold = 2.e-13*unit::metre_cubed_per_second;
                        // // outfile.open("hexagonal_perfusion_quotients.txt", std::ios_base::app);
                        // // outfile << network_name << " " << solver_name << " " << selection_string << " " << heterogeneity_string << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                        // // outfile.close();  
                        PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
                        // QuotientMatrixAggregate[i][KilledVessels] += PerfusionQuotient;  // store the PQ with the associated NHet and number of kills
                        outfile.open("hex_lognormal_individual_perfusion_quotients.txt", std::ios_base::app);
                        outfile << network_name << " " << solver_name << " " << heterogeneity_string << " " << selection_string << " " << kill_string << " " << PerfusionQuotient << " \n"; 
                        outfile.close();

                        // Remove the smallest vessel
                        QLength minimum_radius = inlet_vessel_radius;
                        unsigned int minimum_index = 0;
                        for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)  // for all the segments in the network
                        {                           
                            // // Exclude inlets and outlets
                            if (vessels[vessel_index]->GetStartNode()->GetFlowProperties()->IsInputNode()==0
                            && vessels[vessel_index]->GetStartNode()->GetFlowProperties()->IsOutputNode()==0
                            && vessels[vessel_index]->GetEndNode()->GetFlowProperties()->IsInputNode()==0
                            && vessels[vessel_index]->GetEndNode()->GetFlowProperties()->IsOutputNode()==0)
                            {
                                // Get the current segment's radius
                                QLength current_radius = vessels[vessel_index]->GetRadius();

                                // If the current radius is less than the minimum radius, record the new minimum
                                if (current_radius < minimum_radius)
                                {
                                    minimum_radius = current_radius;
                                    minimum_index = vessel_index;
                                }     
                            }
                        }
                        p_network->RemoveVessel(vessels[minimum_index], true);  // remove the vessel

                        // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                        ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                        ParameterCollection::Instance()->Destroy();
                        BaseUnits::Instance()->Destroy();
                        SimulationTime::Instance()->Destroy();
                    }
                }
            }
        }
         
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // Make a full 2D hexagonal network on a PDE grid with flow, H-splitting, non-inlet/outlet radii set according to a log normal distribution based on biological networks, and kills-based pruning on the non-inlet/outlet vessels (correcting for the unexplained -50 um offset in the other simulation).
    void xTestConstantOffsetBiologicalHexagonalNeighbourhoodWithIndividualPruning2DPaper1()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";
        
        // Set network in filename
        std::string network_name = "TestHexagonalNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("/tmp/narain/testoutput/TestHexagonalNetwork/hex_lognormal_individual_perfusion_quotients.txt");
        outfile.close();

        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Define the key parameters
        // double dimless_domain_size_x = 2050.0;  // x-coordinate of output node
        // double dimless_domain_size_y = 1818.65 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
        // unsigned dimless_vessel_length = 100.0;
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        // double initial_haematocrit = 0.45;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        unsigned max_n_alpha = 2;
        // unsigned n_vessels = 407;  // number of vessels from which to select ones to make thin
        unsigned n_vessels = 386;  // number of non-inlet/outlet vessels from which to select ones to make thin
        // double percToKill = 0.2;  // percentage of vessels to kill
        // unsigned ToBeKilled = (unsigned)(percToKill*n_vessels);  // number to kill
        unsigned ToBeKilled = 200;  // number to kill
        // unsigned ToBeKilled = 2;  // number to kill

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=1; h_solver++)
        {
            // Initialise the simulation
            std::shared_ptr<VesselNetwork<2> > p_network;
            std::vector<std::shared_ptr<Vessel<2> > > vessels;
            VesselPtr<2> p_current_vessel;
            QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
            QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
            std::vector<std::vector<double> > QuotientMatrixMean(max_n_alpha,std::vector<double>(ToBeKilled+1));
            std::vector<std::vector<double> > QuotientMatrixAggregate(max_n_alpha,std::vector<double>(ToBeKilled+1,0.0));
            std::vector<std::vector<unsigned> > Order;
            double PerfusionQuotient;

            // Read the network layout from a file and store it in an array
            VesselNetworkGenerator<2> network_generator;
            std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix_Offset.txt");
            std::vector<std::vector<double> > rEdgesMatrix;
            string line;
            while (std::getline(in, line)) 
            {
                rEdgesMatrix.push_back(std::vector<double>());
                std::stringstream split(line);
                double value;
                while (split >> value)
                {
                    rEdgesMatrix.back().push_back(value);
                }
            }

            // Loop over increasing level of heteregeneity (% of thin vessels)
            for(unsigned n_alpha=0; n_alpha<=max_n_alpha; n_alpha++)
            {
                // Convert to actual mu
                std::vector<std::string> alphas{"22.76", "28.5", "33.64"};
                std::string alpha = alphas[n_alpha];  // alpha determines the mean

                // Loop over different thin vessel layouts
                for(unsigned list_number=1; list_number<=thin_selections; list_number++)
                {                
                    // Choose a radii list file and read the file to an array
                    std::vector<std::vector<double> > radii_array;
                    radii_array.clear();
                    // int list_index=0;
                    string line_1;
                    std::ifstream radius_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/hexagonal_diameter_log_normal_distribution_sigma_8.68/mu_"+alpha+"/radii_list_"+to_string(list_number)+".txt");
                    while (std::getline(radius_list_file, line_1)) 
                    {
                        radii_array.push_back(std::vector<double>());
                        std::stringstream split(line_1);
                        double value;
                        while (split >> value)
                        {
                            radii_array.back().push_back(value);
                        }
                    }

                    // Choose the corresponding vessel ID order file and read the file to an array
                    std::vector<std::vector<double> > id_array;
                    id_array.clear();
                    // int list_index=0;
                    string line_2;
                    std::ifstream id_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/hexagonal_diameter_log_normal_distribution_sigma_8.68/mu_"+alpha+"/id_list_"+to_string(list_number)+".txt");
                    while (std::getline(id_list_file, line_2)) 
                    {
                        id_array.push_back(std::vector<double>());
                        std::stringstream split(line_2);
                        double value;
                        while (split >> value)
                        {
                            id_array.back().push_back(value);
                        }
                    }

                    // Generate the network
                    p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);

                    // Set inlet and outlet nodes (if node (start/end) is at the edge of the domain space, make it input or output depending on which side), and assign each non-inlet/outlet vessel a radius from the list
                    auto p_segment = p_network->GetVesselSegments()[0];
                    p_segment->SetRadius(inlet_vessel_radius);
                    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   
                    vessels = p_network->GetVessels();
                    for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
                    {
                        if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the start node equals one (no other segments at nodes)
                        {
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the start node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);  // set the node as the input node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the start node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the ouput node pressure                                
                            }
                        }
                        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the end node equals one (no other segments at nodes)	
                        {
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the end node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);;  // set the node as the input node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the end node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the output node pressure
                            }
                        }
                    }
                
                    // Assign vessel radii based on imported arrays
                    unsigned int vessel_index = 0;
                    while(vessel_index < n_vessels)
                    {
                        p_current_vessel = p_network->GetVessel(id_array[vessel_index][0]);

                        // Get segment radius from array
                        QLength new_vessel_radius(double(radii_array[vessel_index][0])*unit::microns);

                        p_current_vessel->SetRadius(new_vessel_radius);
                        vessel_index++;
                    }
                    vessels = p_network->GetVessels();

                    // Set the haematocrit solver
                    std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                    std::string solver_name;
                    if (h_solver==1)
                    {
                        solver_name = "ConstantHaematocrit"; 
                        std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==2)
                    {
                        solver_name = "PriesHaematocrit"; 
                        std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                
                    }
                    else if (h_solver==3)
                    {
                        solver_name = "MemoryHaematocrit"; 
                        std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==4)
                    {
                        solver_name = "FungHaematocrit"; 
                        std::cout << "Now using FungHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==5)
                    {
                        solver_name = "YangHaematocrit"; 
                        std::cout << "Now using YangHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();   
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }

                    // Set up the viscosity solver
                    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                    p_viscosity_calculator->SetVesselNetwork(p_network);
                    p_viscosity_calculator->Calculate();
                    
                    // Set up the impedance solver
                    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                    p_impedance_calculator->SetVesselNetwork(p_network);
                    p_impedance_calculator->Calculate();

                    // Set up the flow solver 
                    FlowSolver<2> flow_solver;
                    flow_solver.SetVesselNetwork(p_network);
                    // flow_solver.SetUp();
                    flow_solver.SetUseDirectSolver(true);
                    // flow_solver.Solve();

                    // Set up the grid for the finite difference solver
                    // auto p_grid = RegularGrid<2>::Create();
                    // p_grid->SetSpacing(grid_spacing);
                    // c_vector<unsigned, 3> dimensions;
                    // dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                    // dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                    // dimensions[2] = 1;
                    // p_grid->SetDimensions(dimensions);

                    // Choose the PDE
                    // auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                    
                    // Set the diffusivity and decay terms
                    // p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                    // p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                    // Set up the discrete source
                    // auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                    // QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                    //         GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                    // QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                    //         Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                    // p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                    // p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                    // p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    // p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                    // // Set up the finite difference solver for oxygen (which handles everything)
                    // auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                    // p_oxygen_solver->SetPde(p_oxygen_pde);
                    // p_oxygen_solver->SetLabel("oxygen");
                    // p_oxygen_solver->SetGrid(p_grid);
                    // solver.SetFileName("oxygen_solution_0");

                    // Prune all vessels up to specified dose 
                    for(unsigned KilledVessels=0; KilledVessels < ToBeKilled+1; KilledVessels++)
                    { 
                        // Display status message
                        std::cout << "Now killed " << KilledVessels << " vessels." << std::endl;  

                        // Set filename
                        std::stringstream selection_stream;
                        std::stringstream heterogeneity_stream;
                        std::stringstream kill_stream;
                        selection_stream << std::fixed << std::setprecision(0) << to_string(list_number);                    
                        heterogeneity_stream << std::fixed << std::setprecision(0) << alpha;                              
                        kill_stream << std::fixed << std::setprecision(0) << KilledVessels;           
                        std::string selection_string = selection_stream.str();
                        std::string heterogeneity_string = heterogeneity_stream.str();                    
                        std::string kill_string = kill_stream.str();                    
                        std::string file_name = network_name + solver_name + "/Mu" + heterogeneity_string + "/Selection" + selection_string + "/Kills" + kill_string;
                        std::string str_directory_name = file_name;
                        std::cout << str_directory_name << std::endl;
                        auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                        // Set up flag for broken solver
                        unsigned broken_solver = 0;

                        // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                        // unsigned max_iter = 5;  // 1000 
                        // double tolerance2 = 1.e-10;
                        std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                        // std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                        // for(unsigned idx=0; idx<max_iter; idx++)
                        // {
                            // Run the solvers
                            p_impedance_calculator->Calculate();
                            p_abstract_haematocrit_solver->Calculate();
                            p_viscosity_calculator->Calculate();
                            flow_solver.SetUp();
                            flow_solver.Solve();

                            // Get the residual
                            // double max_difference = 0.0;
                            // double h_for_max = 0.0;
                            // double prev_for_max = 0.0;
                            for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                            {
                                // Set segments with no flow to have no haematocrit
                                if (fabs(segments[jdx]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                                {
                                    segments[jdx]->GetFlowProperties()->SetHaematocrit(0.0);
                                }   
                                // double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                                // double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                                // if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                                // {
                                //     max_difference = difference;
                                //     h_for_max = current_haematocrit;
                                //     prev_for_max = previous_haematocrit[jdx];
                                // }
                                // previous_haematocrit[jdx] = current_haematocrit;
                            }
                        //     std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                        //     // Print the final or intermediary convergence results
                        //     if(max_difference<=tolerance2)  
                        //     {
                        //         std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                        //         broken_solver = 0;
                        //         break;
                        //     }
                        //     else
                        //     {
                        //         if(idx%1==0)
                        //         {
                        //             std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                        //             std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                        //             std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                        //             p_network->Write(output_file);
                        //         }
                        //     }

                        //     // If there is no convergence after all the iterations, print the error message.
                        //     if(idx==max_iter-1)
                        //     {
                        //         std::cout << "Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha << std::endl;
                        //         error_log << "\n Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha; 
                        //         broken_solver = 1;
                        //         break;
                        //         // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                        //     }
                        // }

                        // If solver doesn't converge, move on to next one
                        if (broken_solver == 1)
                        {
                            continue;
                        }

                        // Run the simulation 
                        // SimulationTime::Instance()->SetStartTime(0.0);
                        // SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                        // auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                        // p_microvessel_solver->SetVesselNetwork(p_network);
                        // p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                        // p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                        // p_microvessel_solver->Run();
                        std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                        p_network->Write(output_file);

                        // Print the average oxygenation
                        // std::vector<double> solution = p_oxygen_solver->GetSolution();
                        // double average_oxygen = 0.0;
                        // for(unsigned jdx=0;jdx<solution.size();jdx++)
                        // {
                        //     average_oxygen += solution[jdx];
                        // }
                        // average_oxygen /= double(solution.size());
                        // std::cout << "Average oxygen: " << average_oxygen << std::endl;

                        // // Write the PQs
                        // QFlowRate threshold = 5.e-14*unit::metre_cubed_per_second;
                        // QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;
                        QFlowRate threshold = 3.e-13*unit::metre_cubed_per_second;
                        // // outfile.open("hexagonal_perfusion_quotients.txt", std::ios_base::app);
                        // // outfile << network_name << " " << solver_name << " " << selection_string << " " << heterogeneity_string << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                        // // outfile.close();  
                        PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
                        // QuotientMatrixAggregate[i][KilledVessels] += PerfusionQuotient;  // store the PQ with the associated NHet and number of kills
                        outfile.open("/tmp/narain/testoutput/TestHexagonalNetwork/hex_lognormal_individual_perfusion_quotients.txt", std::ios_base::app);
                        outfile << network_name << " " << solver_name << " " << heterogeneity_string << " " << selection_string << " " << kill_string << " " << PerfusionQuotient << " \n"; 
                        outfile.close();

                        // Remove the smallest vessel
                        vessels = p_network->GetVessels();
                        QLength minimum_radius = inlet_vessel_radius;
                        unsigned int minimum_index = 0;
                        for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)  // for all the segments in the network
                        {                            
                            // // Exclude inlets and outlets
                            // if (vessels[vessel_index]->GetStartNode()->GetFlowProperties()->IsInputNode()==0
                            // && vessels[vessel_index]->GetStartNode()->GetFlowProperties()->IsOutputNode()==0
                            // && vessels[vessel_index]->GetEndNode()->GetFlowProperties()->IsInputNode()==0
                            // && vessels[vessel_index]->GetEndNode()->GetFlowProperties()->IsOutputNode()==0)
                            // {
                                // Get the current segment's radius
                                QLength current_radius = vessels[vessel_index]->GetRadius();

                                // If the current radius is less than the minimum radius, record the new minimum
                                if (current_radius < minimum_radius)
                                {
                                    minimum_radius = current_radius;
                                    minimum_index = vessel_index;
                                }     
                            // }
                        }
                        p_network->RemoveVessel(vessels[minimum_index], true);  // remove the vessel

                        // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                        ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                        ParameterCollection::Instance()->Destroy();
                        BaseUnits::Instance()->Destroy();
                        SimulationTime::Instance()->Destroy();
                    }
                }
            }
        }
         
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

    // // Print the execution time
    // printf("Time taken: %.2fs\n", (double)(clock() - clkStart)/CLOCKS_PER_SEC);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Construction Site: Keep Out
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Make a full 2D hexagonal network on a PDE grid with flow, H-splitting, radii set according to a log normal distribution, and ordered pruning (correcting for the unexplained -50 um offset in the other simulation)
    void xOnlyPerfusionQuotientsTestOffsetHeterogeneousHexagonalNeighbourhoodWithOrderedPruning2D()
    {
        // Initialise error log
        std::ostringstream error_log;
        error_log << "\n The following simulations failed to converge: \n";
        
        // Set network in filename
        std::string network_name = "TestHexagonalNetwork/";

        // Create the output file for the PQs
        std::ofstream outfile;
        outfile.open("heterogeneous_ordered_hexagonal_perfusion_quotients.txt");
        outfile.close();

        // Set up the reference length for the simulation
        QLength reference_length(1.0_um);
        BaseUnits::Instance()->SetReferenceLengthScale(reference_length);

        // Define the key parameters
        // double dimless_domain_size_x = 2050.0;  // x-coordinate of output node
        // double dimless_domain_size_y = 1818.65 + 86.6025;  // y-coordinate of topmost vessel + y-coordinate of lowest vessel (offset from domain edge)
        // unsigned dimless_vessel_length = 100.0;
        QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;
        // double initial_haematocrit = 0.45;
        double tolerance = 0.001;  // for location of inlet/outlet nodes
        max_alpha = 4;
        QLength large_vessel_radius = 37.5_um;
        unsigned n_vessels = 386;  // number of non-inlet/outlet vessels from which to select ones to make thin
        double percToKill = 0.2;  // percentage of vessels to kill
        unsigned ToBeKilled = (unsigned)(percToKill*n_vessels);  // number to kill

        // Run the simulation with different solvers of interest
        for (unsigned h_solver=1; h_solver<=1; h_solver++)
        {
            // Initialise the simulation
            std::shared_ptr<VesselNetwork<2> > p_network;
            std::vector<std::shared_ptr<Vessel<2> > > vessels;
            // string line2;
            VesselPtr<2> p_current_vessel;
            // unsigned ToBeThin;
            // unsigned NThin;
            QLength domain_side_length_x(dimless_domain_size_x * unit::microns);
            QLength domain_side_length_y(dimless_domain_size_y * unit::microns);
            std::vector<std::shared_ptr<Vessel<2> > >::iterator vessel_iterator;
            // std::ofstream outfileMean;
            std::vector<std::vector<double> > QuotientMatrixMean(max_alpha,std::vector<double>(ToBeKilled+1));
            std::vector<std::vector<double> > QuotientMatrixAggregate(max_alpha,std::vector<double>(ToBeKilled+1,0.0));
            std::vector<std::vector<unsigned> > Order;
            double PerfusionQuotient;

            // Read the network layout from a file and store it in an array
            VesselNetworkGenerator<2> network_generator;
            std::ifstream in("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/EdgesMatrix_Offset.txt");
            std::vector<std::vector<double> > rEdgesMatrix;
            string line;
            while (std::getline(in, line)) 
            {
                rEdgesMatrix.push_back(std::vector<double>());
                std::stringstream split(line);
                double value;
                while (split >> value)
                {
                    rEdgesMatrix.back().push_back(value);
                }
            }

            // Loop over increasing level of heteregeneity (% of thin vessels)
            for(unsigned alpha=1; alpha<=max_alpha; alpha++)
            {
                // Loop over different thin vessel layouts
                for(unsigned list_number=1; list_number<=thin_selections; list_number++)
                {                
                    // Choose a radii list file and read the file to an array
                    std::vector<std::vector<double> > radii_array;
                    radii_array.clear();
                    // int list_index=0;
                    string line_1;
                    std::ifstream radius_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/hexagonal_radius_log_normal_distribution/sigma_"+to_string(alpha)+"/radii_list_"+to_string(list_number)+".txt");
                    while (std::getline(radius_list_file, line_1)) 
                    {
                        radii_array.push_back(std::vector<double>());
                        std::stringstream split(line_1);
                        double value;
                        while (split >> value)
                        {
                            radii_array.back().push_back(value);
                        }
                    }

                    // Choose the corresponding vessel ID order file and read the file to an array
                    std::vector<std::vector<double> > id_array;
                    id_array.clear();
                    // int list_index=0;
                    string line_2;
                    std::ifstream id_list_file("/home/narain/Chaste/projects/MicrovesselChaste/test/simulation/flow/"+to_string(dimless_vessel_length)+"VesselLength/hexagonal_radius_log_normal_distribution/sigma_"+to_string(alpha)+"/id_list_"+to_string(list_number)+".txt");
                    while (std::getline(id_list_file, line_2)) 
                    {
                        id_array.push_back(std::vector<double>());
                        std::stringstream split(line_2);
                        double value;
                        while (split >> value)
                        {
                            id_array.back().push_back(value);
                        }
                    }

                    // Generate the network
                    p_network = network_generator.GenerateNetworkFromMatrixWithID(rEdgesMatrix);

                    // Set inlet and outlet nodes (if node (start/end) is at the edge of the domain space, make it input or output depending on which side), and assign each non-inlet/outlet vessel a radius from the list
                    auto p_segment = p_network->GetVesselSegments()[0];
                    p_segment->SetRadius(large_vessel_radius);
                    VesselNetworkPropertyManager<2>::SetSegmentProperties(p_network, p_segment);   
                    vessels = p_network->GetVessels();
                    for (vessel_iterator = vessels.begin(); vessel_iterator != vessels.end(); vessel_iterator++)
                    {
                        if((*vessel_iterator)->GetStartNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the start node equals one (no other segments at nodes)
                        {
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the start node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsInputNode(true);  // set the node as the input node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetStartNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the start node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetStartNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the ouput node pressure
                            }
                        }
                        if((*vessel_iterator)->GetEndNode()->GetNumberOfSegments() == 1)  // if the number of segments attached to the end node equals one (no other segments at nodes)	
                        {
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] < 0.5*dimless_vessel_length+tolerance)  // if the end node is located approx. at the start of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsInputNode(true);;  // set the node as the input node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(3333.0*unit::pascals);  // set the input node pressure
                            }
                            if((*vessel_iterator)->GetEndNode()->rGetLocation().Convert(1_um)[0] > dimless_domain_size_x - tolerance)  // if the end node is located approx. at the end of the domain 
                            {
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetIsOutputNode(true);  // set the node as the output node
                                (*vessel_iterator)->GetEndNode()->GetFlowProperties()->SetPressure(2000.0*unit::pascals);  // set the output node pressure
                            }
                        }
                    }
                
                    // Assign vessel radii based on imported arrays
                    unsigned int vessel_index = 0;
                    while(vessel_index < n_vessels)
                    {
                        p_current_vessel = p_network->GetVessel(id_array[vessel_index][0]);

                        // Get segment radius from array
                        QLength new_vessel_radius(double(radii_array[vessel_index][0])*unit::microns);

                        p_current_vessel->SetRadius(new_vessel_radius);
                        vessel_index++;
                    }
                    vessels = p_network->GetVessels();

                    // Set the haematocrit solver
                    std::shared_ptr<AbstractHaematocritSolver<2>> p_abstract_haematocrit_solver;
                    std::string solver_name;
                    if (h_solver==1)
                    {
                        solver_name = "ConstantHaematocrit"; 
                        std::cout << "Now using ConstantHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();  
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==2)
                    {
                        solver_name = "PriesHaematocrit"; 
                        std::cout << "Now using PriesHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();     
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;                
                    }
                    else if (h_solver==3)
                    {
                        solver_name = "MemoryHaematocrit"; 
                        std::cout << "Now using PriesWithMemoryHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = PriesWithMemoryHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==4)
                    {
                        solver_name = "FungHaematocrit"; 
                        std::cout << "Now using FungHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = BetteridgeHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();    
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }
                    else if (h_solver==5)
                    {
                        solver_name = "YangHaematocrit"; 
                        std::cout << "Now using YangHaematocritSolver..." << std::endl;
                        auto p_haematocrit_solver = YangHaematocritSolver<2>::Create();
                        // auto p_haematocrit_solver = ConstantHaematocritSolver<2>::Create();   
                        p_haematocrit_solver->SetVesselNetwork(p_network);
                        p_haematocrit_solver->SetHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                        p_abstract_haematocrit_solver = p_haematocrit_solver;     
                    }

                    // Set up the viscosity solver
                    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
                    p_viscosity_calculator->SetPlasmaViscosity(viscosity);
                    p_viscosity_calculator->SetVesselNetwork(p_network);
                    p_viscosity_calculator->Calculate();
                    
                    // Set up the impedance solver
                    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
                    p_impedance_calculator->SetVesselNetwork(p_network);
                    p_impedance_calculator->Calculate();

                    // Set up the flow solver 
                    FlowSolver<2> flow_solver;
                    flow_solver.SetVesselNetwork(p_network);
                    // flow_solver.SetUp();
                    flow_solver.SetUseDirectSolver(true);
                    // flow_solver.Solve();

                    // Set up the grid for the finite difference solver
                    auto p_grid = RegularGrid<2>::Create();
                    p_grid->SetSpacing(grid_spacing);
                    c_vector<unsigned, 3> dimensions;
                    dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
                    dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
                    dimensions[2] = 1;
                    p_grid->SetDimensions(dimensions);

                    // // Choose the PDE
                    // auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
                    
                    // // Set the diffusivity and decay terms
                    // p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
                    // p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

                    // Set up the discrete source
                    // auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
                    // QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
                    //         GenericParameters::mpGasConcentrationAtStp->GetValue("User");
                    // QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
                    //         Owen11Parameters::mpReferencePartialPressure->GetValue("User");
                    // p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
                    // p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
                    // p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
                    // p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

                    // Set up the finite difference solver for oxygen (which handles everything)
                    // auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
                    // p_oxygen_solver->SetPde(p_oxygen_pde);
                    // p_oxygen_solver->SetLabel("oxygen");
                    // p_oxygen_solver->SetGrid(p_grid);
                    // solver.SetFileName("oxygen_solution_0");

                    // Prune all vessels up to specified dose 
                    for(unsigned KilledVessels=0; KilledVessels < ToBeKilled+1; KilledVessels++)
                    { 
                        // Display status message
                        std::cout << "Now killed " << KilledVessels << " vessels." << std::endl;  

                        // Set filename
                        std::stringstream selection_stream;
                        std::stringstream heterogeneity_stream;
                        std::stringstream kill_stream;
                        selection_stream << std::fixed << std::setprecision(0) << to_string(list_number);                    
                        heterogeneity_stream << std::fixed << std::setprecision(0) << alpha;                              
                        kill_stream << std::fixed << std::setprecision(0) << KilledVessels;           
                        std::string selection_string = selection_stream.str();
                        std::string heterogeneity_string = heterogeneity_stream.str();                    
                        std::string kill_string = kill_stream.str();                    
                        std::string file_name = network_name + solver_name + "/Sigma" + heterogeneity_string + "/Selection" + selection_string + "/Kills" + kill_string;
                        std::string str_directory_name = file_name;
                        std::cout << str_directory_name << std::endl;
                        auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);

                        // Set up flag for broken solver
                        // unsigned broken_solver = 0;

                        // Set up an iteration to solve the non-linear problem (haematocrit problem is coupled to flow problem via viscosity/impedance)
                        // unsigned max_iter = 1000;  // 1000 
                        // double tolerance2 = 1.e-10;
                        // std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
                        // std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
                        // for(unsigned idx=0; idx<max_iter; idx++)
                        // {
                        //     // Run the solvers
                            p_impedance_calculator->Calculate();
                            p_abstract_haematocrit_solver->Calculate();
                            p_viscosity_calculator->Calculate();
                            flow_solver.SetUp();
                            flow_solver.Solve();

                        //     // Get the residual
                        //     double max_difference = 0.0;
                        //     double h_for_max = 0.0;
                        //     double prev_for_max = 0.0;
                        //     for(unsigned jdx=0;jdx<segments.size();jdx++)  // for all the segments in the network
                        //     {
                        //         // Set segments with no flow to have no haematocrit
                        //         if (fabs(segments[jdx]->GetFlowProperties()->GetFlowRate()) <= 1.e-16 *unit::metre_cubed_per_second)
                        //         {
                        //             segments[jdx]->GetFlowProperties()->SetHaematocrit(0.0);
                        //         }   
                        //         double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();  // get haematocrit
                        //         double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);  // difference in haematocrit
                        //         if(difference>max_difference)  // get the max. diff b/w prev. and current H, the value of H, and the prev. H
                        //         {
                        //             max_difference = difference;
                        //             h_for_max = current_haematocrit;
                        //             prev_for_max = previous_haematocrit[jdx];
                        //         }
                        //         previous_haematocrit[jdx] = current_haematocrit;
                        //     }
                        //     std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;

                        //     // Print the final or intermediary convergence results
                        //     if(max_difference<=tolerance2)  
                        //     {
                        //         std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
                        //         break;
                        //     }
                        //     else
                        //     {
                        //         if(idx%1==0)
                        //         {
                        //             std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                        //             // std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                        //             // std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                        //             // p_network->Write(output_file);
                        //         }
                        //     }

                        //     // If there is no convergence after all the iterations, print the error message.
                        //     if(idx==max_iter-1)
                        //     {
                        //         std::cout << "Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha << std::endl;
                        //         error_log << "\n Problem encountered in hexagonal network with h_solver = " << to_string(h_solver) << " using selection = " << to_string(list_number) << " and alpha = " << alpha; 
                        //         broken_solver = 1;
                        //         break;
                        //         // EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
                        //     }
                        // }

                        // // If solver doesn't converge, move on to next one
                        // if (broken_solver == 1)
                        // {
                        //     continue;
                        // }

                        // // Run the simulation 
                        // SimulationTime::Instance()->SetStartTime(0.0);
                        // SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);  // let's just do 1 time step; will be steady state anyway
                        // auto p_microvessel_solver = MicrovesselSolver<2>::Create();
                        // p_microvessel_solver->SetVesselNetwork(p_network);
                        // p_microvessel_solver->SetOutputFileHandler(p_file_handler);
                        // p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
                        // p_microvessel_solver->Run();
                        // // std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
                        // // p_network->Write(output_file);

                        // // Print the average oxygenation
                        // std::vector<double> solution = p_oxygen_solver->GetSolution();
                        // double average_oxygen = 0.0;
                        // for(unsigned jdx=0;jdx<solution.size();jdx++)
                        // {
                        //     average_oxygen += solution[jdx];
                        // }
                        // average_oxygen /= double(solution.size());
                        // std::cout << "Average oxygen: " << average_oxygen << std::endl;

                        // // Write the PQs
                        QFlowRate threshold = 1.e-13*unit::metre_cubed_per_second;
                        // // outfile.open("hexagonal_perfusion_quotients.txt", std::ios_base::app);
                        // // outfile << network_name << " " << solver_name << " " << selection_string << " " << heterogeneity_string << " " << p_network->GetPerfusionQuotientBeta(threshold) << " \n"; 
                        // // outfile.close();  
                        PerfusionQuotient = p_network->GetPerfusionQuotientBeta(threshold);
                        // QuotientMatrixAggregate[i][KilledVessels] += PerfusionQuotient;  // store the PQ with the associated NHet and number of kills
                        outfile.open("heterogeneous_ordered_hexagonal_perfusion_quotients.txt", std::ios_base::app);
                        outfile << network_name << " " << solver_name << " " << heterogeneity_string << " " << selection_string << " " << kill_string << " " << PerfusionQuotient << " \n"; 
                        outfile.close();

                        // Remove the smallest vessel
                        QLength minimum_radius = large_vessel_radius;
                        unsigned int minimum_index = 0;
                        for(unsigned vessel_index=0; vessel_index<vessels.size(); vessel_index++)  // for all the segments in the network
                        {
                            // Get the current segment's radius
                            QLength current_radius = vessels[vessel_index]->GetRadius();

                            // If the current radius is less than the minimum radius, record the new minimum
                            if (current_radius < minimum_radius)
                            {
                                minimum_radius = current_radius;
                                minimum_index = vessel_index;
                            }
                        }
                        p_network->RemoveVessel(vessels[minimum_index], true);  // remove the vessel

                        // Dump our parameter collection to an xml file and, importantly, clear it for the next test
                        ParameterCollection::Instance()->DumpToFile(p_file_handler->GetOutputDirectoryFullPath() + "parameter_collection.xml");
                        ParameterCollection::Instance()->Destroy();
                        BaseUnits::Instance()->Destroy();
                        SimulationTime::Instance()->Destroy();
                    }
                }
            }
        }
         
        // Print the error log
        std::string error_message = error_log.str();
        std::cout << error_message << std::endl; 
    }

};

#endif /*TESTHAEMATOCRITSOLVERS_HPP_*/
