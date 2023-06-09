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

#ifndef TESTWITHORWITHOUTMEMORY_HPP
#define TESTWITHORWITHOUTMEMORY_HPP

#include <cxxtest/TestSuite.h>
#include <boost/lexical_cast.hpp>
#include "VesselImpedanceCalculator.hpp"
#include "StructuralAdaptationSolver.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "VesselNetworkGenerator.hpp"
#include "FlowSolver.hpp"
#include "SimulationTime.hpp"
#include "PriesHaematocritSolver.hpp"
#include "PriesWithMemoryHaematocritSolver.hpp"
#include "UnitCollection.hpp"
#include "RegularGrid.hpp"
#include "SimulationTime.hpp"
#include "MicrovesselSolver.hpp"
#include "SimpleLinearEllipticFiniteDifferenceSolver.hpp"
#include "VesselBasedDiscreteSource.hpp"
#include "DiscreteContinuumLinearEllipticPde.hpp"
#include "Owen11Parameters.hpp"
#include "ViscosityCalculator.hpp"
#include "Secomb04Parameters.hpp"
#include "GenericParameters.hpp"
#include "StructuralAdaptationSolver.hpp"
#include "VesselNetworkPropertyManager.hpp"
#include "VesselNetworkGeometryCalculator.hpp"

#include <sstream>

#include "PetscAndVtkSetupAndFinalize.hpp"

class TestWithOrWithourMemory : public CxxTest::TestSuite
{


public:
    /** The following is to test that the Pries (without memory) lambda=4 figure can be faithfully reproduced.
        See bottom of RunNoCellsDichotomousWithOrWithoutMemoryEffects() for where this is called.*/
    void  VerifySolutionLambdaEquals4Pries(std::shared_ptr<VesselNetwork<2> > pNetwork, std::vector<double>& rOxygenSolution)
    {
        /* All vessels get 0.45 haematocrit.  This is reasonable because all radius reductions are symmetric.*/
        // (middle)
        TS_ASSERT_DELTA(0.45, pNetwork->GetVesselSegments()[61]->GetFlowProperties()->GetHaematocrit(), 1e-6);
        TS_ASSERT_DELTA(0.45, pNetwork->GetVesselSegments()[120]->GetFlowProperties()->GetHaematocrit(), 1e-6);
        TS_ASSERT_DELTA(0.45, pNetwork->GetVesselSegments()[81]->GetFlowProperties()->GetHaematocrit(), 1e-6);
        TS_ASSERT_DELTA(0.45, pNetwork->GetVesselSegments()[100]->GetFlowProperties()->GetHaematocrit(), 1e-6);

        /* Flow splitting (symmetric) */
        TS_ASSERT_DELTA(5.17974e-10, pNetwork->GetVesselSegments()[124]->GetFlowProperties()->GetFlowRate()/unit::metre_cubed_per_second, 1e-15);
        TS_ASSERT_DELTA(2.58987e-10, pNetwork->GetVesselSegments()[0]->GetFlowProperties()->GetFlowRate()/unit::metre_cubed_per_second, 1e-15);
        TS_ASSERT_DELTA(2.58987e-10, pNetwork->GetVesselSegments()[1]->GetFlowProperties()->GetFlowRate()/unit::metre_cubed_per_second, 1e-15);
        TS_ASSERT_DELTA(1.29493e-10, pNetwork->GetVesselSegments()[4]->GetFlowProperties()->GetFlowRate()/unit::metre_cubed_per_second, 1e-15);
        TS_ASSERT_DELTA(1.29493e-10, pNetwork->GetVesselSegments()[5]->GetFlowProperties()->GetFlowRate()/unit::metre_cubed_per_second, 1e-15);
        TS_ASSERT_DELTA(1.29493e-10, pNetwork->GetVesselSegments()[8]->GetFlowProperties()->GetFlowRate()/unit::metre_cubed_per_second, 1e-15);
        TS_ASSERT_DELTA(1.29493e-10, pNetwork->GetVesselSegments()[9]->GetFlowProperties()->GetFlowRate()/unit::metre_cubed_per_second, 1e-15);

        double average_oxygen = 0.0;
        for(unsigned jdx=0;jdx<rOxygenSolution.size();jdx++)
        {
            average_oxygen += rOxygenSolution[jdx];
        }
        average_oxygen /= rOxygenSolution.size();

        // Average oxygen
        TS_ASSERT_DELTA(16853.2187, average_oxygen, 1.0);
        // Low oxygen
        TS_ASSERT_DELTA(385.507, rOxygenSolution[0], 1.0);
        TS_ASSERT_DELTA(386.01,  rOxygenSolution[1], 1.0);
        // High oxygen
        TS_ASSERT_DELTA(27678.2, rOxygenSolution[16285], 1.0);
        TS_ASSERT_DELTA(27678.2, rOxygenSolution[28585], 1.0);
        TS_ASSERT_DELTA(27678.2, rOxygenSolution[8492], 1.0);
    }

    /** The following is to test that the "with memory" lambda=4 figure can be faithfully reproduced. Or at least that
        the asymmetry in the middle of the network will be visible.
        See bottom of RunNoCellsDichotomousWithOrWithoutMemoryEffects() for where this is called.*/
    void  VerifySolutionLambdaEquals4WithMemory(std::shared_ptr<VesselNetwork<2> > pNetwork, std::vector<double>& rOxygenSolution)
    {
        /* Note that these tests may be too fragile - 6 decimal places is more than enough to give the correct figure.*/
        // Lowest haematocrit segments (in middle at top/bottom)
        TS_ASSERT_DELTA(0.313475, pNetwork->GetVesselSegments()[61]->GetFlowProperties()->GetHaematocrit(), 1e-6);
        TS_ASSERT_DELTA(0.313475, pNetwork->GetVesselSegments()[63]->GetFlowProperties()->GetHaematocrit(), 1e-6);

        TS_ASSERT_DELTA(0.313476, pNetwork->GetVesselSegments()[120]->GetFlowProperties()->GetHaematocrit(), 1e-6);
        TS_ASSERT_DELTA(0.313476, pNetwork->GetVesselSegments()[122]->GetFlowProperties()->GetHaematocrit(), 1e-6);

        // Highest haematocrit segments (in middle and asymmetric)
        TS_ASSERT_DELTA(0.659911, pNetwork->GetVesselSegments()[81]->GetFlowProperties()->GetHaematocrit(), 1e-6);
        TS_ASSERT_DELTA(0.659911, pNetwork->GetVesselSegments()[83]->GetFlowProperties()->GetHaematocrit(), 1e-6);

        TS_ASSERT_DELTA(0.659897, pNetwork->GetVesselSegments()[100]->GetFlowProperties()->GetHaematocrit(), 1e-6);
        TS_ASSERT_DELTA(0.659897, pNetwork->GetVesselSegments()[102]->GetFlowProperties()->GetHaematocrit(), 1e-6);

        /* Flow splitting (non symmetric) */
        TS_ASSERT_DELTA(2.23815e-10, pNetwork->GetVesselSegments()[124]->GetFlowProperties()->GetFlowRate()/unit::metre_cubed_per_second, 1e-15);
        TS_ASSERT_DELTA(1.11895e-10, pNetwork->GetVesselSegments()[0]->GetFlowProperties()->GetFlowRate()/unit::metre_cubed_per_second, 1e-15);
        TS_ASSERT_DELTA(1.1192e-10, pNetwork->GetVesselSegments()[1]->GetFlowProperties()->GetFlowRate()/unit::metre_cubed_per_second, 1e-15);
        TS_ASSERT_DELTA(5.07894e-11, pNetwork->GetVesselSegments()[4]->GetFlowProperties()->GetFlowRate()/unit::metre_cubed_per_second, 1e-15);
        TS_ASSERT_DELTA(6.11309e-11, pNetwork->GetVesselSegments()[5]->GetFlowProperties()->GetFlowRate()/unit::metre_cubed_per_second, 1e-15);
        TS_ASSERT_DELTA(6.11001e-11, pNetwork->GetVesselSegments()[8]->GetFlowProperties()->GetFlowRate()/unit::metre_cubed_per_second, 1e-15);
        TS_ASSERT_DELTA(5.07949e-11, pNetwork->GetVesselSegments()[9]->GetFlowProperties()->GetFlowRate()/unit::metre_cubed_per_second, 1e-15);

        double average_oxygen = 0.0;
        for(unsigned jdx=0;jdx<rOxygenSolution.size();jdx++)
        {
            average_oxygen += rOxygenSolution[jdx];
        }
        average_oxygen /= rOxygenSolution.size();

        // Average oxygen
        TS_ASSERT_DELTA(17097.7, average_oxygen, 1.0);
        // Low oxygen
        TS_ASSERT_DELTA(381.278, rOxygenSolution[0], 1.0);
        TS_ASSERT_DELTA(381.76,  rOxygenSolution[1], 1.0);
        // High oxygen
        TS_ASSERT_DELTA(40581.4, rOxygenSolution[21530], 1.0);
        TS_ASSERT_DELTA(40581.4, rOxygenSolution[21505], 1.0);
    }



void RunNoCellsDichotomousWithOrWithoutMemoryEffects(bool withMemory)
{
    // order of the dichotomous network
    unsigned order=5;

    double dimless_length = 1.0;

    for(unsigned i_aux=1; i_aux<order+1; i_aux++)
    	{
	dimless_length += pow(2.0,-1/3)*sqrt(pow(2.0,-2.0*double(i_aux-1)/3.0)-pow(0.9,2)*pow(2.0, -2.0*double(i_aux-1)));
    	}
    QLength input_radius = 50_um;

    QDynamicViscosity viscosity = 1.e-3*unit::poiseuille;

    double inlet_haematocrit = 0.45;
    double initial_haematocrit = 0.45;

    // Generate the network

    VesselNetworkGenerator<2> network_generator;

   double lambda;
   double twicelambda;

   for (unsigned k_aux=1; k_aux<5; k_aux++)
   {
    lambda = 2.0+double(k_aux)*2.0;
        // lambda is quotient between the length and diameter...in vessel network generator, we use twice this value as an input parameter
    twicelambda = 2.0*lambda;


    // Length of the vertical projection of first-order vessels
    QLength main_vert_length = 0.9*twicelambda*input_radius*pow(2.0,-1.0/3.0);
    // vertical size of the domain
    QLength domain_side_length_y = 4.0*main_vert_length;

    std::ostringstream strs;
    if (withMemory)
    {
        strs << "Dichotomous_NoCorners_NewModel_LambdaEquals" << lambda;
    }
    else
    {
        strs << "Dichotomous_NoCorners_OldModel_LambdaEquals" << lambda;
    }
    std::string str_directory_name = strs.str();
    // horizontal size of the domain
    QLength domain_side_length_x = dimless_length*2.0*twicelambda*input_radius;
    auto p_file_handler = std::make_shared<OutputFileHandler>(str_directory_name, true);
    std::shared_ptr<VesselNetwork<2> > p_network = network_generator.GenerateForkingNetworkNoCorners(order, main_vert_length, input_radius, twicelambda);

    // identify input and output nodes and assign them properties

    VesselNodePtr<2> p_inlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
            Vertex<2>(0.0_um,2.0*main_vert_length));
    VesselNodePtr<2> p_outlet_node = VesselNetworkGeometryCalculator<2>::GetNearestNode(p_network,
            Vertex<2>(domain_side_length_x, 2.0*main_vert_length));
    p_inlet_node->GetFlowProperties()->SetIsInputNode(true);
    p_inlet_node->GetFlowProperties()->SetPressure(3320.0_Pa);
    p_outlet_node->GetFlowProperties()->SetIsOutputNode(true);
    p_outlet_node->GetFlowProperties()->SetPressure(2090.0_Pa);

    //grid for finite difference solver
    auto p_grid = RegularGrid<2>::Create();
    QLength grid_spacing = 10_um;
    p_grid->SetSpacing(grid_spacing);
    c_vector<unsigned, 3> dimensions;
    dimensions[0] = unsigned((domain_side_length_x)/(grid_spacing))+1; // num x
    dimensions[1] = unsigned((domain_side_length_y)/(grid_spacing))+1; // num_y
    dimensions[2] = 1;
    p_grid->SetDimensions(dimensions);

    /**
     * Next set up the PDE for oxygen. Cells will act as a continuum oxygen sink
     */
    auto p_oxygen_pde = DiscreteContinuumLinearEllipticPde<2>::Create();
    p_oxygen_pde->SetIsotropicDiffusionConstant(Owen11Parameters::mpOxygenDiffusivity->GetValue("User"));
    p_oxygen_pde->SetContinuumLinearInUTerm(-1.0*Owen11Parameters::mpCellOxygenConsumptionRate->GetValue("User"));

    /**
    * Vessels release oxygen depending on their haematocrit levels (discrete sinks); see also the description in the supplementary material in our paper on the memory effects
    */
    auto p_vessel_oxygen_source = VesselBasedDiscreteSource<2>::Create();
    QSolubility oxygen_solubility_at_stp = Secomb04Parameters::mpOxygenVolumetricSolubility->GetValue("User") *
            GenericParameters::mpGasConcentrationAtStp->GetValue("User");
    QConcentration vessel_oxygen_concentration = oxygen_solubility_at_stp *
            Owen11Parameters::mpReferencePartialPressure->GetValue("User");
    p_vessel_oxygen_source->SetReferenceConcentration(vessel_oxygen_concentration);
    p_vessel_oxygen_source->SetVesselPermeability(Owen11Parameters::mpVesselOxygenPermeability->GetValue("User"));
    p_vessel_oxygen_source->SetReferenceHaematocrit(Owen11Parameters::mpInflowHaematocrit->GetValue("User"));
    p_oxygen_pde->AddDiscreteSource(p_vessel_oxygen_source);

    /*
    * Set up a finite difference solver and pass it the pde and grid.
    */
    auto p_oxygen_solver = SimpleLinearEllipticFiniteDifferenceSolver<2>::Create();
    p_oxygen_solver->SetPde(p_oxygen_pde);
    p_oxygen_solver->SetLabel("oxygen");
    p_oxygen_solver->SetGrid(p_grid);


    // Switch between solvers for Pries or newer "with memory"
    std::shared_ptr<AbstractHaematocritSolver<2>> p_abs_haematocrit_calculator;
    if (withMemory)
    {
        auto p_haematocrit_calculator = PriesWithMemoryHaematocritSolver<2>::Create();
        p_haematocrit_calculator->SetVesselNetwork(p_network);
        p_haematocrit_calculator->SetHaematocrit(inlet_haematocrit);
        p_abs_haematocrit_calculator = p_haematocrit_calculator;
    }
    else
    {
        auto  p_haematocrit_calculator = PriesHaematocritSolver<2>::Create();
        p_haematocrit_calculator->SetVesselNetwork(p_network);
        p_haematocrit_calculator->SetHaematocrit(inlet_haematocrit);
        p_abs_haematocrit_calculator = p_haematocrit_calculator;
    }
    auto p_impedance_calculator = VesselImpedanceCalculator<2>::Create();
    auto p_viscosity_calculator = ViscosityCalculator<2>::Create();
    p_viscosity_calculator->SetPlasmaViscosity(viscosity);

    p_impedance_calculator->SetVesselNetwork(p_network);
    p_viscosity_calculator->SetVesselNetwork(p_network);
    p_viscosity_calculator->Calculate();
    p_impedance_calculator->Calculate();

    FlowSolver<2> flow_solver;
    flow_solver.SetVesselNetwork(p_network);
    flow_solver.SetUp();

    unsigned max_iter = 1000;
    double tolerance2 = 1.e-10;

    std::vector<VesselSegmentPtr<2> > segments = p_network->GetVesselSegments();
    std::vector<double> previous_haematocrit(segments.size(), double(initial_haematocrit));
    // iteration to solve the nonlinear problem follows (haematocrit problem is coupled to the flow problem via viscosity/impedance)
    for(unsigned idx=0;idx<max_iter;idx++)
    {
        p_impedance_calculator->Calculate();
        flow_solver.SetUp();
        flow_solver.Solve();
        p_abs_haematocrit_calculator->Calculate();
        p_viscosity_calculator->Calculate();
        // Get the residual
        double max_difference = 0.0;
        double h_for_max = 0.0;
        double prev_for_max = 0.0;
        for(unsigned jdx=0;jdx<segments.size();jdx++)
        {
            double current_haematocrit = segments[jdx]->GetFlowProperties()->GetHaematocrit();
            double difference = std::abs(current_haematocrit - previous_haematocrit[jdx]);
            if(difference>max_difference)
            {
                max_difference = difference;
                h_for_max = current_haematocrit;
                prev_for_max = previous_haematocrit[jdx];
            }
            previous_haematocrit[jdx] = current_haematocrit;
        }
        std::cout << "H at max difference: " << h_for_max << ", Prev H at max difference:" << prev_for_max << std::endl;
        if(max_difference<=tolerance2)
        {
            std::cout << "Converged after: " << idx << " iterations. " <<  std::endl;
            break;
        }
        else
        {
            // Output intermediate results
            if(idx%1==0)
            {
                std::cout << "Max Difference at iter: " << idx << " is " << max_difference << std::endl;
                std::string file_suffix = "IntermediateHaematocrit_" + std::to_string(idx) + ".vtp";
                std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append(file_suffix);
                p_network->Write(output_file);
            }
        }

        if(idx==max_iter-1)
        {
            EXCEPTION("Did not converge after " + std::to_string(idx) + " iterations.");
        }
    }

    SimulationTime::Instance()->SetStartTime(0.0);
// Let's just do 1 time step; will be steady state anyway
    SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
    auto p_microvessel_solver = MicrovesselSolver<2>::Create();
    p_microvessel_solver->SetVesselNetwork(p_network);
    p_microvessel_solver->SetOutputFileHandler(p_file_handler);
    p_microvessel_solver->AddDiscreteContinuumSolver(p_oxygen_solver);
    p_microvessel_solver->Run();

    std::vector<double> solution = p_oxygen_solver->GetSolution();
    double average_oxygen = 0.0;
    for(unsigned jdx=0;jdx<solution.size();jdx++)
    {
        average_oxygen += solution[jdx];
    }
    average_oxygen /= double(solution.size());
    std::cout << "Average oxygen: " << average_oxygen << std::endl;
    std::string output_file = p_file_handler->GetOutputDirectoryFullPath().append("FinalHaematocrit.vtp");
    p_network->Write(output_file);
    SimulationTime::Instance()->Destroy();
    
    // Test that the "with memory" lambda=4 figure can be faithfully reproduced.
    if (fabs(lambda-4.0)<1e-1)
    {
        if (withMemory)
        {
            VerifySolutionLambdaEquals4WithMemory(p_network, solution);
        }
        else
        {
            VerifySolutionLambdaEquals4Pries(p_network, solution);
        }
    }


    }
}

    void TestNoCellsDichotomousWithOrWithoutMemoryEffects()
    {
        RunNoCellsDichotomousWithOrWithoutMemoryEffects(false);
        RunNoCellsDichotomousWithOrWithoutMemoryEffects(true);
    }

};

#endif // TESTWITHORWITHOUTMEMORY_HPP
