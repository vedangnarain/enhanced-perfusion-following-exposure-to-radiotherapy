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

#ifndef TESTVESSELIMPEDANCECALCULATOR_HPP_
#define TESTVESSELIMPEDANCECALCULATOR_HPP_

#include <cxxtest/TestSuite.h>
#include "VesselImpedanceCalculator.hpp"
#include "SmartPointers.hpp"
#include "MathsCustomFunctions.hpp"
#include "UnitCollection.hpp"
#include "VesselNetworkPropertyManager.hpp"

#include "PetscAndVtkSetupAndFinalize.hpp"

class TestVesselImpedanceCalculator : public CxxTest::TestSuite
{

    typedef std::shared_ptr<VesselNode<2> > NodePtr2;
    typedef std::shared_ptr<VesselNode<3> > NodePtr3;
    typedef std::shared_ptr<VesselSegment<2> > SegmentPtr2;
    typedef std::shared_ptr<VesselSegment<3> > SegmentPtr3;
    typedef std::shared_ptr<Vessel<2> > VesselPtr2;
    typedef std::shared_ptr<Vessel<3> > VesselPtr3;

public:

    void TestFlowThroughSingleSegment()
    {

        // Make some nodes
        std::vector<ChastePoint<3> > points;
        points.push_back(ChastePoint<3>(0, 0, 0));
        points.push_back(ChastePoint<3>(5, 0, 0));

        std::vector<NodePtr3> nodes;
        nodes.push_back(NodePtr3(VesselNode<3>::Create(0_um)));
        nodes.push_back(NodePtr3(VesselNode<3>::Create(5_um)));

        SegmentPtr3 p_segment(VesselSegment<3>::Create(nodes[0], nodes[1]));
        VesselPtr3 p_vessel(Vessel<3>::Create(p_segment));

        // Generate the network
        std::shared_ptr<VesselNetwork<3> > p_vascular_network(new VesselNetwork<3>());

        p_vascular_network->AddVessel(p_vessel);
        QDynamicViscosity viscosity = 2e-3 * unit::poiseuille;
        QLength radius = 5.0_um;

        p_segment->SetRadius(radius);
        p_segment->GetFlowProperties()->SetViscosity(viscosity);

        VesselNetworkPropertyManager<3>::SetSegmentProperties(p_vascular_network, p_segment);

        VesselImpedanceCalculator<3> calculator;
        calculator.SetVesselNetwork(p_vascular_network);
        calculator.Calculate();

        QVolume term = M_PI * Qpow4(radius)/1_um;
        QFlowImpedance expected_impedance = 8.0 * 5.0 * (viscosity/term);
        TS_ASSERT_DELTA(p_vessel->GetFlowProperties()->GetImpedance()/expected_impedance, 1.0, 1e-6);
        TS_ASSERT_DELTA(p_segment->GetFlowProperties()->GetImpedance()/expected_impedance, 1.0, 1e-6);
    }
};

#endif /*TESTVESSELIMPEDANCECALCULATOR_HPP_*/
