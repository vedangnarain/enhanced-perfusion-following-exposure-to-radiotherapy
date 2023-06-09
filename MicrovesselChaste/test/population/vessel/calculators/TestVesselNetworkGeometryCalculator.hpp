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

#ifndef TESTVESSELNETWORKGEOMETRYCALCULATOR_HPP_
#define TESTVESSELNETWORKGEOMETRYCALCULATOR_HPP_

#include <cxxtest/TestSuite.h>
#include "VesselNode.hpp"
#include "SmartPointers.hpp"
#include "ChastePoint.hpp"
#include "VesselSegment.hpp"
#include "VesselNetwork.hpp"
#include "OutputFileHandler.hpp"
#include "UblasIncludes.hpp"
#include "VesselNetworkGenerator.hpp"
#include "UnitCollection.hpp"
#include "VesselNetworkGeometryCalculator.hpp"
#include "VesselNetworkPropertyManager.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestVesselNetworkGeometryCalculator : public CxxTest::TestSuite
{
public:

    void TestSimpleNetwork()
    {
        // Make some nodes
        std::vector<std::shared_ptr<VesselNode<3> > > nodes;
        nodes.push_back(VesselNode<3>::Create(1.0_um, 2.0_um, 6.0_um));
        nodes.push_back(VesselNode<3>::Create(3.0_um, 4.0_um, 7.0_um));
        nodes.push_back(VesselNode<3>::Create(3.0_um, 4.0_um, 8.0_um));
        nodes.push_back(VesselNode<3>::Create(3.0_um, 4.0_um, 9.0_um));

        // Make some vessels
        std::shared_ptr<Vessel<3> > pVessel1(Vessel<3>::Create(nodes[0], nodes[1]));
        std::shared_ptr<Vessel<3> > pVessel2(Vessel<3>::Create(nodes[1], nodes[2]));
        std::shared_ptr<Vessel<3> > pVessel3(Vessel<3>::Create(nodes[2], nodes[3]));

        std::vector<std::shared_ptr<Vessel<3> > > vessels;
        vessels.push_back(pVessel1);
        vessels.push_back(pVessel2);
        vessels.push_back(pVessel3);

        // Make a network
        std::shared_ptr<VesselNetwork<3> > p_network = VesselNetwork<3>::Create();
        p_network->AddVessels(vessels);

        std::shared_ptr<VesselNetworkGeometryCalculator<3> > p_calculator = VesselNetworkGeometryCalculator<3>::Create();
        TS_ASSERT_DELTA(p_calculator->GetAverageInterSegmentDistance(p_network)/1_m, 1.24402e-06, 1.e-3);
        TS_ASSERT_DELTA(p_calculator->GetAverageVesselLength(p_network)/1_m, 1.66667e-06, 1.e-8);
        TS_ASSERT_DELTA(p_calculator->GetTotalLength(p_network)/1_m, 5.e-6, 1.e-8);
        TS_ASSERT_DELTA(p_calculator->GetTotalSurfaceArea(p_network)/(1.0*unit::metres_squared), 3.14159e-10, 1.e-12);
        TS_ASSERT_DELTA(p_calculator->GetTotalVolume(p_network)/(1.0*unit::metres_cubed), 1.5708e-15, 1.e-16);
    }
};

#endif /*TESTVESSELNETWORKGEOMETRYCALCULATOR_HPP_*/
