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

#include "VesselSegment.hpp"
#include "UnitCollection.hpp"
#include "Vertex.hpp"
#include "DistanceMap.hpp"
#include "PetscTools.hpp"
#include "PetscVecTools.hpp"
#include "ReplicatableVector.hpp"

template<unsigned DIM>
DistanceMap<DIM>::DistanceMap()
    :   AbstractRegularGridDiscreteContinuumSolver<DIM>(),
        mUseSegmentRadii(false)
{
    this->mLabel = "Distance Map";
}

template<unsigned DIM>
std::shared_ptr<DistanceMap<DIM> > DistanceMap<DIM>::Create()
{
    return std::make_shared<DistanceMap>();
}

template<unsigned DIM>
DistanceMap<DIM>::~DistanceMap()
{

}

template<unsigned DIM>
void DistanceMap<DIM>::SetUseSegmentRadii(bool useRadii)
{
    this->mUseSegmentRadii = useRadii;
}

template<unsigned DIM>
void DistanceMap<DIM>::Solve()
{
    // Set up the vtk solution on the regular grid
    this->Setup();

    unsigned num_points = this->mpDensityMap->GetGridCalculator()->GetGrid()->GetNumberOfPoints();
    std::vector<double> distances(num_points, 0.0);

    std::vector<std::shared_ptr<VesselSegment<DIM> > > segments;
    segments = this->mpDensityMap->GetVesselNetwork()->GetVesselSegments();
    QLength ref_length = this->mpDensityMap->GetGridCalculator()->GetGrid()->GetReferenceLengthScale();

    for(unsigned idx=0;idx<num_points; idx++)
    {
        Vertex<DIM> location = this->mpDensityMap->GetGridCalculator()->GetGrid()->GetPoint(idx);
        QLength min_distance = DBL_MAX * unit::metres;
        for (unsigned jdx = 0; jdx <  segments.size(); jdx++)
        {
            QLength seg_dist = segments[jdx]->GetDistance(location);
            if(this->mUseSegmentRadii && seg_dist<=segments[jdx]->GetRadius())
            {
                seg_dist = 0_m;
            }
            if(seg_dist < min_distance)
            {
                min_distance = seg_dist;
            }
        }
        distances[idx] = min_distance/ref_length;
    }

    this->UpdateSolution(distances);
    if (this->mWriteSolution)
    {
        this->Write();
    }
}

// Explicit instantiation
template class DistanceMap<2> ;
template class DistanceMap<3> ;
