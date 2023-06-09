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

#include "VesselNetworkGeometryCalculator.hpp"
#include "VesselNetworkCellPopulationInteractor.hpp"

template<unsigned DIM>
VesselNetworkCellPopulationInteractor<DIM>::VesselNetworkCellPopulationInteractor() :
        mpNetwork()
{

}

template<unsigned DIM>
VesselNetworkCellPopulationInteractor<DIM>::~VesselNetworkCellPopulationInteractor()
{

}

template<unsigned DIM>
void VesselNetworkCellPopulationInteractor<DIM>::LabelVesselsInCellPopulation(AbstractCellPopulation<DIM>& rCellPopulation,  QLength cellLengthScale,
                                  boost::shared_ptr<AbstractCellMutationState> pTipMutationState,
                                  boost::shared_ptr<AbstractCellMutationState> pStalkState,
                                  double threshold)
{
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();cell_iter != rCellPopulation.End();++cell_iter)
    {
        Vertex<DIM> cell_location = Vertex<DIM>(rCellPopulation.GetLocationOfCellCentre(*cell_iter), cellLengthScale);
        std::shared_ptr<VesselNode<DIM> > p_nearest_node = VesselNetworkGeometryCalculator<DIM>::GetNearestNode(mpNetwork, cell_location);
        QLength node_distance = p_nearest_node->GetDistance(cell_location);
        std::pair<std::shared_ptr<VesselSegment<DIM> >, QLength > segment_distance_pair =
                VesselNetworkGeometryCalculator<DIM>::GetNearestSegment(mpNetwork, cell_location);
        if (segment_distance_pair.second < threshold * cellLengthScale || node_distance < threshold * cellLengthScale)
        {
            if(p_nearest_node->IsMigrating())
            {
                cell_iter->SetMutationState(pTipMutationState);
            }
            else
            {
                cell_iter->SetMutationState(pStalkState);
            }
        }
    }
}

template<unsigned DIM>
void VesselNetworkCellPopulationInteractor<DIM>::PartitionNetworkOverCells(AbstractCellPopulation<DIM>& rCellPopulation,  QLength cellLengthScale,
                                                                           double threshold)
{
    // Loop through each cell add get its location. If the centre falls along a vessel, split the vessel at that location

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();cell_iter != rCellPopulation.End();++cell_iter)
    {
        Vertex<DIM> cell_location = Vertex<DIM>(rCellPopulation.GetLocationOfCellCentre(*cell_iter), cellLengthScale);
        std::shared_ptr<VesselNode<DIM> > p_nearest_node = VesselNetworkGeometryCalculator<DIM>::GetNearestNode(mpNetwork, cell_location);
        QLength node_distance = p_nearest_node->GetDistance(cell_location);

        std::pair<std::shared_ptr<VesselSegment<DIM> >, QLength > segment_distance_pair =
                VesselNetworkGeometryCalculator<DIM>::GetNearestSegment(mpNetwork, cell_location);
        if (segment_distance_pair.second < threshold*cellLengthScale || node_distance < threshold*cellLengthScale)
        {
            if (node_distance >= threshold*cellLengthScale)
            {
                std::shared_ptr<Vessel<DIM> > pVessel = segment_distance_pair.first->GetVessel();
                pVessel->DivideSegment(cell_location);
                pVessel->UpdateNodes();
            }
            mpNetwork->UpdateNodes();
            mpNetwork->UpdateSegments();
            mpNetwork->UpdateVesselNodes();
        }
    }
}

template<unsigned DIM>
void VesselNetworkCellPopulationInteractor<DIM>::KillNonVesselOverlappingCells(AbstractCellPopulation<DIM>& rCellPopulation,  QLength cellLengthScale,
                                                                               double threshold)
{
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();cell_iter != rCellPopulation.End();++cell_iter)
    {
        Vertex<DIM> cell_location = Vertex<DIM>(rCellPopulation.GetLocationOfCellCentre(*cell_iter), cellLengthScale);
        std::shared_ptr<VesselNode<DIM> > p_nearest_node = VesselNetworkGeometryCalculator<DIM>::GetNearestNode(mpNetwork, cell_location);
        QLength node_distance = p_nearest_node->GetDistance(cell_location);

        std::pair<std::shared_ptr<VesselSegment<DIM> >, QLength > segment_distance_pair = VesselNetworkGeometryCalculator<DIM>::GetNearestSegment(mpNetwork, cell_location);
        if (segment_distance_pair.second > threshold*cellLengthScale and node_distance > threshold*cellLengthScale)
        {
            cell_iter->Kill();
        }
    }
    rCellPopulation.RemoveDeadCells();
}

template<unsigned DIM>
void VesselNetworkCellPopulationInteractor<DIM>::KillOverlappingVesselCells(AbstractCellPopulation<DIM>& rCellPopulation,  QLength cellLengthScale,
                                                                            double threshold)
{
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();cell_iter != rCellPopulation.End();++cell_iter)
    {
        Vertex<DIM> cell_location = Vertex<DIM>(rCellPopulation.GetLocationOfCellCentre(*cell_iter), cellLengthScale);
        std::shared_ptr<VesselNode<DIM> > p_nearest_node = VesselNetworkGeometryCalculator<DIM>::GetNearestNode(mpNetwork, cell_location);
        QLength node_distance = p_nearest_node->GetDistance(cell_location);

        std::pair<std::shared_ptr<VesselSegment<DIM> >, QLength > segment_distance_pair = VesselNetworkGeometryCalculator<DIM>::GetNearestSegment(mpNetwork, cell_location);
        if (segment_distance_pair.second < threshold*cellLengthScale or node_distance < threshold*cellLengthScale)
        {
            cell_iter->Kill();
        }
    }
    rCellPopulation.RemoveDeadCells();
}

template<unsigned DIM>
void VesselNetworkCellPopulationInteractor<DIM>::SetVesselNetwork(std::shared_ptr<VesselNetwork<DIM> > pNetwork)
{
    mpNetwork = pNetwork;
}

// Explicit instantiation
template class VesselNetworkCellPopulationInteractor<2> ;
template class VesselNetworkCellPopulationInteractor<3> ;
