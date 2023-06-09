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

#include "GeometryTools.hpp"
#include "LatticeBasedMigrationRule.hpp"
#include "RandomNumberGenerator.hpp"
#include "RegularGrid.hpp"

template<unsigned DIM>
LatticeBasedMigrationRule<DIM>::LatticeBasedMigrationRule()
    : AbstractMigrationRule<DIM>(),
      mMovementProbability(0.01)
{

}

template <unsigned DIM>
std::shared_ptr<LatticeBasedMigrationRule<DIM> > LatticeBasedMigrationRule<DIM>::Create()
{
    return std::make_shared<LatticeBasedMigrationRule<DIM> >();
}

template<unsigned DIM>
LatticeBasedMigrationRule<DIM>::~LatticeBasedMigrationRule()
{

}

template<unsigned DIM>
std::vector<double> LatticeBasedMigrationRule<DIM>::GetNeighbourMovementProbabilities(std::shared_ptr<VesselNode<DIM> > pNode,
                                                       std::vector<unsigned> neighbourIndices, unsigned gridIndex)
{
    std::vector<double> probability_of_moving(neighbourIndices.size(), 0.0);
    for(unsigned idx=0; idx<neighbourIndices.size(); idx++)
    {
        // Make sure that tip cell does not try to move into a location already occupied by the vessel that it comes from
        Vertex<DIM> neighbour_location = this->mpGridCalculator->GetGrid()->GetGlobalCellLocation(neighbourIndices[idx]);

        bool already_attached = false;
        for (unsigned seg_index = 0; seg_index < pNode->GetNumberOfSegments(); seg_index++)
        {
            if(pNode->GetSegment(seg_index)->GetOppositeNode(pNode)->IsCoincident(neighbour_location))
            {
                already_attached = true;
                break;
            }
        }

        if(already_attached)
        {
            continue;
        }

        // Simple rule, equal probability for all directions
        probability_of_moving[idx] = mMovementProbability * SimulationTime::Instance()->GetTimeStep();
    }
    return probability_of_moving;
}

template<unsigned DIM>
int LatticeBasedMigrationRule<DIM>::GetNeighbourMovementIndex(std::vector<double> movementProbabilities,
                                                                   std::vector<unsigned> neighbourIndices)
{
    int location_index = -1;

    // Check that the cumulative movement probability is less than one, otherwise our time-step is too large
    std::vector<double> cumulativeProbabilityVector(movementProbabilities.size());
    std::partial_sum(movementProbabilities.begin(), movementProbabilities.end(), cumulativeProbabilityVector.begin());

    if (cumulativeProbabilityVector.back() > 1.0)
    {
        EXCEPTION("Cumulative probability of tip cell moving is greater than one");
    }

    // Use roulette-wheel style selection to select which location the tip will move into
    double cumulativeProbability = cumulativeProbabilityVector.back();
    double random_number = RandomNumberGenerator::Instance()->ranf();

    // If we move, choose a node to go to
    if(random_number < cumulativeProbability)
    {
        for (unsigned ind = 0; ind < cumulativeProbabilityVector.size(); ind++)
        {
            if (random_number <= cumulativeProbabilityVector[ind])
            {
                location_index = neighbourIndices[ind];
                break;
            }
        }
    }
    return location_index;
}

template<unsigned DIM>
void LatticeBasedMigrationRule<DIM>::SetMovementProbability(double movementProbability)
{
    mMovementProbability = movementProbability;
}

template<unsigned DIM>
std::vector<int> LatticeBasedMigrationRule<DIM>::GetIndices(const std::vector<std::shared_ptr<VesselNode<DIM> > >& rNodes)
{
    if(!this->mpGridCalculator)
    {
        EXCEPTION("A regular grid is required for this type of migration rule.");
    }

    if(!this->mpVesselNetwork)
    {
        EXCEPTION("A vessel network is required for this type of migration rule.");
    }

    // Set up the output indices vector
    std::vector<int> indices(rNodes.size(), -1);

    // Get the point-node map from the regular grid
    std::vector<std::vector<std::shared_ptr<VesselNode<DIM> > > > point_node_map = this->mpGridCalculator->rGetVesselNodeMap();

    // Need a regular grid for this rule
    std::shared_ptr<RegularGrid<DIM> > p_regular_grid =
            std::dynamic_pointer_cast<RegularGrid<DIM> >(this->mpGridCalculator->GetGrid());

    if(!p_regular_grid)
    {
        EXCEPTION("Can't cast to regular grid");
    }

    // Get the neighbour data from the regular grid
    std::vector<std::vector<unsigned> > neighbour_indices;
    if(this->mUseMooreNeighbourhood)
    {
        neighbour_indices = p_regular_grid->rGetMooreNeighbourData();
    }
    else
    {
        neighbour_indices = p_regular_grid->rGetNeighbourData();
    }

    // Loop over all nodes, if they can move set the index
    for(unsigned idx = 0; idx < rNodes.size(); idx++)
    {
        // Get the grid index of the node
        unsigned grid_index = p_regular_grid->GetNearestCellIndex(rNodes[idx]->rGetLocation());

        // Get the probability of moving into each of the neighbour sites
        std::vector<double> probability_of_moving = GetNeighbourMovementProbabilities(rNodes[idx], neighbour_indices[grid_index], grid_index);

        // Get the index of the neighbour to move into
        double sum = std::fabs(std::accumulate(probability_of_moving.begin(), probability_of_moving.end(), 0.0));
        if(sum > 0.0)
        {
            indices[idx] = GetNeighbourMovementIndex(probability_of_moving, neighbour_indices[grid_index]);
        }
    }
    return indices;
}

// Explicit instantiation
template class LatticeBasedMigrationRule<2> ;
template class LatticeBasedMigrationRule<3> ;
