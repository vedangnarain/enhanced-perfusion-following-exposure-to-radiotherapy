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

#include "AbstractMigrationRule.hpp"

template<unsigned DIM>
AbstractMigrationRule<DIM>::AbstractMigrationRule()
     :mpSolver(),
      mpVesselNetwork(),
      mIsSprouting(false),
      mpCellPopulation(),
      mpGridCalculator(),
      mpBoundingDomain(),
      mUseMooreNeighbourhood(false)
{

}

template<unsigned DIM>
AbstractMigrationRule<DIM>::~AbstractMigrationRule()
{

}

template <unsigned DIM>
std::shared_ptr<AbstractMigrationRule<DIM> > AbstractMigrationRule<DIM>::Create()
{
    return std::make_shared<AbstractMigrationRule<DIM> >();

}

template <unsigned DIM>
void AbstractMigrationRule<DIM>::SetUseMooreNeighbourhood(bool useMooreNeighbourhood)
{
    mUseMooreNeighbourhood = useMooreNeighbourhood;
}

template <unsigned DIM>
std::vector<Vertex<DIM> > AbstractMigrationRule<DIM>::GetDirections(const std::vector<VesselNodePtr<DIM> >& rNodes)
{
    return std::vector<Vertex<DIM> >();
}

template <unsigned DIM>
std::vector<int> AbstractMigrationRule<DIM>::GetIndices(const std::vector<VesselNodePtr<DIM> >& rNodes)
{
    return std::vector<int>();
}

template<unsigned DIM>
void AbstractMigrationRule<DIM>::SetBoundingDomain(PartPtr<DIM> pPart)
{
    mpBoundingDomain = pPart;
}

template<unsigned DIM>
void AbstractMigrationRule<DIM>::SetIsSprouting(bool isSprouting)
{
    mIsSprouting = isSprouting;
}

template<unsigned DIM>
void AbstractMigrationRule<DIM>::SetGridCalculator(std::shared_ptr<GridCalculator<DIM> > pGrid)
{
    mpGridCalculator = pGrid;
}

template<unsigned DIM>
void AbstractMigrationRule<DIM>::SetDiscreteContinuumSolver(std::shared_ptr<AbstractDiscreteContinuumSolver<DIM> > pSolver)
{
    mpSolver = pSolver;
}

template<unsigned DIM>
void AbstractMigrationRule<DIM>::SetNetwork(VesselNetworkPtr<DIM> pNetwork)
{
    mpVesselNetwork = pNetwork;
}

template<unsigned DIM>
void AbstractMigrationRule<DIM>::SetCellPopulation(std::shared_ptr<AbstractCellPopulation<DIM> > pPopulation)
{
    mpCellPopulation = pPopulation;
}

// Explicit instantiation
template class AbstractMigrationRule<2>;
template class AbstractMigrationRule<3>;
