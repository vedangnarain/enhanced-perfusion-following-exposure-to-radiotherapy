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

#include "SolutionDependentDiscreteSource.hpp"

template<unsigned DIM>
SolutionDependentDiscreteSource<DIM>::SolutionDependentDiscreteSource()
    :   DiscreteSource<DIM>(),
        mpSolution(),
        mConstantInUSinkRatePerSolutionQuantity(0.0*unit::per_second),
        mLinearInUSinkRatePerSolutionQuantity(0.0*unit::metre_cubed_per_mole_per_second)
{

}

template<unsigned DIM>
SolutionDependentDiscreteSource<DIM>::~SolutionDependentDiscreteSource()
{

}

template<unsigned DIM>
std::shared_ptr<SolutionDependentDiscreteSource<DIM> > SolutionDependentDiscreteSource<DIM>::Create()
{
    return std::make_shared<SolutionDependentDiscreteSource<DIM> >();

}

template<unsigned DIM>
std::vector<QConcentrationFlowRate > SolutionDependentDiscreteSource<DIM>::GetConstantInUValues()
{
    if(!this->mpDensityMap->GetGridCalculator())
    {
        EXCEPTION("A regular grid is required for this type of source");
    }

    unsigned num_points = this->mpDensityMap->GetGridCalculator()->GetGrid()->GetNumberOfCells();
    std::vector<QConcentrationFlowRate > values(num_points, 0.0*unit::mole_per_metre_cubed_per_second);
    if(mpSolution.size() != num_points)
    {
        EXCEPTION("A solution sampled on the grid is required for this type of source");
    }
    for(unsigned idx=0; idx<num_points; idx++)
    {
        values[idx] = mpSolution[idx]*mConstantInUSinkRatePerSolutionQuantity;
    }

    return values;
}

template<unsigned DIM>
std::vector<QRate > SolutionDependentDiscreteSource<DIM>::GetLinearInUValues()
{
    if(!this->mpDensityMap->GetGridCalculator())
    {
        EXCEPTION("A regular grid is required for this type of source");
    }

    unsigned num_points = this->mpDensityMap->GetGridCalculator()->GetGrid()->GetNumberOfCells();
    std::vector<QRate > values(num_points, 0.0*unit::per_second);
    if(mpSolution.size() != num_points)
    {
        EXCEPTION("A solution sampled on the grid is required for this type of source");
    }
    for(unsigned idx=0; idx<num_points; idx++)
    {
        values[idx] = mpSolution[idx]*mLinearInUSinkRatePerSolutionQuantity;
    }

    return values;
}

template<unsigned DIM>
void SolutionDependentDiscreteSource<DIM>::SetSolution(std::vector<QConcentration > solution)
{
    mpSolution = solution;
}

template<unsigned DIM>
void SolutionDependentDiscreteSource<DIM>::SetConstantInUSinkRatePerSolutionQuantity(QRate value)
{
    mConstantInUSinkRatePerSolutionQuantity = value;
}

template<unsigned DIM>
void SolutionDependentDiscreteSource<DIM>::SetLinearInUSinkRatePerSolutionQuantity(QRatePerConcentration value)
{
    mLinearInUSinkRatePerSolutionQuantity = value;
}

// Explicit instantiation
template class SolutionDependentDiscreteSource<2>;
template class SolutionDependentDiscreteSource<3>;
