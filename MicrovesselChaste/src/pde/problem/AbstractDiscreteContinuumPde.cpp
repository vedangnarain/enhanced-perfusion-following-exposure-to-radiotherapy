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

#include <algorithm>
#include "AbstractDiscreteContinuumPde.hpp"
#include "BaseUnits.hpp"
#include "DiscreteSource.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractDiscreteContinuumPde<ELEMENT_DIM, SPACE_DIM>::AbstractDiscreteContinuumPde() :
            mDiffusivity(1.0 * unit::metre_squared_per_second),
            mDiscreteSources(),
            mReferenceConcentration(BaseUnits::Instance()->GetReferenceConcentrationScale()),
            mReferenceLengthScale(BaseUnits::Instance()->GetReferenceLengthScale()),
            mReferenceTimeScale(BaseUnits::Instance()->GetReferenceTimeScale())
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractDiscreteContinuumPde<ELEMENT_DIM, SPACE_DIM>::~AbstractDiscreteContinuumPde()
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractDiscreteContinuumPde<ELEMENT_DIM, SPACE_DIM>::AddDiscreteSource(std::shared_ptr<DiscreteSource<SPACE_DIM> > pDiscreteSource)
{
    mDiscreteSources.push_back(pDiscreteSource);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
QDiffusivity AbstractDiscreteContinuumPde<ELEMENT_DIM, SPACE_DIM>::ComputeIsotropicDiffusionTerm()
{
    return mDiffusivity;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<std::shared_ptr<DiscreteSource<SPACE_DIM> > > AbstractDiscreteContinuumPde<ELEMENT_DIM, SPACE_DIM>::GetDiscreteSources()
{
    return mDiscreteSources;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractDiscreteContinuumPde<ELEMENT_DIM, SPACE_DIM>::SetIsotropicDiffusionConstant(QDiffusivity diffusivity)
{
    mDiffusivity = diffusivity;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractDiscreteContinuumPde<ELEMENT_DIM, SPACE_DIM>::SetReferenceConcentration(QConcentration referenceConcentration)
{
    mReferenceConcentration = referenceConcentration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractDiscreteContinuumPde<ELEMENT_DIM, SPACE_DIM>::SetReferenceLength(QLength referenceLength)
{
    mReferenceLengthScale = referenceLength;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractDiscreteContinuumPde<ELEMENT_DIM, SPACE_DIM>::SetReferenceTime(QTime referenceTime)
{
    mReferenceTimeScale = referenceTime;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractDiscreteContinuumPde<ELEMENT_DIM, SPACE_DIM>::UpdateDiscreteSourceStrengths()
{

}

// Explicit instantiation
template class AbstractDiscreteContinuumPde<2>;
template class AbstractDiscreteContinuumPde<3>;
