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



#include "Owen11Parameters.hpp"
#include "MetabolicStimulusCalculator.hpp"

template<unsigned DIM>
MetabolicStimulusCalculator<DIM>::MetabolicStimulusCalculator() : AbstractVesselNetworkCalculator<DIM>(),
        mQRef(Owen11Parameters::mpReferenceFlowRateForMetabolicStimulus->GetValue("MetabolicStimulusCalculator")),
        mKm(Owen11Parameters::mpBasalMetabolicStimulus->GetValue("MetabolicStimulusCalculator")),
        mMaxStimulus(1.e12 * unit::per_second)
{

}

template<unsigned DIM>
MetabolicStimulusCalculator<DIM>::~MetabolicStimulusCalculator()
{

}

template <unsigned DIM>
std::shared_ptr<MetabolicStimulusCalculator<DIM> > MetabolicStimulusCalculator<DIM>::Create()
{
    return std::make_shared<MetabolicStimulusCalculator<DIM> >();

}


template<unsigned DIM>
QFlowRate MetabolicStimulusCalculator<DIM>::GetQRef()
{
    return mQRef;
}

template<unsigned DIM>
QRate MetabolicStimulusCalculator<DIM>::GetKm()
{
    return mKm;
}

template<unsigned DIM>
QRate MetabolicStimulusCalculator<DIM>::GetMaxStimulus()
{
    return mMaxStimulus;
}

template<unsigned DIM>
void MetabolicStimulusCalculator<DIM>::SetQRef(QFlowRate qRef)
{
    mQRef = Qabs(qRef);
}

template<unsigned DIM>
void MetabolicStimulusCalculator<DIM>::SetKm(QRate km)
{
    mKm = Qabs(km);
}

template<unsigned DIM>
void MetabolicStimulusCalculator<DIM>::SetMaxStimulus(QRate maxStimulus)
{
    mMaxStimulus = Qabs(maxStimulus);
}

template<unsigned DIM>
void MetabolicStimulusCalculator<DIM>::Calculate()
{
    std::vector<std::shared_ptr<VesselSegment<DIM> > > segments = this->mpNetwork->GetVesselSegments();
    for (unsigned idx = 0; idx < segments.size(); idx++)
    {
        QRate metabolic_stimulus;
        QDimensionless haematocrit = segments[idx]->GetFlowProperties()->GetHaematocrit();
        QFlowRate flow_rate = Qabs(segments[idx]->GetFlowProperties()->GetFlowRate());

        if (flow_rate > 0.0 * unit::metre_cubed_per_second)
        {
            if (haematocrit > 0.0)
            {
                metabolic_stimulus = mKm * log10(mQRef / (flow_rate * haematocrit) + 1.0);
            }
            else
            {
                metabolic_stimulus = mKm;
            }
        }
        else
        {
            metabolic_stimulus = 0.0 * unit::per_second;
        }
        segments[idx]->GetFlowProperties()->SetGrowthStimulus(segments[idx]->GetFlowProperties()->GetGrowthStimulus() + metabolic_stimulus);
    }
}

// Explicit instantiation
template class MetabolicStimulusCalculator<2>;
template class MetabolicStimulusCalculator<3>;
