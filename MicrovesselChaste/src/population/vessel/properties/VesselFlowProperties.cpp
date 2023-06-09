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

#include "VesselFlowProperties.hpp"

template<unsigned DIM>
VesselFlowProperties<DIM>::VesselFlowProperties() : AbstractVesselNetworkComponentFlowProperties<DIM>(),
    mSegments(),
    mUndergoingRegression(false),
    mRemoveViaRegression(false),
    mRegressionTime(DBL_MAX*unit::seconds)
{
}

template<unsigned DIM>
VesselFlowProperties<DIM>::~VesselFlowProperties()
{
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::CheckSegments() const
{
    if(mSegments.size() == 0)
    {
        EXCEPTION("No vessel segments have been set for vessel flow property calculation.");
    }
}

template<unsigned DIM>
QConcentration VesselFlowProperties<DIM>::GetAntiAngiogenicDrugConcentration() const
{
    this->CheckSegments();

    QConcentration value = 0.0*unit::mole_per_metre_cubed;
    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        value = value +  mSegments[i]->GetFlowProperties()->GetAntiAngiogenicDrugConcentration() / double(mSegments.size());
    }
    return value;
}

template<unsigned DIM>
QDimensionless VesselFlowProperties<DIM>::GetHaematocrit() const
{
    this->CheckSegments();

    QDimensionless value = 0.0;
    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        value = value +  mSegments[i]->GetFlowProperties()->GetHaematocrit() / double(mSegments.size());
    }
    return value;
}

template<unsigned DIM>
QFlowRate VesselFlowProperties<DIM>::GetFlowRate() const
{
    this->CheckSegments();

    QFlowRate value = 0.0 * unit::metre_cubed_per_second;
    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        value = value +  mSegments[i]->GetFlowProperties()->GetFlowRate() / double(mSegments.size());
    }
    return value;
}

template<unsigned DIM>
QFlowImpedance VesselFlowProperties<DIM>::GetImpedance() const
{
    this->CheckSegments();

    QFlowImpedance value = 0.0 * unit::pascal_second_per_metre_cubed;
    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        value = value +  mSegments[i]->GetFlowProperties()->GetImpedance();
    }
    return value;
}

template<unsigned DIM>
QDynamicViscosity VesselFlowProperties<DIM>::GetViscosity() const
{
    this->CheckSegments();

    QDynamicViscosity value = 0.0 * unit::poiseuille;
    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        value = value +  mSegments[i]->GetFlowProperties()->GetViscosity()/ double(mSegments.size());
    }
    return value;
}

template<unsigned DIM>
QPressure VesselFlowProperties<DIM>::GetWallShearStress() const
{
    this->CheckSegments();

    QPressure value = 0.0 * unit::pascals;
    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        value = value +  mSegments[i]->GetFlowProperties()->GetWallShearStress()/ double(mSegments.size());
    }
    return value;
}

template<unsigned DIM>
QRate VesselFlowProperties<DIM>::GetGrowthStimulus() const
{
    this->CheckSegments();

    QRate value = 0.0 * unit::per_second;
    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        value = value + mSegments[i]->GetFlowProperties()->GetGrowthStimulus()/ double(mSegments.size());
    }
    return value;
}

template<unsigned DIM>
std::map<std::string, double> VesselFlowProperties<DIM>::GetOutputData() const
{
    std::map<std::string, double> output_data;
    output_data["Vessel Impedance kg/m^4/s"] = this->GetImpedance() / unit::pascal_second_per_metre_cubed;
    output_data["Vessel Haematocrit"] = this->GetHaematocrit();
    output_data["Vessel Flow Rate m^3/s"] = this->GetFlowRate() / unit::metre_cubed_per_second;
    output_data["Absolute Vessel Flow Rate m^3/s"] = fabs(this->GetFlowRate()) / unit::metre_cubed_per_second;
    output_data["Vessel Viscosity Pa.s"] = this->GetViscosity() / unit::poiseuille;
    output_data["Vessel Wall Shear Stress Pa"] = this->GetWallShearStress()  / unit::pascals;
    output_data["Vessel Growth Stimulus s^-1"] = this->GetGrowthStimulus() / unit::per_second;
    output_data["Vessel Time Until Regression s"] = this->mRegressionTime / unit::seconds;
    return output_data;
}

template<unsigned DIM>
QTime VesselFlowProperties<DIM>::GetRegressionTime() const
{
    return mRegressionTime;
}

template<unsigned DIM>
bool VesselFlowProperties<DIM>::HasVesselRegressed(QTime simulationReferenceTime)
{
    if (SimulationTime::Instance()->GetTime()*simulationReferenceTime >= this->mRegressionTime)
    {
        assert(mUndergoingRegression);
        mRemoveViaRegression = true;
        return mRemoveViaRegression;
    }
    else
    {
        mRemoveViaRegression = false;
        return mRemoveViaRegression;
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetHaematocrit(QDimensionless haematocrit)
{
    this->CheckSegments();

    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        mSegments[i]->GetFlowProperties()->SetHaematocrit(haematocrit);
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetFlowRate(QFlowRate haematocrit)
{
    this->CheckSegments();

    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        mSegments[i]->GetFlowProperties()->SetFlowRate(haematocrit);
    }
}


template<unsigned DIM>
void VesselFlowProperties<DIM>::SetImpedance(QFlowImpedance value)
{
    this->CheckSegments();

    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        mSegments[i]->GetFlowProperties()->SetImpedance(value/double(mSegments.size()));
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetViscosity(QDynamicViscosity value)
{
    this->CheckSegments();

    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        mSegments[i]->GetFlowProperties()->SetViscosity(value);
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetWallShearStress(QPressure value)
{
    this->CheckSegments();

    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        mSegments[i]->GetFlowProperties()->SetWallShearStress(value);
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetGrowthStimulus(QRate value)
{
    this->CheckSegments();

    for (unsigned i = 0; i < mSegments.size(); i++)
    {
        mSegments[i]->GetFlowProperties()->SetGrowthStimulus(value);
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetTimeUntilRegression(QTime time, QTime simulationReferenceTime)
{
    assert(!mRemoveViaRegression);

    if (HasRegressionTimerStarted())
    {
        EXCEPTION("SetTimeUntilRegression(time) called when already undergoing regression");
    }

    mUndergoingRegression = true;
    this->mRegressionTime = SimulationTime::Instance()->GetTime()*simulationReferenceTime + time;
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetRegressionTime(QTime time)
{
    this->mRegressionTime = time;
    if(time<DBL_MAX*unit::seconds)
    {
        mUndergoingRegression = true;
    }
}

template<unsigned DIM>
bool VesselFlowProperties<DIM>::HasRegressionTimerStarted()
{
    return mUndergoingRegression;
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::ResetRegressionTimer()
{
    this->mRegressionTime = DBL_MAX*unit::seconds;
    mUndergoingRegression = false;
    mRemoveViaRegression = false;
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::UpdateSegments(std::vector<std::shared_ptr<VesselSegment<DIM> > > segments)
{
    mSegments = segments;
}

// Explicit instantiation
template class VesselFlowProperties<2> ;
template class VesselFlowProperties<3> ;
