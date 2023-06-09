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

#include "CellStateDependentDiscreteSource.hpp"
#include "AbstractCellPopulation.hpp"

template<unsigned DIM>
CellStateDependentDiscreteSource<DIM>::CellStateDependentDiscreteSource()
    :   DiscreteSource<DIM>(),
        mStateRateMap(),
        mStateRateThresholdMap(),
        mConsumptionRatePerUnitConcentration(1.0*unit::metre_cubed_per_second)
{

}

template<unsigned DIM>
CellStateDependentDiscreteSource<DIM>::~CellStateDependentDiscreteSource()
{

}

template<unsigned DIM>
std::shared_ptr<CellStateDependentDiscreteSource<DIM> > CellStateDependentDiscreteSource<DIM>::Create()
{
    return std::make_shared<CellStateDependentDiscreteSource<DIM> >();

}

template<unsigned DIM>
std::vector<QConcentrationFlowRate > CellStateDependentDiscreteSource<DIM>::GetConstantInUValues()
{
    if(!this->mpDensityMap->GetGridCalculator())
    {
        EXCEPTION("A regular grid is required for this type of source");
    }

    std::vector<QConcentrationFlowRate > values(this->mpDensityMap->GetGridCalculator()->GetGrid()->GetNumberOfCells(),
            0.0*unit::mole_per_metre_cubed_per_second);

    std::shared_ptr<ApoptoticCellProperty> apoptotic_property(new ApoptoticCellProperty);
    unsigned apoptotic_label = apoptotic_property->GetColour();

    // Loop through all points
    std::vector<std::vector<CellPtr> > point_cell_map = this->mpDensityMap->GetGridCalculator()->rGetCellMap();
    for(unsigned idx=0; idx<point_cell_map.size(); idx++)
    {
        for(unsigned jdx=0; jdx<point_cell_map[idx].size(); jdx++)
        {
            // If a mutation specific consumption rate has been specified
            if(mStateRateMap.size()>0)
            {
                std::map<unsigned, QConcentrationFlowRate >::iterator it;
                // If the cell is apoptotic
                if (point_cell_map[idx][jdx]->template HasCellProperty<ApoptoticCellProperty>())
                {
                    it = mStateRateMap.find(apoptotic_label);
                    if (it != mStateRateMap.end())
                    {
                        values[idx] = values[idx] + it->second;
                    }
                }
                else
                {
                    it = mStateRateMap.find(point_cell_map[idx][jdx]->GetMutationState()->GetColour());
                    if (it != mStateRateMap.end())
                    {
                        if(!this->mLabel.empty())
                        {
                            // Get a threshold value if it has been set, use the label to determine the field from which the label
                            // value is obtained.
                            QConcentration threshold = 0.0 * unit::mole_per_metre_cubed;
                            if(mStateRateThresholdMap.size()>0)
                            {
                                std::map<unsigned, QConcentration >::iterator it_threshold;
                                it_threshold = mStateRateThresholdMap.find(point_cell_map[idx][jdx]->GetMutationState()->GetColour());
                                if (it_threshold != mStateRateThresholdMap.end())
                                {
                                    threshold = it_threshold->second;
                                }
                            }
//                            QConcentrationFlowRate conversion_factor = mConsumptionRatePerUnitConcentration*unit::mole_per_metre_cubed/grid_volume;

                            if(threshold>0.0* unit::mole_per_metre_cubed)
                            {
                                if(point_cell_map[idx][jdx]->GetCellData()->GetItem(this->mLabel)>threshold.GetValue())
                                {
                                    values[idx] = values[idx] + it->second;
                                }
                            }
                            else
                            {
                                values[idx] = values[idx] + it->second;
                            }
                        }
                        else
                        {
                            values[idx] = values[idx] + it->second;
                        }
                    }
                }
            }
            else
            {
                values[idx] = values[idx] + this->mConstantInUValue;
            }
        }
    }
    return values;
}

template<unsigned DIM>
std::vector<QRate > CellStateDependentDiscreteSource<DIM>::GetLinearInUValues()
{
    return std::vector<QRate >(this->mpDensityMap->GetGridCalculator()->GetGrid()->GetNumberOfCells(), 0.0*unit::per_second);
}

template<unsigned DIM>
void CellStateDependentDiscreteSource<DIM>::SetStateRateMap(std::map<unsigned, QConcentrationFlowRate > stateRateMap)
{
    mStateRateMap = stateRateMap;
}

template<unsigned DIM>
void CellStateDependentDiscreteSource<DIM>::SetStateRateThresholdMap(std::map<unsigned, QConcentration > stateThresholdMap)
{
    mStateRateThresholdMap = stateThresholdMap;
}

// Explicit instantiation
template class CellStateDependentDiscreteSource<2>;
template class CellStateDependentDiscreteSource<3>;
