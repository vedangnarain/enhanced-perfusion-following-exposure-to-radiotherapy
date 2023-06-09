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

#include "AbstractVesselNetworkComponent.hpp"
#include "Exception.hpp"

template<unsigned DIM>
AbstractVesselNetworkComponent<DIM>::AbstractVesselNetworkComponent() :
        mOutputData(),
        mId(0),
        mRadius(10.0*unit::microns)
{

}

template<unsigned DIM>
AbstractVesselNetworkComponent<DIM>::~AbstractVesselNetworkComponent()
{
}

template<unsigned DIM>
unsigned AbstractVesselNetworkComponent<DIM>::GetId() const
{
    return mId;
}

template<unsigned DIM>
double AbstractVesselNetworkComponent<DIM>::GetOutputDataValue(const std::string& rKey)
{
    // Make sure the output data is up to date
    GetOutputData();

    std::map<std::string,double>::const_iterator it = mOutputData.find(rKey);
    if (it != mOutputData.end())
    {
        return it->second;
    }
    else
    {
        EXCEPTION("Requested output data key not found");
    }
}

template<unsigned DIM>
std::vector<std::string> AbstractVesselNetworkComponent<DIM>::GetOutputDataKeys()
{
    // Make sure the output data is up to date
    GetOutputData();

    std::vector<std::string> keys;
    for(std::map<std::string,double>::const_iterator it = mOutputData.begin(); it != mOutputData.end(); ++it)
    {
        keys.push_back(it->first);
    }
    return keys;
}

template<unsigned DIM>
QLength AbstractVesselNetworkComponent<DIM>::GetRadius() const
{
    return mRadius;
}

template<unsigned DIM>
void AbstractVesselNetworkComponent<DIM>::SetOutputData(const std::string& rKey, double value)
{
    mOutputData[rKey] = value;
}

template<unsigned DIM>
void AbstractVesselNetworkComponent<DIM>::SetId(unsigned id)
{
    mId = id;
}

template<unsigned DIM>
void AbstractVesselNetworkComponent<DIM>::SetRadius(QLength radius)
{
    mRadius = radius;
}

// Explicit instantiation
template class AbstractVesselNetworkComponent<2> ;
template class AbstractVesselNetworkComponent<3> ;
