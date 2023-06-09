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



#ifndef CSVVESSELNETWORKREADER_HPP_
#define CSVVESSELNETWORKREADER_HPP_

#include <string>
#include "SmartPointers.hpp"
#include "VesselNetwork.hpp"

/**
 * This class reads vessel network skeleton information in the format used by the ImageJ Skeletonize3D/Analyze Skeleton
 * plugins and creates a vessel network.
 */
template<unsigned DIM>
class CsvVesselNetworkReader
{
    /**
     * The input file name
     */
    std::string mFileName;

    /**
     * The label for radii
     */
    std::string mRadiusLabel;

    QLength mReferenceLength;

public:

    /**
     * Constructor
     */
    CsvVesselNetworkReader();

    /**
     * Destructor
     */
    ~CsvVesselNetworkReader();

    /**
     * Construct a new instance of the class and return a shared pointer to it.
     * @return a new reader
     */
    static std::shared_ptr<CsvVesselNetworkReader<DIM> > Create();

    /**
     * Set the input file name
     * @param rFileName the input file name
     */
    void SetFileName(const std::string& rFileName);

    /**
     * Do the read and return the vessel network
     * @return the vessel network
     */
    std::shared_ptr<VesselNetwork<DIM> > Read();

};

#endif /* CSVVESSELNETWORKREADER_HPP_ */
