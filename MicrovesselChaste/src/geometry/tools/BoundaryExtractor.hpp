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

#ifndef BOUNDARYEXTRACTOR_HPP_
#define BOUNDARYEXTRACTOR_HPP_

#include <memory>
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning
#include <vtkSmartPointer.h>

/**
 * Forward declare VTK members
 */
class vtkPolyData;

/**
* This class extracts exterior boundaries from 2D surfaces defined as vtk polydata. It includes methods
* for smoothing the resulting line boundary using spline resampling of the boundary.
*/
class BoundaryExtractor
{
    /**
     *  The input surface
     */
    vtkSmartPointer<vtkPolyData> mpInputSurface;

    /**
     *  The output surface
     */
    vtkSmartPointer<vtkPolyData> mpOutputSurface;

    /**
     *  A characteristic length for the spline resampling, a greater length leads to more smoothing.
     */
    double mSmoothingLength;

    /**
     *  Whether or not to do smoothing
     */
    bool mDoSmoothing;

    /**
     *  Whether to remove disconnected regions
     */
    bool mRemoveDisconnected;

public:

    /**
     * Constructor
     */
    BoundaryExtractor();

    /**
     * Destructor
     */
    ~BoundaryExtractor();

    /**
     * Factory constructor method
     * @return a pointer to a class instance
     */
    static std::shared_ptr<BoundaryExtractor> Create();

    /**
     * Get the boundary edges
     * @return the boundary edges
     */
    vtkSmartPointer<vtkPolyData> GetOutput();

    /**
     * A 2D surface from which boundary edges are to be extracted
     * @param pInputSurface the surface
     */
    void SetInput(vtkSmartPointer<vtkPolyData> pInputSurface);

    /**
     * A 2D surface from which boundary edges are to be extracted. Raw pointer
     * for Python interface.
     * @param pInputSurface the surface
     */
    void SetInputRaw(vtkPolyData* pInputSurface);

    /**
     * A length for the spline resampling filter
     * @param value the resampling length
     */
    void SetSmoothingLength(double value);

    /**
     * Whether to do the smoothing
     * @param doSmoothing whether to do the smoothing
     */
    void SetDoSmoothing(bool doSmoothing);

    /**
     * Only retain the largest region
     * @param removeDisconnected only retain the largest region
     */
    void SetRemoveDisconnected(bool removeDisconnected);

    /**
     * Run the tool
     */
    void Update();
};

#endif /*BoundaryExtractor_HPP_*/
