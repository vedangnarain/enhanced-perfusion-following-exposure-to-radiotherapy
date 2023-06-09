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

#ifndef IMAGETOSURFACE_HPP_
#define IMAGETOSURFACE_HPP_
#include "SmartPointers.hpp"
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

/**
* Extract a vtkpolydata surface from an image using thresholding.
* If marching cubes are used the result is an isosurface on the 'threshold' value. Otherwise
* the surface is composed of all regions either above or below the threshold value. Surfaces
* may need 'cleaning' before further processing. This can be done with a 'SurfaceCleaner'.
*/
class ImageToSurface
{
    /**
     * The image
     */
    vtkSmartPointer<vtkImageData> mpImage;

    /**
     * Whether to segment above a certain threshold
     */
    bool mSegmentAboveThreshold;

    /**
     * The segmentation threshold
     */
    double mThreshold;

    /**
     * The output surface
     */
    vtkSmartPointer<vtkPolyData> mpSurface;

    /**
     * Whether to use marching cubes
     */
    bool mUseMarchingCubes;

    /**
     * Whether to remove disconnected small regions
     */
    bool mRemoveDisconnected;

public:

    /**
     * Constructor
     */
    ImageToSurface();

    /**
     * Destructor
     */
    ~ImageToSurface();

    /**
     * Factory constructor method
     * @return a pointer to the converter
     */
    static std::shared_ptr<ImageToSurface> Create();

    /**
     * Return the filter output
     * @return the filter output
     */
    vtkSmartPointer<vtkPolyData> GetOutput();

    /**
     * Set the image
     * @param pImage the input image
     */
    void SetInput(vtkSmartPointer<vtkImageData> pImage);

    /**
     * Set the image
     * @param pImage the input image
     */
    void SetInputRaw(vtkImageData* pImage);

    /**
     * Set the threshold value for segmentation
     * @param threshold the threshold value for segmentation
     * @param segmentAboveThreshold whether to segment above the threshold
     */
    void SetThreshold(double threshold, bool segmentAboveThreshold);

    /**
     * Set whether to use marching cubes for the segmentation
     * @param useMarchingCubes whether to use marching cubes for the segmentation
     */
    void SetUseMarchingCubes(bool useMarchingCubes);

    /**
     * Only retain the largest region
     * @param removeDisconnected only retain the largest region
     */
    void SetRemoveDisconnected(bool removeDisconnected);

    /**
     * Run the update
     */
    void Update();

};

#endif /*ImageToSurface_HPP_*/
