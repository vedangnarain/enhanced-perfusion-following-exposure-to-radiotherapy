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



#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkUnsignedCharArray.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>
#include <vtkGlyph2D.h>
#include <vtkCubeAxesActor2D.h>
#include <vtkImageData.h>
#include <vtkGeometryFilter.h>
#include <vtkTubeFilter.h>
#include <vtkExtractEdges.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vtkPolygon.h>
#include <vtkIdList.h>
#include <vtkFeatureEdges.h>
#include "UblasIncludes.hpp"
#include "UblasVectorInclude.hpp"
#include "Exception.hpp"
#include "BaseUnits.hpp"
#include "VesselNetworkWriter.hpp"

#include "PartActorGenerator.hpp"


template<unsigned DIM>
PartActorGenerator<DIM>::PartActorGenerator()
    : AbstractActorGenerator<DIM>(),
      mpPart()
{

}

template<unsigned DIM>
PartActorGenerator<DIM>::~PartActorGenerator()
{

}

template<unsigned DIM>
void PartActorGenerator<DIM>::AddActor(vtkSmartPointer<vtkRenderer> pRenderer)
{
    if(mpPart)
    {
        vtkSmartPointer<vtkPolyData> p_polydata = mpPart->GetVtk();

        // Add the volume
        if(this->mShowVolume)
        {
            vtkSmartPointer<vtkPolyDataMapper> p_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            #if VTK_MAJOR_VERSION <= 5
                p_mapper->SetInput(p_polydata);
            #else
                p_mapper->SetInputData(p_polydata);
            #endif

            vtkSmartPointer<vtkActor> p_actor = vtkSmartPointer<vtkActor>::New();
            p_actor->GetProperty()->SetColor(this->mVolumeColor[0]/255.0, this->mVolumeColor[1]/255.0, this->mVolumeColor[2]/255.0);
            p_actor->GetProperty()->SetOpacity(this->mVolumeOpacity);
            p_actor->SetMapper(p_mapper);
            pRenderer->AddActor(p_actor);
        }

        // Add the points
        if(this->mShowPoints)
        {
            vtkSmartPointer<vtkSphereSource> p_spheres = vtkSmartPointer<vtkSphereSource>::New();
            p_spheres->SetRadius(this->mPointSize);
            p_spheres->SetPhiResolution(16);
            p_spheres->SetThetaResolution(16);

            vtkSmartPointer<vtkGlyph3D> p_glyph = vtkSmartPointer<vtkGlyph3D>::New();
            #if VTK_MAJOR_VERSION <= 5
                p_glyph->SetInput(p_polydata);
                p_glyph->SetSource(p_spheres->GetOutput());
            #else
                p_glyph->SetInputData(p_polydata);
                p_glyph->SetSourceConnection(p_spheres->GetOutputPort());
            #endif
            p_glyph->ClampingOff();
            p_glyph->SetScaleModeToScaleByScalar();
            p_glyph->SetScaleFactor(1.0);
            p_glyph->Update();

            vtkSmartPointer<vtkPolyDataMapper> p_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            #if VTK_MAJOR_VERSION <= 5
                p_mapper->SetInput(p_glyph->GetOutput());
            #else
                p_mapper->SetInputData(p_glyph->GetOutput());
            #endif
            p_mapper->ScalarVisibilityOn();

            vtkSmartPointer<vtkActor> p_actor = vtkSmartPointer<vtkActor>::New();
            p_actor->SetMapper(p_mapper);
            p_actor->GetProperty()->SetColor(this->mPointColor[0]/255.0, this->mPointColor[1]/255.0, this->mPointColor[2]/255.0);
            pRenderer->AddActor(p_actor);
        }

        // Add the edges
        if(this->mShowEdges)
        {
            vtkSmartPointer<vtkFeatureEdges> p_edges = vtkSmartPointer<vtkFeatureEdges>::New();
            #if VTK_MAJOR_VERSION <= 5
                p_edges->SetInput(p_polydata);
            #else
                p_edges->SetInputData(p_polydata);
            #endif
            p_edges->SetFeatureEdges(false);
            p_edges->SetBoundaryEdges(true);
            p_edges->SetManifoldEdges(true);
            p_edges->SetNonManifoldEdges(false);

            vtkSmartPointer<vtkTubeFilter> p_voronoi_tubes = vtkSmartPointer<vtkTubeFilter>::New();
            p_voronoi_tubes->SetInputConnection(p_edges->GetOutputPort());
            p_voronoi_tubes->SetRadius(this->mEdgeSize);
            p_voronoi_tubes->SetNumberOfSides(12);

            vtkSmartPointer<vtkPolyDataMapper> p_voronoi_tube_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            p_voronoi_tube_mapper->SetInputConnection(p_voronoi_tubes->GetOutputPort());

            vtkSmartPointer<vtkActor> p_voronoi_tube_actor = vtkSmartPointer<vtkActor>::New();
            p_voronoi_tube_actor->SetMapper(p_voronoi_tube_mapper);
            p_voronoi_tube_actor->GetProperty()->SetColor(this->mEdgeColor[0]/255.0, this->mEdgeColor[1]/255.0, this->mEdgeColor[2]/255.0);
            pRenderer->AddActor(p_voronoi_tube_actor);
        }

    }
}

template<unsigned DIM>
void PartActorGenerator<DIM>::SetPart(PartPtr<DIM> pPart)
{
    this->mpPart = pPart;
}

template class PartActorGenerator<2>;
template class PartActorGenerator<3>;
