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

#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning for now (gcc4.3)
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkProbeFilter.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include "AbstractDiscreteContinuumPde.hpp"
#include "AbstractUnstructuredGridDiscreteContinuumSolver.hpp"

template<unsigned DIM>
AbstractUnstructuredGridDiscreteContinuumSolver<DIM>::AbstractUnstructuredGridDiscreteContinuumSolver()
    :   AbstractDiscreteContinuumSolver<DIM>(),
        mpMesh()
{

}

template<unsigned DIM>
AbstractUnstructuredGridDiscreteContinuumSolver<DIM>::~AbstractUnstructuredGridDiscreteContinuumSolver()
{

}

template<unsigned DIM>
std::vector<QConcentration > AbstractUnstructuredGridDiscreteContinuumSolver<DIM>::GetConcentrationsAtCentroids()
{
    return this->GetConcentrations(this->mpDensityMap->GetGridCalculator()->GetGrid()->GetCellLocations());
}

template<unsigned DIM>
void AbstractUnstructuredGridDiscreteContinuumSolver<DIM>::Setup()
{
    if(!this->mpDensityMap)
    {
        EXCEPTION("Mesh needed before Setup can be called.");
    }

    mpMesh = std::dynamic_pointer_cast<DiscreteContinuumMesh<DIM> >(this->mpDensityMap->GetGridCalculator()->GetGrid());
    if(!mpMesh)
    {
        EXCEPTION("Can't cast to mesh during Setup");
    }

    // Set up the VTK solution
    this->mSolution = std::vector<double>(this->mpDensityMap->GetGridCalculator()->GetGrid()->GetNumberOfPoints(), 0.0);
}

template<unsigned DIM>
void AbstractUnstructuredGridDiscreteContinuumSolver<DIM>::UpdateSolution(const std::vector<double>& data)
{
    if(this->mSolution.size()==0)
    {
        this->Setup();
    }

    this->mpDensityMap->GetGridCalculator()->GetGrid()->AddPointData(data, this->GetLabel());

    // Note, if the data vector being passed in is mPointSolution, then it will be overwritten with zeros.
    this->mSolution = std::vector<double>(data.size(), 0.0);
    for (unsigned i = 0; i < data.size(); i++)
    {
        this->mSolution[i] = data[i];
    }
    this->mConcentrations = std::vector<QConcentration >(data.size(), 0.0*this->mReferenceConcentration);
    for (unsigned i = 0; i < data.size(); i++)
    {
        this->mConcentrations[i] = data[i]*this->mReferenceConcentration;
    }
}

template<unsigned DIM>
void AbstractUnstructuredGridDiscreteContinuumSolver<DIM>::UpdateElementSolution(const std::vector<double>& data)
{
    if(this->mSolution.size()==0)
    {
        this->Setup();
    }
    this->mpDensityMap->GetGridCalculator()->GetGrid()->AddCellData(data, this->GetLabel());
}

template<unsigned DIM>
void AbstractUnstructuredGridDiscreteContinuumSolver<DIM>::UpdateSolution(const std::vector<QConcentration >& data)
{
    if(this->mSolution.size()==0)
    {
        this->Setup();
    }

    // Note, if the data vector being passed in is mPointSolution, then it will be overwritten with zeros.
    this->mSolution = std::vector<double>(data.size(), 0.0);
    for (unsigned i = 0; i < data.size(); i++)
    {
        this->mSolution[i] = data[i]/this->mReferenceConcentration;
    }

    this->mpDensityMap->GetGridCalculator()->GetGrid()->AddPointData(this->mSolution, this->GetLabel());

    // Note, if the data vector being passed in is mPointSolution, then it will be overwritten with zeros.
    this->mConcentrations = std::vector<QConcentration >(data.size(), 0.0*this->mReferenceConcentration);
    for (unsigned i = 0; i < data.size(); i++)
    {
        this->mConcentrations[i] = data[i];
    }
}

template<unsigned DIM>
void AbstractUnstructuredGridDiscreteContinuumSolver<DIM>::UpdateCellData()
{
    if(this->mSolution.size()==0)
    {
        this->Setup();
    }
    std::vector<std::vector<CellPtr> > point_cell_map = this->mpDensityMap->GetGridCalculator()->rGetCellMap();
    for(unsigned idx=0; idx<point_cell_map.size(); idx++)
    {
        for(unsigned jdx=0; jdx<point_cell_map[idx].size(); jdx++)
        {
            point_cell_map[idx][jdx]->GetCellData()->SetItem(this->mLabel, this->mSolution[idx]);
        }
    }
}

template<unsigned DIM>
void AbstractUnstructuredGridDiscreteContinuumSolver<DIM>::Update()
{

}

template<unsigned DIM>
void AbstractUnstructuredGridDiscreteContinuumSolver<DIM>::Write()
{
    if(!this->mpOutputFileHandler)
    {
        EXCEPTION("An output file handler has not been set for the DiscreteContinuum solver.");
    }

    if(PetscTools::AmMaster())
    {
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> p_writer1 = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        if(!this->mFilename.empty())
        {
            p_writer1->SetFileName((this->mpOutputFileHandler->GetOutputDirectoryFullPath() + "/" + this->mFilename+".vtu").c_str());
        }
        else
        {
            p_writer1->SetFileName((this->mpOutputFileHandler->GetOutputDirectoryFullPath() + "/solution.vtu").c_str());
        }
        #if VTK_MAJOR_VERSION <= 5
            p_writer1->SetInput(vtkUnstructuredGrid::SafeDownCast(this->GetVtkSolution()));
        #else
            p_writer1->SetInputData(vtkUnstructuredGrid::SafeDownCast(this->GetVtkSolution()));
        #endif
        p_writer1->Update();
        p_writer1->Write();
    }
}

// Explicit instantiation
template class AbstractUnstructuredGridDiscreteContinuumSolver<2> ;
template class AbstractUnstructuredGridDiscreteContinuumSolver<3> ;
