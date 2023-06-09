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

#ifndef TESTSURFACETOOLS_HPP_
#define TESTSURFACETOOLS_HPP_

#include <cxxtest/TestSuite.h>
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning
#include "SmartPointers.hpp"
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include "BoundaryExtractor.hpp"
#include "GeometryWriter.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "UnitCollection.hpp"

#include "PetscAndVtkSetupAndFinalize.hpp"

class TestSurfaceTools : public CxxTest::TestSuite
{
public:

    void TestExtractBoundaryWithRemove()
    {
        // Read the image from file
        OutputFileHandler file_handler1 = OutputFileHandler("TestSurfaceTools/");
        FileFinder finder = FileFinder("projects/MicrovesselChaste/test/data/surface.vtp", RelativeTo::ChasteSourceRoot);

        vtkSmartPointer<vtkXMLPolyDataReader> p_reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        p_reader->SetFileName(finder.GetAbsolutePath().c_str());
        p_reader->Update();

        auto p_extractor = BoundaryExtractor::Create();
        TS_ASSERT_THROWS_THIS(p_extractor->GetOutput(), "No output set. Did you run 'Update()' ?")

        p_extractor->SetInput(p_reader->GetOutput());
        p_extractor->SetDoSmoothing(false);
        p_extractor->Update();

        auto p_writer = GeometryWriter::Create();
        p_writer->SetFileName((file_handler1.GetOutputDirectoryFullPath()+"boundary.vtp").c_str());
        p_writer->AddInput(p_extractor->GetOutput());
        p_writer->Write();
        p_writer->ClearInputs();

        p_extractor->SetDoSmoothing(true);
        p_extractor->SetRemoveDisconnected(true);
        p_extractor->SetSmoothingLength(10.0);
        p_extractor->Update();
        p_writer->SetFileName((file_handler1.GetOutputDirectoryFullPath()+"boundary_smoothed_remove.vtp").c_str());
        p_writer->AddInput(p_extractor->GetOutput());
        p_writer->Write();
    }

    void TestExtractBoundary()
    {
        // Read the image from file
        OutputFileHandler file_handler1 = OutputFileHandler("TestSurfaceTools/");
        FileFinder finder = FileFinder("projects/MicrovesselChaste/test/data/surface.vtp", RelativeTo::ChasteSourceRoot);

        vtkSmartPointer<vtkXMLPolyDataReader> p_reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        p_reader->SetFileName(finder.GetAbsolutePath().c_str());
        p_reader->Update();

        auto p_extractor = BoundaryExtractor::Create();
        TS_ASSERT_THROWS_THIS(p_extractor->Update(), "No input set.")
        p_extractor->SetInput(p_reader->GetOutput());
        p_extractor->SetDoSmoothing(false);
        p_extractor->Update();

        auto p_writer = GeometryWriter::Create();
        p_writer->SetFileName((file_handler1.GetOutputDirectoryFullPath()+"boundary.vtp").c_str());
        p_writer->AddInput(p_extractor->GetOutput());
        p_writer->Write();
        p_writer->ClearInputs();

        p_extractor->SetDoSmoothing(true);
        p_extractor->SetSmoothingLength(10.0);
        p_extractor->Update();
        p_writer->SetFileName((file_handler1.GetOutputDirectoryFullPath()+"boundary_smoothed.vtp").c_str());
        p_writer->AddInput(p_extractor->GetOutput());
        p_writer->Write();
    }
};
#endif // TESTSURFACETOOLS_HPP_
