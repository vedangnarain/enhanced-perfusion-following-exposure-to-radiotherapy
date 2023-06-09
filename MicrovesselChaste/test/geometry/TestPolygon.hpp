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

#ifndef TESTPOLYGON_HPP_
#define TESTPOLYGON_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning
#include <vtkSmartPointer.h>
#include <vtkPolygon.h>
#include <vtkPoints.h>
#include <vtkPlane.h>
#include "CheckpointArchiveTypes.hpp"
#include "ArchiveLocationInfo.hpp"
#include "OutputFileHandler.hpp"
#include "UblasIncludes.hpp"
#include "SmartPointers.hpp"
#include "Vertex.hpp"
#include "Polygon.hpp"

#include "PetscAndVtkSetupAndFinalize.hpp"


std::vector<VertexPtr<3> > SetUpTriangle()
{
    std::vector<VertexPtr<3> > vertices;
    vertices.push_back(Vertex<3>::Create(0.0_um, 0.0_um));
    vertices.push_back(Vertex<3>::Create(1.0_um, 0.0_um));
    vertices.push_back(Vertex<3>::Create(1.0_um, 1.0_um));
    return vertices;
}

template<unsigned DIM>
std::vector<VertexPtr<DIM> > SetUpSquare()
{
    std::vector<VertexPtr<DIM> > vertices;
    vertices.push_back(Vertex<DIM>::Create(0.0_um, 0.0_um));
    vertices.push_back(Vertex<DIM>::Create(1.0_um, 0.0_um));
    vertices.push_back(Vertex<DIM>::Create(1.0_um, 1.0_um));
    vertices.push_back(Vertex<DIM>::Create(0.0_um, 1.0_um));
    return vertices;
}

class TestPolygon : public CxxTest::TestSuite
{
public:

    void TestConstructor()
    {
        std::vector<VertexPtr<3> > vertices = SetUpTriangle();
        Polygon<3> polygon1(vertices);
        Polygon<3> polygon2(vertices[0]);
        TS_ASSERT_EQUALS(polygon1.rGetVertices().size(), 3u);
        TS_ASSERT_EQUALS(polygon2.rGetVertices().size(), 1u);
    }

    void TestFactoryConstructor()
    {
        std::vector<VertexPtr<3> > vertices = SetUpTriangle();
        auto p_polygon1 = Polygon<3>::Create(vertices);
        auto p_polygon2 = Polygon<3>::Create(vertices[0]);
        TS_ASSERT_EQUALS(p_polygon1->rGetVertices().size(), 3u);
        TS_ASSERT_EQUALS(p_polygon2->rGetVertices().size(), 1u);
    }

    void TestAddingVertices()
    {
        std::vector<VertexPtr<3> > vertices = SetUpTriangle();

        std::vector<VertexPtr<3> > new_vertices;
        new_vertices.push_back(Vertex<3>::Create(1.0_um, 1.0_um));
        new_vertices.push_back(Vertex<3>::Create(1.0_um, 2.0_um));

        auto p_polygon = Polygon<3>::Create(vertices);
        p_polygon->AddVertices(new_vertices);
        TS_ASSERT_EQUALS(p_polygon->rGetVertices().size(), 5u);

        auto p_polygon2 = Polygon<3>::Create(vertices);
        p_polygon2->AddVertex(new_vertices[0]);
        TS_ASSERT_EQUALS(p_polygon2->rGetVertices().size(), 4u);

        TS_ASSERT_THROWS_THIS(p_polygon2->GetVertex(100), "Requested vertex index out of range");
        TS_ASSERT_THROWS_THIS(p_polygon2->ReplaceVertex(100, vertices[0]), "Requested vertex index out of range");
    }

    void TestVtkMethods()
    {
        std::vector<VertexPtr<3> > vertices = SetUpSquare<3>();
        auto p_polygon = Polygon<3>::Create(vertices);

        c_vector<double, 3> centroid = p_polygon->GetCentroid().Convert(1_um);
        TS_ASSERT_DELTA(centroid[0], 0.5, 1.e-6);
        TS_ASSERT_DELTA(centroid[1], 0.5, 1.e-6);
        TS_ASSERT_DELTA(centroid[2], 0.0, 1.e-6);

        c_vector<double, 3> normal = p_polygon->GetNormal();
        TS_ASSERT_DELTA(normal[0], 0.0, 1.e-6);
        TS_ASSERT_DELTA(normal[1], 0.0, 1.e-6);
        TS_ASSERT_DELTA(normal[2], 1.0, 1.e-6);

        TS_ASSERT(p_polygon->ContainsPoint(Vertex<3>(0.75_um, 0.75_um, 0.0_um)));
        TS_ASSERT(!p_polygon->ContainsPoint(Vertex<3>(1.25_um, 0.75_um, 0.0_um)));
        TS_ASSERT(!p_polygon->ContainsPoint(Vertex<3>(0.75_um, 0.75_um, 1.0_um)));

        // Cover all bounding box cases for point containment
        c_vector<double, 3> rotation_axis = unit_vector<double>(3, 0);
        p_polygon->RotateAboutAxis(rotation_axis, M_PI/2.0);
        TS_ASSERT(!p_polygon->ContainsPoint(Vertex<3>(10_um, 10_um, 10_um)));
        p_polygon->RotateAboutAxis(rotation_axis, -M_PI/2.0);

        rotation_axis = unit_vector<double>(3, 1);
        p_polygon->RotateAboutAxis(rotation_axis, M_PI/2.0);
        TS_ASSERT(!p_polygon->ContainsPoint(Vertex<3>(10_um, 10_um, 10_um)));

        vtkSmartPointer<vtkPolygon> p_vtk_polygon =  p_polygon->GetVtkPolygon();
        std::pair<vtkSmartPointer<vtkPoints>, vtkSmartPointer<vtkIdTypeArray> > vertex_pair = p_polygon->GetVtkVertices();
    }

    void TestVtkMethods2d()
    {
        std::vector<VertexPtr<2> > vertices = SetUpSquare<2>();

        auto p_polygon = Polygon<2>::Create(vertices);

        c_vector<double, 2> centroid = p_polygon->GetCentroid().Convert(1_um);
        TS_ASSERT_DELTA(centroid[0], 0.5, 1.e-6);
        TS_ASSERT_DELTA(centroid[1], 0.5, 1.e-6);

        c_vector<double, 3> normal = p_polygon->GetNormal();
        TS_ASSERT_DELTA(normal[0], 0.0, 1.e-6);
        TS_ASSERT_DELTA(normal[1], 0.0, 1.e-6);

        TS_ASSERT(p_polygon->ContainsPoint(Vertex<2>(0.75_um, 0.75_um, 0.0_um)));
        TS_ASSERT(!p_polygon->ContainsPoint(Vertex<2>(1.25_um, 0.75_um, 0.0_um)));

        vtkSmartPointer<vtkPolygon> p_vtk_polygon =  p_polygon->GetVtkPolygon();
        std::pair<vtkSmartPointer<vtkPoints>, vtkSmartPointer<vtkIdTypeArray> > vertex_pair = p_polygon->GetVtkVertices();
    }

    void TestTransforms()
    {
        std::vector<VertexPtr<3> > vertices = SetUpSquare<3>();
        auto p_polygon = Polygon<3>::Create(vertices);

        Vertex<3> translation_vector(2.0_um, 2.0_um, 0.0_um);
        Vertex<3> new_position = *vertices[1] + translation_vector;

        p_polygon->Translate(translation_vector);
        TS_ASSERT_DELTA(p_polygon->rGetVertices()[1]->Convert(1_um)[0], 3.0, 1.e-6);

        c_vector<double, 3> rotation_axis = unit_vector<double>(3, 2);
        p_polygon->RotateAboutAxis(rotation_axis, M_PI/2.0);
        TS_ASSERT_DELTA(p_polygon->rGetVertices()[1]->Convert(1_um)[1], 3.0, 1.e-6);
    }

    void TestGeometryOperations()
    {
        std::vector<VertexPtr<3> > vertices = SetUpSquare<3>();
        auto p_polygon = Polygon<3>::Create(vertices);

        std::array<QLength, 6> bbox = p_polygon->GetBoundingBox();
        TS_ASSERT_DELTA(bbox[0]/1_um, 0.0, 1.e-6);
        TS_ASSERT_DELTA(bbox[1]/1_um, 1.0, 1.e-6);
        TS_ASSERT_DELTA(bbox[2]/1_um, 0.0, 1.e-6);
        TS_ASSERT_DELTA(bbox[3]/1_um, 1.0, 1.e-6);
        TS_ASSERT_DELTA(bbox[4]/1_um, 0.0, 1.e-6);
        TS_ASSERT_DELTA(bbox[5]/1_um, 0.0, 1.e-6);

        TS_ASSERT_DELTA(p_polygon->GetDistance(Vertex<3>(0.5_um, 0.5_um, 0.5_um))/1_um, 0.5, 1.e-6);
        TS_ASSERT_DELTA(p_polygon->GetPlane()->GetNormal()[0], 0.0, 1.e-6);
        TS_ASSERT_DELTA(p_polygon->GetPlane()->GetNormal()[1], 0.0, 1.e-6);
        TS_ASSERT_DELTA(std::abs(p_polygon->GetPlane()->GetNormal()[2]), 1.0, 1.e-6);

        TS_ASSERT_DELTA(p_polygon->GetDistanceToEdges(Vertex<3>(1.5_um, 0.5_um, 0.0_um))/1_um, 0.5, 1.e-6);
    }

    void TestLabelling()
    {
        std::vector<VertexPtr<3> > vertices = SetUpSquare<3>();
        auto p_polygon = Polygon<3>::Create(vertices);

        std::string test = "Test";
        bool edge_found = p_polygon->AddAttributeToEdgeIfFound(Vertex<3>(0.5_um), test, 2.0);
        std::vector<std::map<std::string, double> > edge_attributes = p_polygon->rGetEdgeAttributes();
        TS_ASSERT(edge_found);
        TS_ASSERT_DELTA(edge_attributes[0]["Test"], 2.0, 1.e-6);
        TS_ASSERT(p_polygon->EdgeHasAttribute(Vertex<3>(0.5_um), "Test"));
        TS_ASSERT(!p_polygon->EdgeHasAttribute(Vertex<3>(20_um), "Test"));

        std::string poly_label = "TestPoly";
        p_polygon->AddAttribute(poly_label, 2.0);
    }

    void TestArchiving()
    {
#if BOOST_VERSION >= 105600
        // Test Archiving
        OutputFileHandler handler("archive", false);
        ArchiveLocationInfo::SetArchiveDirectory(handler.FindFile(""));
        std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath("Polygon.arch");

        std::vector<VertexPtr<3> > vertices;
        vertices.push_back(Vertex<3>::Create(0.0_um, 0.0_um));
        vertices.push_back(Vertex<3>::Create(1.0_um, 0.0_um));
        vertices.push_back(Vertex<3>::Create(1.0_um, 1.0_um));
        auto p_polygon1 = Polygon<3>::Create(vertices);

        // Save archive
        {
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_polygon1;
        }

        // Load archive
        {
            std::shared_ptr<Polygon<3> > p_polygon_from_archive;

            // Read from this input file
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // restore from the archive
            input_arch >> p_polygon_from_archive;

            // Check that we remember the reference length
            TS_ASSERT_EQUALS(p_polygon_from_archive->rGetVertices().size(), 3u);
        }
#endif
    }
};

#endif /*TESTPOLYGON_HPP_*/
