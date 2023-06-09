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

#ifndef TESTMAPPABLEGRIDGENERATOR_HPP_
#define TESTMAPPABLEGRIDGENERATOR_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#include <string>
#include "CheckpointArchiveTypes.hpp"
#include "ArchiveLocationInfo.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "MappableGridGenerator.hpp"
#include "Part.hpp"
#include "DiscreteContinuumMesh.hpp"
#include "DiscreteContinuumMeshGenerator.hpp"
#include "Vertex.hpp"
#include "PetscTools.hpp"
#include "MultiFormatMeshWriter.hpp"

#include "PetscAndVtkSetupAndFinalize.hpp"

class TestMappableGridGenerator : public CxxTest::TestSuite
{

public:

    void TestMakePlane()
    {
        std::string output_directory = "TestMappableGridGenerator";
        if(PetscTools::IsParallel())
        {
            output_directory += "Parallel";
        }
        OutputFileHandler output_file_handler(output_directory);

        // Make one with and without end caps
        MappableGridGenerator<3> generator;
        std::shared_ptr<Part<3> > p_part = generator.GeneratePlane(10, 10);
        TS_ASSERT_EQUALS(p_part->GetVertices().size(), 200u);
        TS_ASSERT_EQUALS(p_part->GetPolygons().size(), 198u);

        TS_ASSERT_THROWS_THIS(generator.GeneratePlane(0, 10), "The number of points in X and Y must be greater than 0.");

        std::shared_ptr<Part<3> > p_part_no_caps = generator.GeneratePlane(10, 10, false, false);
        TS_ASSERT_EQUALS(p_part_no_caps->GetVertices().size(), 200u);
        TS_ASSERT_EQUALS(p_part_no_caps->GetPolygons().size(), 180u);

        // Make sure the resulting part can be meshed
        auto p_mesh_generator = DiscreteContinuumMeshGenerator<3>::Create();
        p_mesh_generator->SetDomain(p_part);
        p_mesh_generator->Update();

        MultiFormatMeshWriter<3> mesh_writer;
        mesh_writer.SetMesh(p_mesh_generator->GetMesh());
        mesh_writer.SetFileName(output_file_handler.GetOutputDirectoryFullPath()+"plane");
        mesh_writer.Write();
    }

    void TestMakeCylinder()
    {
        std::string output_directory = "TestMappableGridGenerator";
        if(PetscTools::IsParallel())
        {
            output_directory += "Parallel";
        }
        OutputFileHandler output_file_handler(output_directory, false);

        // Make one closed cylinder, one open cylinder and one with too large an angle.
        MappableGridGenerator<3> generator;
        QLength radius = 1.5_m;
        QLength thickness = 0.1_m;
        QLength height = 5.0_m;

        TS_ASSERT_THROWS_THIS(generator.GenerateCylinder(0_m, thickness, height, 10, 10),
                "The cylinder radius and height must be greater than 0.0");

        std::shared_ptr<Part<3> > p_part = generator.GenerateCylinder(radius, thickness, height, 10, 10);
        std::shared_ptr<Part<3> > p_part_open = generator.GenerateCylinder(radius, thickness, height, 10, 10, M_PI);
        TS_ASSERT_THROWS_ANYTHING(generator.GenerateCylinder(radius, thickness, height, 10, 10, 2.1*M_PI));

        // Make sure the vertices are in the expected locations
        std::vector<std::shared_ptr<Vertex<3> > > vertices = p_part->GetVertices();
        for(unsigned idx=0; idx<vertices.size(); idx++)
        {
            double loc_x = vertices[idx]->Convert(1.0*unit::metres)[0];
            double loc_z = vertices[idx]->Convert(1.0*unit::metres)[2];
            double distance = std::sqrt(loc_x*loc_x + loc_z*loc_z);
            bool is_inside = (distance < 1.5  + 1.e-6) && (distance > 1.4  - 1.e-6);
            TS_ASSERT(is_inside);
        }

        // Make sure the part can be meshed
        std::shared_ptr<DiscreteContinuumMeshGenerator<3> > p_mesh_generator = DiscreteContinuumMeshGenerator<3>::Create();
        p_mesh_generator->SetDomain(p_part);
        p_mesh_generator->Update();
        MultiFormatMeshWriter<3> mesh_writer;
        mesh_writer.SetMesh(p_mesh_generator->GetMesh());
        mesh_writer.SetFileName(output_file_handler.GetOutputDirectoryFullPath()+"cylinder");
        mesh_writer.Write();
    }

    void TestMakeCylinderShell()
    {
        std::string output_directory = "TestMappableGridGenerator";
        if(PetscTools::IsParallel())
        {
            output_directory += "Parallel";
        }
        OutputFileHandler output_file_handler(output_directory, false);

        // Make one closed cylinder, one open cylinder and one with too large an angle.
        MappableGridGenerator<3> generator;
        QLength radius = 1.5*unit::metres;
        QLength thickness = 0.0*unit::metres;
        QLength height = 5.0*unit::metres;
        std::shared_ptr<Part<3> > p_part = generator.GenerateCylinder(radius, thickness, height, 10, 10);
        p_part->Write(output_file_handler.GetOutputDirectoryFullPath()+"cylinder_shell.vtp");
    }

    void TestMakeHemisphere()
    {
        std::string output_directory = "TestMappableGridGenerator";
        if(PetscTools::IsParallel())
        {
            output_directory += "Parallel";
        }
        OutputFileHandler output_file_handler(output_directory, false);

        // Make one good and two 'bad' hemispheres
        MappableGridGenerator<3> generator;
        QLength radius = 1.5*unit::metres;
        QLength thickness = 0.1*unit::metres;

        TS_ASSERT_THROWS_THIS(generator.GenerateHemisphere(0_m, thickness, 10, 10, M_PI, 0.5*M_PI),
                "The sphere radius must be greater than 0.0");

        std::shared_ptr<Part<3> > p_part = generator.GenerateHemisphere(radius, thickness, 10, 10, M_PI, 0.5*M_PI);
        TS_ASSERT_THROWS_ANYTHING(generator.GenerateHemisphere(radius, thickness, 10, 10, 2.0*M_PI, 0.5*M_PI));
        TS_ASSERT_THROWS_ANYTHING(generator.GenerateHemisphere(radius, thickness, 10, 10, M_PI, 1.0*M_PI));

        // Make sure the vertices are in the expected locations
        std::vector<std::shared_ptr<Vertex<3> > > vertices = p_part->GetVertices();
        for(unsigned idx=0; idx<vertices.size(); idx++)
        {
            bool is_inside = (vertices[idx]->GetNorm2()/(1.0*unit::metres) < 1.5  + 1.e-6) &&
                    (vertices[idx]->GetNorm2()/(1.0*unit::metres) > 1.4  - 1.e-6);
            TS_ASSERT(is_inside);
        }

        // Make sure the part can be meshed
        std::shared_ptr<DiscreteContinuumMeshGenerator<3> > p_mesh_generator = DiscreteContinuumMeshGenerator<3>::Create();
        p_mesh_generator->SetDomain(p_part);
        p_mesh_generator->Update();

        MultiFormatMeshWriter<3> mesh_writer;
        mesh_writer.SetMesh(p_mesh_generator->GetMesh());
        mesh_writer.SetFileName(output_file_handler.GetOutputDirectoryFullPath()+"hemisphere");
        mesh_writer.Write();
    }

    void TestMakeHemisphereShell()
    {
        std::string output_directory = "TestMappableGridGenerator";
        if(PetscTools::IsParallel())
        {
            output_directory += "Parallel";
        }
        OutputFileHandler output_file_handler(output_directory, false);

        // Make one good and two 'bad' hemispheres
        MappableGridGenerator<3> generator;
        QLength radius = 1.5*unit::metres;
        QLength thickness = 0.0*unit::metres;

        std::shared_ptr<Part<3> > p_part = generator.GenerateHemisphere(radius, thickness, 10, 10, M_PI, 0.5*M_PI);
        p_part->Write(output_file_handler.GetOutputDirectoryFullPath()+"hemisphere_shell.vtp");
    }

    void TestArchiving()
    {
        #if BOOST_VERSION >= 105600
        // Test Archiving
        OutputFileHandler handler("archive", false);
        ArchiveLocationInfo::SetArchiveDirectory(handler.FindFile(""));
        std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath("MappableGridGenerator.arch");

        BaseUnits::Instance()->SetReferenceLengthScale(2.0*unit::metres);

        // Save archive
        {
            std::shared_ptr<MappableGridGenerator<3> > p_generator = MappableGridGenerator<3>::Create();
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_generator;
        }

        // Load archive
        {
            std::shared_ptr<MappableGridGenerator<3> > p_generator_from_archive;

            // Read from this input file
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // restore from the archive
            input_arch >> p_generator_from_archive;

            // Check that we remember the reference length
            BaseUnits::Instance()->SetReferenceLengthScale(3.0*unit::metres);
            TS_ASSERT_DELTA(p_generator_from_archive->GetReferenceLengthScale()/(1.0*unit::metres), 2.0, 1.e-6);
        }
        #endif
    }
};

#endif /*TESTMAPPABLEGRIDGENERATORHPP_*/
