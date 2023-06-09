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

#ifndef TESTUNITCOLLECTION_HPP
#define TESTUNITCOLLECTION_HPP

#include <cxxtest/TestSuite.h>
#include <SmartPointers.hpp>
#include "UnitCollection.hpp"
#include "VectorUnitCollection.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
 * Check units are instantiated correctly and that they stream out the correct amount.
 */
class TestUnitCollection : public CxxTest::TestSuite
{

public:

    void TestAllUnits()
    {
        // Time
        QTime s1(3.0*unit::seconds);
        QTime m1(3.0*unit::minutes);
        QTime h1(3.0*unit::hours);
        QTime d1(3.0*unit::days);
        std::stringstream s1s;
        s1s << s1;
        std::stringstream m1s;
        m1s << m1;
        std::stringstream h1s;
        h1s << h1;
        std::stringstream d1s;
        d1s << d1;
        TS_ASSERT_EQUALS(s1s.str(), "3");
        TS_ASSERT_EQUALS(m1s.str(), "180");
        TS_ASSERT_EQUALS(h1s.str(), "10800");
        TS_ASSERT_EQUALS(d1s.str(), "259200");
    }

    void TestVectorUnits()
    {
        VecQLength<2> loc(1.0);
        std::cout << loc.Convert(1_um)[0] << std::endl;
        std::cout << loc[0] << std::endl;

        VecQLength<2> origin = QOrigin2;
        std::cout << origin.Convert(1_m)[0] << std::endl;
    }

};

#endif // TESTUNITCOLLECTION_HPP
