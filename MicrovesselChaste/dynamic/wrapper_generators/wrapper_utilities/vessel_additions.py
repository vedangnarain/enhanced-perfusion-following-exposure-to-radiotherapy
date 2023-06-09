#!/usr/bin/env python

"""Copyright (c) 2005-2016, University of Oxford.
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
"""

import sys
from pyplusplus import module_builder
from pyplusplus.module_builder import call_policies
from pygccxml import parser

def update_builder(builder, class_collection):


    # Default template args problem
#     builder.class_('VesselNetworkGenerator<2>').calldefs().use_default_arguments=False    
#     builder.class_('VesselNetworkGenerator<3>').calldefs().use_default_arguments=False    

    # The VesselSegment and Vessel classes need factory constructors as they have private constructor methods
    builder.class_('VesselSegment<3>').add_declaration_code('boost::shared_ptr<VesselSegment<3> > (*VS3_Nodes)(boost::shared_ptr<VesselNode<3> >, boost::shared_ptr<VesselNode<3> >) = &VesselSegment<3>::Create;')
    builder.class_('VesselSegment<3>').add_declaration_code('boost::shared_ptr<VesselSegment<3> > (*VS3_Copy)(boost::shared_ptr<VesselSegment<3> >) = &VesselSegment<3>::Create;')
    builder.class_('Vessel<3>').add_declaration_code('boost::shared_ptr<Vessel<3> > (*V3_SingleSegment)(boost::shared_ptr<VesselSegment<3> >) = &Vessel<3>::Create;')
    builder.class_('Vessel<3>').add_declaration_code('boost::shared_ptr<Vessel<3> > (*V3_MultiSegment)(std::vector<boost::shared_ptr<VesselSegment<3> > >) = &Vessel<3>::Create;')
    builder.class_('Vessel<3>').add_declaration_code('boost::shared_ptr<Vessel<3> > (*V3_Nodes)(std::vector<boost::shared_ptr<VesselNode<3> > >) = &Vessel<3>::Create;')
    builder.class_('VesselSegment<3>').add_registration_code('def("__init__", bp::make_constructor(VS3_Nodes))')
    builder.class_('VesselSegment<3>').add_registration_code('def("__init__", bp::make_constructor(VS3_Copy))')
    builder.class_('Vessel<3>').add_registration_code('def("__init__", bp::make_constructor(V3_SingleSegment))')
    builder.class_('Vessel<3>').add_registration_code('def("__init__", bp::make_constructor(V3_MultiSegment))')
    builder.class_('Vessel<3>').add_registration_code('def("__init__", bp::make_constructor(V3_Nodes))')

    return builder, class_collection