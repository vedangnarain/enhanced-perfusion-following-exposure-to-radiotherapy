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

#include "SmartPointers.hpp"
#include "SegmentFlowProperties.hpp"
#include "VesselNetworkPropertyManager.hpp"
#include "VesselNetworkGeometryCalculator.hpp"

template <unsigned DIM>
VesselNetworkPropertyManager<DIM>::VesselNetworkPropertyManager()
{

}

template <unsigned DIM>
VesselNetworkPropertyManager<DIM>::~VesselNetworkPropertyManager()
{

}

template <unsigned DIM>
std::shared_ptr<VesselNetworkPropertyManager<DIM> > VesselNetworkPropertyManager<DIM>::Create()
{
    return std::make_shared<VesselNetworkPropertyManager<DIM> >();

}

template <unsigned DIM>
void VesselNetworkPropertyManager<DIM>::AssignInflows(std::shared_ptr<VesselNetwork<DIM> > pNetwork,
        Vertex<DIM> location, QLength searchRadius)
{
    if(pNetwork->GetNodes().size()>0)
    {
        std::vector<std::shared_ptr<VesselNode<DIM> > > inside_nodes = VesselNetworkGeometryCalculator<DIM>::GetNodesInSphere(
                pNetwork, location, searchRadius);
        for(unsigned idx=0;idx<inside_nodes.size();idx++)
        {
            inside_nodes[idx]->GetFlowProperties()->SetIsInputNode(true);
        }
    }
}

template <unsigned DIM>
void VesselNetworkPropertyManager<DIM>::AssignOutflows(std::shared_ptr<VesselNetwork<DIM> > pNetwork,
        Vertex<DIM> location, QLength searchRadius)
{
    if(pNetwork->GetNodes().size()>0)
    {
        std::vector<std::shared_ptr<VesselNode<DIM> > > outside_nodes = VesselNetworkGeometryCalculator<DIM>::GetNodesInSphere(
                pNetwork, location, searchRadius);
        for(unsigned idx=0;idx<outside_nodes.size();idx++)
        {
            outside_nodes[idx]->GetFlowProperties()->SetIsOutputNode(true);
        }
    }
}

template <unsigned DIM>
void VesselNetworkPropertyManager<DIM>::SetInflowPressures(std::shared_ptr<VesselNetwork<DIM> > pNetwork,
        QPressure pressure)
{
    std::vector<std::shared_ptr<VesselNode<DIM> > > nodes = pNetwork->GetNodes();
    for(unsigned idx=0; idx<nodes.size(); idx++)
    {
        if(nodes[idx]->GetFlowProperties()->IsInputNode())
        {
            nodes[idx]->GetFlowProperties()->SetPressure(pressure);
        }
    }
}

template <unsigned DIM>
void VesselNetworkPropertyManager<DIM>::SetOutflowPressures(std::shared_ptr<VesselNetwork<DIM> > pNetwork, QPressure pressure)
{
    std::vector<std::shared_ptr<VesselNode<DIM> > > nodes = pNetwork->GetNodes();
    for(unsigned idx=0; idx<nodes.size(); idx++)
    {
        if(nodes[idx]->GetFlowProperties()->IsOutputNode())
        {
            nodes[idx]->GetFlowProperties()->SetPressure(pressure);
        }
    }
}

template <unsigned DIM>
void VesselNetworkPropertyManager<DIM>::CopySegmentFlowProperties(std::shared_ptr<VesselNetwork<DIM> > pNetwork, unsigned index)
{
    std::shared_ptr<SegmentFlowProperties<DIM> > properties = pNetwork->GetVesselSegments()[index]->GetFlowProperties();
    std::vector<std::shared_ptr<VesselSegment<DIM> > > segments = pNetwork->GetVesselSegments();
    typename std::vector<std::shared_ptr<VesselSegment<DIM> > >::iterator it;
    for(it = segments.begin(); it != segments.end(); it++)
    {
        (*it)->SetFlowProperties(*properties);
    }
}


template <unsigned DIM>
void VesselNetworkPropertyManager<DIM>::SetNodeRadiiFromSegments(std::shared_ptr<VesselNetwork<DIM> > pNetwork)
{
    std::vector<std::shared_ptr<VesselNode<DIM> > > nodes = pNetwork->GetNodes();
    for(unsigned idx=0; idx<nodes.size(); idx++)
    {
        QLength av_radius = 0.0 * unit::metres;
        for(unsigned jdx=0; jdx<nodes[idx]->GetNumberOfSegments(); jdx++)
        {
            av_radius = av_radius + nodes[idx]->GetSegment(jdx)->GetRadius();
        }
        av_radius = av_radius / double(nodes[idx]->GetNumberOfSegments());
        nodes[idx]->SetRadius(av_radius);
    }
    pNetwork->Modified(false, false, false);
}

template <unsigned DIM>
std::vector<std::shared_ptr<VesselNode<DIM> > > VesselNetworkPropertyManager<DIM>::GetInflowNodes(std::shared_ptr<VesselNetwork<DIM> > pNetwork)
{
    std::vector<std::shared_ptr<VesselNode<DIM> > > nodes = pNetwork->GetNodes();
    std::vector<std::shared_ptr<VesselNode<DIM> > > inflow_nodes;
    for(unsigned idx=0;idx<nodes.size();idx++)
    {
        if(nodes[idx]->GetFlowProperties()->IsInputNode())
        {
            inflow_nodes.push_back(nodes[idx]);
        }
    }
    return inflow_nodes;
}

template <unsigned DIM>
std::vector<std::shared_ptr<VesselNode<DIM> > > VesselNetworkPropertyManager<DIM>::GetOutflowNodes(std::shared_ptr<VesselNetwork<DIM> > pNetwork)
{
    std::vector<std::shared_ptr<VesselNode<DIM> > > nodes = pNetwork->GetNodes();
    std::vector<std::shared_ptr<VesselNode<DIM> > > outflow_nodes;
    for(unsigned idx=0;idx<nodes.size();idx++)
    {
        if(nodes[idx]->GetFlowProperties()->IsOutputNode())
        {
            outflow_nodes.push_back(nodes[idx]);
        }
    }
    return outflow_nodes;
}

template <unsigned DIM>
void VesselNetworkPropertyManager<DIM>::SetSegmentProperties(std::shared_ptr<VesselNetwork<DIM> > pNetwork,
        std::shared_ptr<VesselSegment<DIM> >  prototype)
{
    std::vector<std::shared_ptr<VesselSegment<DIM> > > segments = pNetwork->GetVesselSegments();
    typename std::vector<std::shared_ptr<VesselSegment<DIM> > >::iterator it;
    for(it = segments.begin(); it != segments.end(); it++)
    {
        (*it)->SetRadius(prototype->GetRadius());
        (*it)->GetFlowProperties()->SetImpedance(prototype->GetFlowProperties()->GetImpedance());
        (*it)->GetFlowProperties()->SetHaematocrit(prototype->GetFlowProperties()->GetHaematocrit());
        (*it)->GetFlowProperties()->SetFlowRate(prototype->GetFlowProperties()->GetFlowRate());
        (*it)->GetFlowProperties()->SetViscosity(prototype->GetFlowProperties()->GetViscosity());
    }
}

template <unsigned DIM>
void VesselNetworkPropertyManager<DIM>::SetNodeRadii(std::shared_ptr<VesselNetwork<DIM> > pNetwork,
        QLength radius)
{
    std::vector<std::shared_ptr<VesselNode<DIM> > > nodes = pNetwork->GetNodes();

    for(unsigned idx=0; idx<nodes.size();idx++)
    {
        nodes[idx]->SetRadius(radius);
    }
}

template <unsigned DIM>
void VesselNetworkPropertyManager<DIM>::SetSegmentRadii(std::shared_ptr<VesselNetwork<DIM> > pNetwork,
        QLength radius)
{
    std::vector<std::shared_ptr<VesselSegment<DIM> > > segments = pNetwork->GetVesselSegments();

    for(unsigned idx=0; idx<segments.size();idx++)
    {
        segments[idx]->SetRadius(radius);
    }
}

template <unsigned DIM>
void VesselNetworkPropertyManager<DIM>::SetSegmentViscosity(std::shared_ptr<VesselNetwork<DIM> > pNetwork,
        QDynamicViscosity viscosity)
{
    std::vector<std::shared_ptr<VesselSegment<DIM> > > segments = pNetwork->GetVesselSegments();
    for(unsigned idx=0; idx<segments.size(); idx++)
    {
        segments[idx]->GetFlowProperties()->SetViscosity(viscosity);
    }
}

// Explicit instantiation
template class VesselNetworkPropertyManager<2>;
template class VesselNetworkPropertyManager<3>;

