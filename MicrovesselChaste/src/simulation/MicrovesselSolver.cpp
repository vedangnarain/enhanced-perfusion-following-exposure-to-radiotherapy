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

#include <boost/lexical_cast.hpp>
#include "UblasIncludes.hpp"
#include "VesselSegment.hpp"
#include "VesselNode.hpp"
#include "MicrovesselSolver.hpp"
#include "VesselNetworkWriter.hpp"
#include "SolutionDependentDiscreteSource.hpp"
#include "AbstractDiscreteContinuumPde.hpp"

template<unsigned DIM>
MicrovesselSolver<DIM>::MicrovesselSolver() :
        mpNetwork(),
        mOutputFrequency(1),
        mpOutputFileHandler(),
        mDiscreteContinuumSolvers(),
        mpStructuralAdaptationSolver(),
        mpAngiogenesisSolver(),
        mpRegressionSolver(),
        mDiscreteContinuumSolversHaveCompatibleGridIndexing(false),
        mUpdatePdeEachSolve(true),
        mMicrovesselModifiers(),
        mpCellPopulation()
{

}

template<unsigned DIM>
MicrovesselSolver<DIM>::~MicrovesselSolver()
{

}

template<unsigned DIM>
std::shared_ptr<MicrovesselSolver<DIM> > MicrovesselSolver<DIM>::Create()
{
    return std::make_shared<MicrovesselSolver<DIM> >();

}

template<unsigned DIM>
void MicrovesselSolver<DIM>::SetVesselNetwork(std::shared_ptr<VesselNetwork<DIM> > pNetwork)
{
    mpNetwork = pNetwork;
}

template<unsigned DIM>
void MicrovesselSolver<DIM>::SetUpdatePdeEachSolve(bool doUpdate)
{
    mUpdatePdeEachSolve = doUpdate;
}

template<unsigned DIM>
void MicrovesselSolver<DIM>::AddDiscreteContinuumSolver(std::shared_ptr<AbstractDiscreteContinuumSolver<DIM> > pDiscreteContinuumSolver)
{
    mDiscreteContinuumSolvers.push_back(pDiscreteContinuumSolver);
}

template<unsigned DIM>
void MicrovesselSolver<DIM>::AddMicrovesselModifier(std::shared_ptr<AbstractMicrovesselModifier<DIM> > pMicrovesselModifier)
{
    mMicrovesselModifiers.push_back(pMicrovesselModifier);
}

template<unsigned DIM>
std::vector<std::shared_ptr<AbstractDiscreteContinuumSolver<DIM> > > MicrovesselSolver<DIM>::GetDiscreteContinuumSolvers()
{
    return mDiscreteContinuumSolvers;
}
template<unsigned DIM>
void MicrovesselSolver<DIM>::Increment()
{
    unsigned num_steps = SimulationTime::Instance()->GetTimeStepsElapsed();
    if(num_steps==0)
    {
        for(unsigned idx=0;idx<mMicrovesselModifiers.size(); idx++)
        {
            mMicrovesselModifiers[idx]->SetVesselNetwork(mpNetwork);
            mMicrovesselModifiers[idx]->SetCellPopulation(mpCellPopulation);
            for(unsigned jdx=0; jdx<mDiscreteContinuumSolvers.size(); jdx++)
            {
                mMicrovesselModifiers[idx]->AddDiscreteContinuumSolver(mDiscreteContinuumSolvers[jdx]);
            }
            mMicrovesselModifiers[idx]->SetupSolve(mpOutputFileHandler->GetOutputDirectoryFullPath());
        }
    }

    // If there is a structural adaptation or flow problem solve it
    if(mpStructuralAdaptationSolver)
    {
        mpStructuralAdaptationSolver->UpdateFlowSolver(true);
        mpStructuralAdaptationSolver->Solve();
    }

    // If there are PDEs solve them
    if(mDiscreteContinuumSolvers.size()>0)
    {
        for(unsigned idx=0; idx<mDiscreteContinuumSolvers.size(); idx++)
        {
            if(num_steps>0 and !mUpdatePdeEachSolve)
            {
                break;
            }
            mDiscreteContinuumSolvers[idx]->Update();
            mDiscreteContinuumSolvers[idx]->SetFileName("/" + mDiscreteContinuumSolvers[idx]->GetLabel() +"_solution_" + boost::lexical_cast<std::string>(num_steps));

            // Transfer PDE solutions for coupled problems
            if(idx>0 and mDiscreteContinuumSolvers[idx]->GetPde())
            {
                for(unsigned jdx=0; jdx<mDiscreteContinuumSolvers[idx]->GetPde()->GetDiscreteSources().size(); jdx++)
                {
                    std::shared_ptr<SolutionDependentDiscreteSource<DIM> > p_solution_dep_source =
                            std::dynamic_pointer_cast<SolutionDependentDiscreteSource<DIM> >(mDiscreteContinuumSolvers[idx]->GetPde()->GetDiscreteSources()[jdx]);
                    if(p_solution_dep_source)
                    {
                        if(mDiscreteContinuumSolversHaveCompatibleGridIndexing)
                        {
                            p_solution_dep_source->SetSolution(mDiscreteContinuumSolvers[idx-1]->GetConcentrations());
                        }
                        else
                        {
                            EXCEPTION("Solution dependent PDEs only work with compatible grids at the moment.");
                        }
                    }
                }
            }
            if(mOutputFrequency > 0 && num_steps % mOutputFrequency == 0)
            {
                mDiscreteContinuumSolvers[idx]->SetWriteSolution(true);
            }
            else
            {
                mDiscreteContinuumSolvers[idx]->SetWriteSolution(false);
            }
            mDiscreteContinuumSolvers[idx]->Solve();
        }
    }

    // Do angiogenesis if the is a network and solver
    if(this->mpNetwork && mpAngiogenesisSolver)
    {
        mpAngiogenesisSolver->Increment();
    }
    // Do regression if the is a network and solver
    if(this->mpNetwork && mpRegressionSolver)
    {
        mpRegressionSolver->Increment();
    }

    // Manage vessel network output
    if (this->mpNetwork)
    {
        mpNetwork->UpdateAll();
        if (mOutputFrequency > 0 && num_steps % mOutputFrequency == 0)
        {
            std::shared_ptr<VesselNetworkWriter<DIM> > p_network_writer = VesselNetworkWriter<DIM>::Create();
            p_network_writer->SetVesselNetwork(mpNetwork);
            p_network_writer->SetFileName(
                    mpOutputFileHandler->GetOutputDirectoryFullPath() + "/VesselNetwork_inc_"
                            + boost::lexical_cast<std::string>(num_steps + 1) + ".vtp");
            p_network_writer->Write();
        }
    }
    for(unsigned idx=0;idx<mMicrovesselModifiers.size(); idx++)
    {
        mMicrovesselModifiers[idx]->UpdateAtEndOfTimeStep();
    }
}

template<unsigned DIM>
void MicrovesselSolver<DIM>::Run()
{
    if (this->mpNetwork)
    {
        std::shared_ptr<VesselNetworkWriter<DIM> > p_network_writer = VesselNetworkWriter<DIM>::Create();
        mpNetwork->UpdateAll(true);
        p_network_writer->SetVesselNetwork(mpNetwork);
        p_network_writer->SetFileName(mpOutputFileHandler->GetOutputDirectoryFullPath() + "/VesselNetwork_inc_0.vtp");
        p_network_writer->Write();
    }

    Setup();

    while (!SimulationTime::Instance()->IsFinished())
    {
        Increment();
        SimulationTime::Instance()->IncrementTimeOneStep();
    }
}

template<unsigned DIM>
void MicrovesselSolver<DIM>::SetAngiogenesisSolver(std::shared_ptr<AngiogenesisSolver<DIM> > pAngiogenesisSolver)
{
    mpAngiogenesisSolver = pAngiogenesisSolver;
}

template<unsigned DIM>
void MicrovesselSolver<DIM>::SetCellPopulation(std::shared_ptr<AbstractCellPopulation<DIM,DIM> > pCellPopulation)
{
    mpCellPopulation = pCellPopulation;
}

template<unsigned DIM>
void MicrovesselSolver<DIM>::SetOutputFileHandler(std::shared_ptr<OutputFileHandler> pFileHandler)
{
    mpOutputFileHandler = pFileHandler;
}

template<unsigned DIM>
void MicrovesselSolver<DIM>::SetOutputFrequency(unsigned frequency)
{
    mOutputFrequency = frequency;
}

template<unsigned DIM>
void MicrovesselSolver<DIM>::SetupFromModifier(AbstractCellPopulation<DIM,DIM>& rCellPopulation,
                                               QLength cellReferenceLength,
                                               QConcentration cellReferenceConcentration,
                                               const std::string& rDirectory)
{
    // Set up an output file handler
    mpOutputFileHandler = std::shared_ptr<OutputFileHandler>(new OutputFileHandler(rDirectory, false));
    std::shared_ptr<DensityMap<DIM> > p_default_density_map;

    // Set up all the DiscreteContinuum solvers
    for(unsigned idx=0; idx<mDiscreteContinuumSolvers.size(); idx++)
    {
        mDiscreteContinuumSolvers[idx]->SetCellPopulation(rCellPopulation, cellReferenceLength, cellReferenceConcentration);
        mDiscreteContinuumSolvers[idx]->SetFileHandler(mpOutputFileHandler);
        if(mpNetwork)
        {
            mDiscreteContinuumSolvers[idx]->SetVesselNetwork(mpNetwork);
        }
        mDiscreteContinuumSolvers[idx]->Setup();
    }

    Setup();
}

template<unsigned DIM>
void MicrovesselSolver<DIM>::SetDiscreteContinuumSolversHaveCompatibleGridIndexing(bool compatibleIndexing)
{
    mDiscreteContinuumSolversHaveCompatibleGridIndexing = compatibleIndexing;
}

template<unsigned DIM>
void MicrovesselSolver<DIM>::Setup()
{
    // Set up all the DiscreteContinuum solvers
    std::shared_ptr<DensityMap<DIM> > p_default_density_map;
    for(unsigned idx=0; idx<mDiscreteContinuumSolvers.size(); idx++)
    {
        mDiscreteContinuumSolvers[idx]->SetFileHandler(mpOutputFileHandler);
        if(mpNetwork)
        {
            mDiscreteContinuumSolvers[idx]->SetVesselNetwork(mpNetwork);
        }
        mDiscreteContinuumSolvers[idx]->Setup();
    }

    // Set up the flow and structural adaptation solvers
    if(mpStructuralAdaptationSolver)
    {
        mpStructuralAdaptationSolver->SetVesselNetwork(mpNetwork);
    }

    if(mpRegressionSolver)
    {
        mpRegressionSolver->SetVesselNetwork(mpNetwork);
    }

}

template<unsigned DIM>
void MicrovesselSolver<DIM>::SetStructuralAdaptationSolver(std::shared_ptr<StructuralAdaptationSolver<DIM> > pStructuralAdaptationSolver)
{
    mpStructuralAdaptationSolver = pStructuralAdaptationSolver;
}

template<unsigned DIM>
void MicrovesselSolver<DIM>::UpdateCellData(std::vector<std::string> labels)
{
    if(labels.size()==0)
    {
        //update everything
        for(unsigned jdx=0; jdx<mDiscreteContinuumSolvers.size(); jdx++)
        {
            mDiscreteContinuumSolvers[jdx]->UpdateCellData();
        }
    }
    else
    {
        for(unsigned idx=0; idx<labels.size(); idx++)
        {
            for(unsigned jdx=0; jdx<mDiscreteContinuumSolvers.size(); jdx++)
            {
                if(labels[idx]==mDiscreteContinuumSolvers[idx]->GetLabel())
                {
                    mDiscreteContinuumSolvers[jdx]->UpdateCellData();
                }
            }
        }
    }
}

template<unsigned DIM>
void MicrovesselSolver<DIM>::SetRegressionSolver(std::shared_ptr<RegressionSolver<DIM> > pRegressionSolver)
{
    mpRegressionSolver = pRegressionSolver;
}

// Explicit instantiation
template class MicrovesselSolver<2> ;
template class MicrovesselSolver<3> ;
