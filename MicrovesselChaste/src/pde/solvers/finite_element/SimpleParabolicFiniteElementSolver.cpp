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
#include "AbstractDiscreteContinuumParabolicPde.hpp"
#include "Exception.hpp"
#include "SimpleParabolicFiniteElementSolver.hpp"
#include "SimpleLinearParabolicSolverWithStorage.hpp"

template<unsigned DIM>
SimpleParabolicFiniteElementSolver<DIM>::SimpleParabolicFiniteElementSolver()
    : AbstractFiniteElementSolverBase<DIM>(),
      mIntermediateSolutionCollection(),
      mIntermediateSolutionFrequency(1),
      mStoreIntermediate(false),
      mWriteIntermediate(false),
      mTimeIncrement(0.001),
      mSolveStartTime(0.0),
      mSolveEndTime(1.0),
      mInitialGuess()
{

}

template<unsigned DIM>
SimpleParabolicFiniteElementSolver<DIM>::~SimpleParabolicFiniteElementSolver()
{

}

template <unsigned DIM>
std::shared_ptr<SimpleParabolicFiniteElementSolver<DIM> > SimpleParabolicFiniteElementSolver<DIM>::Create()
{
    return std::make_shared<SimpleParabolicFiniteElementSolver<DIM> >();

}

template <unsigned DIM>
void SimpleParabolicFiniteElementSolver<DIM>::SetTargetTimeIncrement(double targetIncrement)
{
    mTimeIncrement = targetIncrement;
}

template <unsigned DIM>
void SimpleParabolicFiniteElementSolver<DIM>::SetStartTime(double startTime)
{
    mSolveStartTime = startTime;
}

template <unsigned DIM>
void SimpleParabolicFiniteElementSolver<DIM>::SetEndTime(double endTime)
{
    mSolveEndTime = endTime;
}

template <unsigned DIM>
void SimpleParabolicFiniteElementSolver<DIM>::SetInitialGuess(const std::vector<double>& rInitialGuess)
{
    mInitialGuess = rInitialGuess;
}

template <unsigned DIM>
const std::vector<std::pair<std::vector<double>, double> >& SimpleParabolicFiniteElementSolver<DIM>::rGetIntermediateSolutions()
{
    return mIntermediateSolutionCollection;
}

template <unsigned DIM>
void SimpleParabolicFiniteElementSolver<DIM>::SetStoreIntermediateSolutions(bool store, unsigned frequency)
{
    mStoreIntermediate = store;
    mIntermediateSolutionFrequency = frequency;
}

template <unsigned DIM>
void SimpleParabolicFiniteElementSolver<DIM>::SetWriteIntermediateSolutions(bool write, unsigned frequency)
{
    mWriteIntermediate = write;
    mStoreIntermediate = write;
    mIntermediateSolutionFrequency = frequency;
}

template<unsigned DIM>
void SimpleParabolicFiniteElementSolver<DIM>::Solve()
{
    AbstractFiniteElementSolverBase<DIM>::Solve();

    // Set up the boundary conditions in the Chaste format
    std::shared_ptr<BoundaryConditionsContainer<DIM, DIM, 1> > p_bcc =
            std::shared_ptr<BoundaryConditionsContainer<DIM, DIM, 1> >(new BoundaryConditionsContainer<DIM, DIM, 1> );

    for(unsigned idx=0; idx<this->mBoundaryConditions.size(); idx++)
    {
        this->mBoundaryConditions[idx]->SetGridCalculator(this->mpDensityMap->GetGridCalculator());
        this->mBoundaryConditions[idx]->UpdateBoundaryConditions(p_bcc);
    }

    // Check the type of pde
    if(std::shared_ptr<AbstractDiscreteContinuumParabolicPde<DIM, DIM> > p_parabolic_pde =
            std::dynamic_pointer_cast<AbstractDiscreteContinuumParabolicPde<DIM, DIM> >(this->mpPde))
    {
        Vec initial_guess = PetscTools::CreateAndSetVec(this->mpMesh->GetNumNodes(), 0.0);
        SimpleLinearParabolicSolverWithStorage<DIM, DIM> solver(this->mpMesh.get(), p_parabolic_pde.get(), p_bcc.get());

        /* The interface is exactly the same as the `SimpleLinearParabolicSolver`. */
        solver.SetTimeStep(mTimeIncrement);
        solver.SetTimes(mSolveStartTime, mSolveEndTime);
        solver.SetInitialCondition(initial_guess);
        solver.SetStoreIntermediateSolutions(mStoreIntermediate, mIntermediateSolutionFrequency);

        Vec result = solver.Solve();
        ReplicatableVector solution_repl(result);

        this->mSolution = std::vector<double>(solution_repl.GetSize());
        this->mConcentrations = std::vector<QConcentration >(solution_repl.GetSize());
        for(unsigned idx = 0; idx < solution_repl.GetSize(); idx++)
        {
            this->mSolution[idx] = solution_repl[idx];
            this->mConcentrations[idx] = solution_repl[idx]*this->mReferenceConcentration;
        }
        this->UpdateSolution(this->mSolution);

        if(mStoreIntermediate)
        {
            mIntermediateSolutionCollection = solver.rGetIntermediateSolutions();
        }

        // Tidy up
        PetscTools::Destroy(initial_guess);
        PetscTools::Destroy(result);
    }
    else
    {
        EXCEPTION("PDE Type could not be identified, did you set a PDE?");
    }

    if(this->mWriteSolution)
    {
        this->Write();
    }

    if(mWriteIntermediate)
    {
        std::string base_file_name;
        std::string original_file_name;
        if(this->mFilename.empty())
        {
            original_file_name = "solution";
        }
        else
        {
            original_file_name = this->mFilename;
        }
        base_file_name = original_file_name + "_intermediate_t_";

        for(unsigned idx=0;idx<mIntermediateSolutionCollection.size();idx++)
        {
            this->UpdateSolution(mIntermediateSolutionCollection[idx].first);
            this->mFilename = base_file_name + boost::lexical_cast<std::string>(unsigned(100.0*mIntermediateSolutionCollection[idx].second/mTimeIncrement));
            this->Write();
        }
        this->mFilename = original_file_name;
    }
}

// Explicit instantiation
template class SimpleParabolicFiniteElementSolver<2>;
template class SimpleParabolicFiniteElementSolver<3>;
