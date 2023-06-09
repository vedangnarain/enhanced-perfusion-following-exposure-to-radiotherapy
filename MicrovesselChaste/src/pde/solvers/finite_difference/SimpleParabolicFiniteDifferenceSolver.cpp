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

#include <petscts.h>
#include <petscvec.h>
#include <boost/lexical_cast.hpp>
#include "ReplicatableVector.hpp"
#include "SimpleParabolicFiniteDifferenceSolver.hpp"
#include "BaseUnits.hpp"
#include "AbstractDiscreteContinuumParabolicPde.hpp"

// Parabolic solve method interfaces, needed later.
template<unsigned DIM>
PetscErrorCode ParabolicFiniteDifferenceSolver_RHSFunction(TS ts, PetscReal t, Vec currentSolution, Vec dUdt, void* pContext);

#if ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=5 )
template<unsigned DIM>
PetscErrorCode ParabolicFiniteDifferenceSolver_ComputeJacobian(TS ts, PetscReal t, Vec currentSolution, Mat pGlobalJacobian,
                                                               Mat pPreconditioner, void *pContext);
#else
template<unsigned DIM>
PetscErrorCode ParabolicFiniteDifferenceSolver_ComputeJacobian(TS ts, PetscReal t, Vec currentSolution ,Mat* pJacobian ,
        Mat* pPreconditioner, MatStructure* pMatStructure ,void* pContext);
#endif


template<unsigned DIM>
SimpleParabolicFiniteDifferenceSolver<DIM>::SimpleParabolicFiniteDifferenceSolver()
    :   AbstractFiniteDifferenceSolverBase<DIM>(),
        mIntermediateSolutionCollection(),
        mIntermediateSolutionFrequency(1),
        mStoreIntermediate(false),
        mWriteIntermediate(false),
        mTimeIncrement(0.001),
        mSolveStartTime(0.0),
        mSolveEndTime(1.0)
{

}

template<unsigned DIM>
SimpleParabolicFiniteDifferenceSolver<DIM>::~SimpleParabolicFiniteDifferenceSolver()
{

}

template <unsigned DIM>
std::shared_ptr<SimpleParabolicFiniteDifferenceSolver<DIM> > SimpleParabolicFiniteDifferenceSolver<DIM>::Create()
{
    return std::make_shared<SimpleParabolicFiniteDifferenceSolver<DIM> >();

}

template<unsigned DIM>
void SimpleParabolicFiniteDifferenceSolver<DIM>::ComputeRHSFunction(const Vec currentGuess, Vec dUdt, TS ts)
{
    this->SetVectorToAssemble(dUdt);
    this->SetCurrentSolution(currentGuess);
    this->AssembleVector();

    // Apply BCs
    std::vector<unsigned> bc_indices;
    unsigned lo = this->mpRegularGrid->GetDistributedVectorFactory()->GetLow();
    unsigned hi = this->mpRegularGrid->GetDistributedVectorFactory()->GetHigh();
    for(unsigned idx=lo; idx<hi; idx++)
    {
        if((*(this->mpBoundaryConditions))[idx].first)
        {
            PetscVecTools::SetElement(this->mVectorToAssemble, idx, 0.0);
        }
    }
    PetscVecTools::Finalise(this->mVectorToAssemble);

    #if ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=6 )
    if(mStoreIntermediate)
    {
        PetscInt totalSteps;
#if ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=8 )
        TSGetStepNumber(ts, &totalSteps);
#else
        TSGetTotalSteps(ts, &totalSteps);
#endif
        if(totalSteps%mIntermediateSolutionFrequency==0)
        {
            // Get the current time
            PetscReal currentTime;
            TSGetTime(ts, &currentTime);

            ReplicatableVector soln_guess_repl(currentGuess);
            std::vector<double> soln(soln_guess_repl.GetSize(), 0.0);
            for (unsigned row = 0; row < soln_guess_repl.GetSize(); row++)
            {
                soln[row] = soln_guess_repl[row];
            }
            mIntermediateSolutionCollection.push_back(std::pair<std::vector<double>, double>(soln, double(currentTime)));
        }
    }
    #endif
}

template<unsigned DIM>
void SimpleParabolicFiniteDifferenceSolver<DIM>::ComputeJacobian(const Vec currentGuess, Mat* pJacobian, TS ts)
{
    this->SetMatrixToAssemble(*pJacobian);
    this->SetCurrentSolution(currentGuess);
    this->AssembleMatrix();

    // Apply the boundary conditions
    std::vector<unsigned> bc_indices;
    unsigned lo = this->mpRegularGrid->GetDistributedVectorFactory()->GetLow();
    unsigned hi = this->mpRegularGrid->GetDistributedVectorFactory()->GetHigh();
    for(unsigned idx=lo; idx<hi; idx++)
    {
        if((*this->mpBoundaryConditions)[idx].first)
        {
            bc_indices.push_back(idx);
        }
    }

    PetscMatTools::ZeroRowsWithValueOnDiagonal(this->mMatrixToAssemble, bc_indices, 1.0);
    PetscMatTools::Finalise(this->mMatrixToAssemble);
}

template <unsigned DIM>
void SimpleParabolicFiniteDifferenceSolver<DIM>::SetTargetTimeIncrement(double targetIncrement)
{
    mTimeIncrement = targetIncrement;
}

template <unsigned DIM>
void SimpleParabolicFiniteDifferenceSolver<DIM>::SetStartTime(double startTime)
{
    mSolveStartTime = startTime;
}

template <unsigned DIM>
void SimpleParabolicFiniteDifferenceSolver<DIM>::SetEndTime(double endTime)
{
    mSolveEndTime = endTime;
}

template <unsigned DIM>
const std::vector<std::pair<std::vector<double>, double> >& SimpleParabolicFiniteDifferenceSolver<DIM>::rGetIntermediateSolutions()
{
    return mIntermediateSolutionCollection;
}

template <unsigned DIM>
void SimpleParabolicFiniteDifferenceSolver<DIM>::SetStoreIntermediateSolutions(bool store, unsigned frequency)
{
    mStoreIntermediate = store;
    mIntermediateSolutionFrequency = frequency;
}

template <unsigned DIM>
void SimpleParabolicFiniteDifferenceSolver<DIM>::SetWriteIntermediateSolutions(bool write, unsigned frequency)
{
    mWriteIntermediate = write;
    mStoreIntermediate = write;
    mIntermediateSolutionFrequency = frequency;
}

template<unsigned DIM>
void SimpleParabolicFiniteDifferenceSolver<DIM>::AssembleMatrix()
{
    c_vector<unsigned, 6> extents = this->mpRegularGrid->GetExtents();
    c_vector<unsigned, 3> dimensions = this->mpRegularGrid->GetDimensions();
    QTime reference_time = BaseUnits::Instance()->GetReferenceTimeScale();
    QLength spacing = this->mpRegularGrid->GetSpacing();

    std::shared_ptr<AbstractDiscreteContinuumParabolicPde<DIM> > p_parabolic_pde =
                std::dynamic_pointer_cast<AbstractDiscreteContinuumParabolicPde<DIM> >(this->GetPde());

    double diffusion_term = (p_parabolic_pde->ComputeIsotropicDiffusionTerm() / (spacing * spacing))*reference_time;
    PetscMatTools::Zero(this->mMatrixToAssemble);
    PetscMatTools::SwitchWriteMode(this->mMatrixToAssemble);
    // From PETSc 3.6 current solution VEC is read only, which is incompatible with
    // PetscVecTools::GetElement. Copy it over instead.
    ReplicatableVector soln_guess_repl(this->mCurrentSolution);

    for (unsigned i = extents[4]; i <= extents[5]; i++) // Z
    {
        for (unsigned j = extents[2]; j <= extents[3]; j++) // Y
        {
            for (unsigned k = extents[0]; k <= extents[1]; k++) // X
            {
                unsigned grid_index = this->mpRegularGrid->GetGlobalGridIndex(k, j, i);
                double current_solution = soln_guess_repl[grid_index];

                QConcentration  current_dimensional_solution =
                        current_solution*this->GetReferenceConcentration();
                QRate sink_terms =
                        p_parabolic_pde->ComputeSourceTermPrime(grid_index, current_dimensional_solution);
                double nondim_sink_terms = sink_terms*reference_time;
                PetscMatTools::AddToElement(this->mMatrixToAssemble, grid_index, grid_index, nondim_sink_terms);
                PetscMatTools::AddToElement(this->mMatrixToAssemble, grid_index, grid_index,-2.0*double(DIM)*diffusion_term);

                // No flux at x = 0 boundary
                if (k>0)
                {
                    PetscMatTools::AddToElement(this->mMatrixToAssemble,grid_index, grid_index-1, diffusion_term);
                }
                else
                {
                    PetscMatTools::AddToElement(this->mMatrixToAssemble,grid_index, grid_index+1, diffusion_term);
                }

                // No flux at x != 0 boundary
                if (k<dimensions[0]-1 )
                {
                    PetscMatTools::AddToElement(this->mMatrixToAssemble,grid_index, grid_index+1, diffusion_term);
                }
                else
                {
                    PetscMatTools::AddToElement(this->mMatrixToAssemble,grid_index, grid_index-1, diffusion_term);
                }

                // No flux at y = 0 boundary
                if (j>0)
                {
                    PetscMatTools::AddToElement(this->mMatrixToAssemble, grid_index, grid_index - dimensions[0], diffusion_term);
                }
                else
                {
                    PetscMatTools::AddToElement(this->mMatrixToAssemble, grid_index, grid_index + dimensions[0], diffusion_term);
                }

                // No flux at y != 0 boundary
                if (j<dimensions[1]-1)
                {
                    PetscMatTools::AddToElement(this->mMatrixToAssemble, grid_index, grid_index + dimensions[0], diffusion_term);
                }
                else
                {
                    PetscMatTools::AddToElement(this->mMatrixToAssemble, grid_index, grid_index - dimensions[0], diffusion_term);
                }

                if(DIM>2)
                {
                    // No flux at z = 0 boundary
                    if (i>0)
                    {
                        PetscMatTools::AddToElement(this->mMatrixToAssemble, grid_index, grid_index - dimensions[0]*dimensions[1], diffusion_term);
                    }
                    else
                    {
                        PetscMatTools::AddToElement(this->mMatrixToAssemble, grid_index, grid_index + dimensions[0]*dimensions[1], diffusion_term);
                    }

                    // No flux at z != 0 boundary
                    if (i<dimensions[2]-1)
                    {
                        PetscMatTools::AddToElement(this->mMatrixToAssemble, grid_index, grid_index + dimensions[0]*dimensions[1], diffusion_term);
                    }
                    else
                    {
                        PetscMatTools::AddToElement(this->mMatrixToAssemble, grid_index, grid_index - dimensions[0]*dimensions[1], diffusion_term);
                    }
                }
            }
        }
    }

    PetscMatTools::Finalise(this->mMatrixToAssemble);
}

template<unsigned DIM>
void SimpleParabolicFiniteDifferenceSolver<DIM>::AssembleVector()
{
    c_vector<unsigned, 6> extents = this->mpRegularGrid->GetExtents();
    c_vector<unsigned, 3> dimensions = this->mpRegularGrid->GetDimensions();
    QTime reference_time = BaseUnits::Instance()->GetReferenceTimeScale();
    QLength spacing = this->mpRegularGrid->GetSpacing();

    std::shared_ptr<AbstractDiscreteContinuumParabolicPde<DIM, DIM> > p_parabolic_pde =
                std::dynamic_pointer_cast<AbstractDiscreteContinuumParabolicPde<DIM, DIM> >(this->GetPde());
    double diffusion_term = (p_parabolic_pde->ComputeIsotropicDiffusionTerm() / (spacing * spacing))*reference_time;

    // compute function value, given current guess
    PetscVecTools::Zero(this->mVectorToAssemble);
    // From PETSc 3.6 current solution VEC is read only, which is incompatible with
    // PetscVecTools::GetElement. Copy it over instead.
    ReplicatableVector soln_guess_repl(this->mCurrentSolution);

    unsigned num_points_xy = dimensions[0]*dimensions[1];
    for (unsigned i = extents[4]; i <= extents[5]; i++) // Z
    {
        for (unsigned j = extents[2]; j <= extents[3]; j++) // Y
        {
            for (unsigned k = extents[0]; k <= extents[1]; k++) // X
            {
                unsigned grid_index = this->mpRegularGrid->GetGlobalGridIndex(k, j, i);
                double current_solution = soln_guess_repl[grid_index];
                QConcentration  current_dimensional_solution =
                        current_solution*this->GetReferenceConcentration();
                QConcentrationFlowRate sink_terms =
                        p_parabolic_pde->ComputeSourceTerm(grid_index, current_dimensional_solution);
                double nondim_sink_terms = sink_terms*(reference_time/this->GetReferenceConcentration());
                PetscVecTools::AddToElement(this->mVectorToAssemble, grid_index, nondim_sink_terms);
                PetscVecTools::AddToElement(this->mVectorToAssemble, grid_index,-2.0*double(DIM)*diffusion_term*current_solution);

                // No flux at x bottom
                if (k > 0)
                {
                    double nbr_solution = soln_guess_repl[grid_index - 1];
                    PetscVecTools::AddToElement(this->mVectorToAssemble, grid_index, diffusion_term*nbr_solution);
                }
                else
                {
                    double nbr_solution = soln_guess_repl[grid_index + 1];
                    PetscVecTools::AddToElement(this->mVectorToAssemble, grid_index, diffusion_term*nbr_solution);
                }

                // No flux at x top
                if (k < dimensions[0] - 1)
                {
                    double nbr_solution = soln_guess_repl[grid_index + 1];
                    PetscVecTools::AddToElement(this->mVectorToAssemble, grid_index, diffusion_term*nbr_solution);
                }
                else
                {
                    double nbr_solution = soln_guess_repl[grid_index - 1];
                    PetscVecTools::AddToElement(this->mVectorToAssemble, grid_index, diffusion_term*nbr_solution);
                }

                // No flux at y bottom
                if (j > 0)
                {
                    double nbr_solution = soln_guess_repl[grid_index - dimensions[0]];
                    PetscVecTools::AddToElement(this->mVectorToAssemble, grid_index, diffusion_term*nbr_solution);
                }
                else
                {
                    double nbr_solution = soln_guess_repl[grid_index + dimensions[0]];
                    PetscVecTools::AddToElement(this->mVectorToAssemble, grid_index, diffusion_term*nbr_solution);
                }

                // No flux at y top
                if (j < dimensions[1] - 1)
                {
                    double nbr_solution = soln_guess_repl[grid_index + dimensions[0]];
                    PetscVecTools::AddToElement(this->mVectorToAssemble, grid_index, diffusion_term*nbr_solution);
                }
                else
                {
                    double nbr_solution = soln_guess_repl[grid_index - dimensions[0]];
                    PetscVecTools::AddToElement(this->mVectorToAssemble, grid_index, diffusion_term*nbr_solution);
                }

                if(DIM>2)
                {
                    // No flux at z bottom
                    if (i > 0)
                    {
                        double nbr_solution = soln_guess_repl[grid_index - num_points_xy];
                        PetscVecTools::AddToElement(this->mVectorToAssemble, grid_index, diffusion_term*nbr_solution);
                    }
                    else
                    {
                        double nbr_solution = soln_guess_repl[grid_index + num_points_xy];
                        PetscVecTools::AddToElement(this->mVectorToAssemble, grid_index, diffusion_term*nbr_solution);
                    }

                    // No flux at z top
                    if (i < dimensions[2] - 1)
                    {
                        double nbr_solution = soln_guess_repl[grid_index + num_points_xy];
                        PetscVecTools::AddToElement(this->mVectorToAssemble, grid_index, diffusion_term*nbr_solution);
                    }
                    else
                    {
                        double nbr_solution = soln_guess_repl[grid_index - num_points_xy];
                        PetscVecTools::AddToElement(this->mVectorToAssemble, grid_index, diffusion_term*nbr_solution);
                    }
                }
            }
        }
    }
    PetscVecTools::Finalise(this->mVectorToAssemble);
}

template<unsigned DIM>
void SimpleParabolicFiniteDifferenceSolver<DIM>::Solve()
{
    if(!this->IsSetupForSolve)
    {
        this->Setup();
    }
    this->Update();

    c_vector<unsigned, 3> dimensions = this->mpRegularGrid->GetDimensions();
    unsigned number_of_points = dimensions[0]*dimensions[1]*dimensions[2];

    // Initialize solution
    if(this->mSolution.size()==0)
    {
        this->mSolution = std::vector<double>(number_of_points, 0.0);
    }

    // Apply BCs
    std::vector<unsigned> bc_indices;
    unsigned lo = this->mpRegularGrid->GetDistributedVectorFactory()->GetLow();
    unsigned hi = this->mpRegularGrid->GetDistributedVectorFactory()->GetHigh();
    for(unsigned idx=lo; idx<hi; idx++)
    {
        if((*(this->mpBoundaryConditions))[idx].first)
        {
            this->mSolution[idx] = (*(this->mpBoundaryConditions))[idx].second/this->GetReferenceConcentration();
        }
    }

    Vec previous_solution = PetscTools::CreateVec(number_of_points);
    for(unsigned idx=0; idx<this->mSolution.size(); idx++)
    {
        PetscVecTools::SetElement(previous_solution, idx, this->mSolution[idx]);
    }

    TS ts; // time stepper
    SNES snes; // nonlinear solver
    TSCreate(PETSC_COMM_WORLD,&ts);
    TSSetType(ts,TSBEULER);
    TSSetProblemType(ts,TS_NONLINEAR);

    TSSetRHSFunction(ts,NULL,ParabolicFiniteDifferenceSolver_RHSFunction<DIM>, this);

    Mat jacobian; // Jacobian Matrix
    PetscInt N; // number of elements
    // Get the size of the Jacobian from the residual
    VecGetSize(previous_solution, &N);
    // 8 is the maximum number of non-zero entries per row
    PetscTools::SetupMat(jacobian, N, N, 8, PETSC_DECIDE, PETSC_DECIDE, true, false);
    TSSetRHSJacobian(ts, jacobian, jacobian, ParabolicFiniteDifferenceSolver_ComputeJacobian<DIM>,this);

    // Time stepping and SNES settings
    PetscInt time_steps_max = 1e7;
    PetscReal time_total_max = mSolveEndTime;
    TSGetSNES(ts,&snes);
    PetscReal abstol = 1.0e-50;
    PetscReal reltol = 1.0e-10;
    PetscReal stol = 1.0e-10;
    PetscInt maxits = 5000;
    PetscInt maxf = 100000;
    PetscInt maxFails = -1;

    SNESSetTolerances(snes, abstol, reltol, stol, maxits, maxf);
    SNESSetMaxLinearSolveFailures(snes, maxFails);
#if ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=8 )
    TSSetMaxSteps(ts, time_steps_max);
    TSSetMaxTime(ts, time_total_max);
#else
    TSSetDuration(ts, time_steps_max, time_total_max);
#endif
    TSSetExactFinalTime(ts,TS_EXACTFINALTIME_INTERPOLATE);

    // Set initial timestep
    PetscReal dt;
    dt = mTimeIncrement;
#if ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=8 )
    TSSetTimeStep(ts, dt);
#else
    TSSetInitialTimeStep(ts, mSolveStartTime, dt);
#endif

    // Do the solve
    TSSolve(ts, previous_solution);

    ReplicatableVector soln_repl(previous_solution);
    // Populate the solution vector
    std::vector<double> solution = std::vector<double>(number_of_points, 0.0);
    for (unsigned row = 0; row < number_of_points; row++)
    {
        solution[row] = soln_repl[row];
    }

    this->UpdateSolution(solution);

    if (this->mWriteSolution)
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

    TSDestroy(&ts);
    PetscTools::Destroy(jacobian);
    PetscTools::Destroy(previous_solution);
}

template<unsigned DIM>
PetscErrorCode ParabolicFiniteDifferenceSolver_RHSFunction(TS ts, PetscReal t, Vec currentSolution, Vec dUdt, void* pContext)
{
    // extract solver from void so that we can use locally set values from the calculator object
    SimpleParabolicFiniteDifferenceSolver<DIM>* p_solver = (SimpleParabolicFiniteDifferenceSolver<DIM>*) pContext;
    p_solver->ComputeRHSFunction(currentSolution, dUdt, ts);
    return 0;
}

#if ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=5 )
template<unsigned DIM>
PetscErrorCode ParabolicFiniteDifferenceSolver_ComputeJacobian(TS ts, PetscReal t, Vec currentSolution,
                                                               Mat pGlobalJacobian, Mat pPreconditioner, void *pContext)
{
    SimpleParabolicFiniteDifferenceSolver<DIM>* p_solver =
            (SimpleParabolicFiniteDifferenceSolver<DIM>*) pContext;
    p_solver->ComputeJacobian(currentSolution, &pGlobalJacobian, ts);

#else
template<unsigned DIM>
PetscErrorCode ParabolicFiniteDifferenceSolver_ComputeJacobian(TS ts, PetscReal t, Vec currentSolution, Mat* pJacobian,
        Mat* pPreconditioner, MatStructure* pMatStructure, void* pContext)
{
    SimpleParabolicFiniteDifferenceSolver<DIM>* p_solver =
            (SimpleParabolicFiniteDifferenceSolver<DIM>*) pContext;
    p_solver->ComputeJacobian(currentSolution, pJacobian, ts);
#endif
    return 0;
}

// Explicit instantiation
template class SimpleParabolicFiniteDifferenceSolver<2>;
template class SimpleParabolicFiniteDifferenceSolver<3>;
