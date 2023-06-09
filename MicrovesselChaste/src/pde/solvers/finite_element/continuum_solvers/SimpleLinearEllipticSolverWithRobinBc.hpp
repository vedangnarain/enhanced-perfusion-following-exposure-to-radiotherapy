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

#ifndef SimpleLinearEllipticSolverWithRobinBc_HPP_
#define SimpleLinearEllipticSolverWithRobinBc_HPP_

#include "AbstractAssemblerSolverHybrid.hpp"
#include "AbstractStaticLinearPdeSolver.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "RobinConditionsSurfaceTermAssembler.hpp"

/**
 * SimpleLinearEllipticSolverWithRobinBc.
 *
 * Solver for solving AbstractLinearEllipticPdes.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class SimpleLinearEllipticSolverWithRobinBc
    : public AbstractAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, 1, NORMAL>,
      public AbstractStaticLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, 1>
{
protected:

    /** The PDE to be solved. */
    AbstractLinearEllipticPde<ELEMENT_DIM,SPACE_DIM>* mpEllipticPde;

    RobinConditionsSurfaceTermAssembler<ELEMENT_DIM,SPACE_DIM,1> mRobinConditionsSurfaceTermAssembler;

    /** Boundary conditions container */
    BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM, 1>* mpBoundaryConditions;

    /**
     * @return the term to be added to the element stiffness matrix - see AbstractFeVolumeIntegralAssembler
     *
     * grad_phi[row] . ( pde_diffusion_term * grad_phi[col])
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    virtual c_matrix<double, 1*(ELEMENT_DIM+1), 1*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double,1,SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement);

    /**
     * @return the term arising from boundary conditions to be added to the element
     * stiffness vector - see AbstractFeVolumeIntegralAssembler
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    virtual c_vector<double,1*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double,1,SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement);


    // Note: does not have to provide a ComputeVectorSurfaceTerm for surface integrals,
    // the parent AbstractAssemblerSolverHybrid assumes natural Neumann BCs and uses a
    // NaturalNeumannSurfaceTermAssembler for assembling this part of the vector.

    /**
     * Delegate to AbstractAssemblerSolverHybrid::SetupGivenLinearSystem.
     *
     * @param currentSolution The current solution which can be used in setting up
     *   the linear system if needed (NULL if there isn't a current solution)
     * @param computeMatrix Whether to compute the LHS matrix of the linear system
     *   (mainly for dynamic solves).
     */
    void SetupLinearSystem(Vec currentSolution, bool computeMatrix)
    {
        assert(this->mpLinearSystem->rGetLhsMatrix() != NULL);
        assert(this->mpLinearSystem->rGetRhsVector() != NULL);

        // Assemble the matrix and vector calling methods on AbstractFeVolumeIntegralAssembler
        this->SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix());
        this->SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), true);

        if (currentSolution != NULL)
        {
            this->SetCurrentSolution(currentSolution);
        }

        if (computeMatrix)
        {
            this->Assemble();
        }
        else
        {
            this->AssembleVector();
        }

        // Add the Robin boundary conditions
        mRobinConditionsSurfaceTermAssembler.SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix(), false);
        mRobinConditionsSurfaceTermAssembler.SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), false);
        mRobinConditionsSurfaceTermAssembler.Assemble();

        this->mpLinearSystem->FinaliseRhsVector();
        this->mpLinearSystem->SwitchWriteModeLhsMatrix();

        // add Dirichlet BCs
        mpBoundaryConditions->ApplyDirichletToLinearProblem(*(this->mpLinearSystem), true);

    //// #2033 - see Test2dHeatEquationWithPeriodicBcs in TestSimpleLinearEllipticSolver.hpp
        //mpBoundaryConditions->ApplyPeriodicBcsToLinearProblem(*pLinearSystem, true);

        this->mpLinearSystem->FinaliseRhsVector();
        this->mpLinearSystem->FinaliseLhsMatrix();
    }

public:

    /**
     * Constructor.
     *
     * @param pMesh pointer to the mesh
     * @param pPde pointer to the PDE
     * @param pBoundaryConditions pointer to the boundary conditions
     */
    SimpleLinearEllipticSolverWithRobinBc(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                               AbstractLinearEllipticPde<ELEMENT_DIM,SPACE_DIM>* pPde,
                               BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions);

    /**
     * Overloaded InitaliseForSolve() which just calls the base class but also
     * sets the matrix as symmetric and sets Conjugate Gradients as the solver
     *
     * @param initialSolution initialSolution (used in base class version of this method)
     */
    void InitialiseForSolve(Vec initialSolution = NULL);
};

#endif /*SimpleLinearEllipticSolverWithRobinBc_HPP_*/
