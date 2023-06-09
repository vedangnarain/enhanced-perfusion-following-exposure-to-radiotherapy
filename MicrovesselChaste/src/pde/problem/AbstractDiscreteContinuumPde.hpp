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

#ifndef ABSTRACTDISCRETECONTINUUMPDE_HPP_
#define ABSTRACTDISCRETECONTINUUMPDE_HPP_

#include <memory>
#include <string>
#include "UblasIncludes.hpp"
#include "UblasVectorInclude.hpp"
#include "GeometryTools.hpp"
#include "UnitCollection.hpp"
#include "DiscreteSource.hpp"

/**
* Discrete source forward declared because it can contain its own DiscreteContinuumSolvers
*/
template<unsigned DIM>
class DiscreteSource;

/**
 * This class specifies interfaces that PDEs being used in Discrete Continuum solvers should
 * implement, mostly related to interaction with discrete sinks and sources.
 * It can be combined with a range of abstract PDEs in Chaste or used to build
 * more general PDE descriptions.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractDiscreteContinuumPde
{

protected:

    /**
     * The diffusion constant for isotropic diffusion
     */
    QDiffusivity mDiffusivity;

    /**
     * The collection of discrete sources for addition to the continuum terms
     */
    std::vector<std::shared_ptr<DiscreteSource<SPACE_DIM> > > mDiscreteSources;

    /**
     * This is used internally to scale concentrations before and after linear system solves, reads and writes.
     * Since those functions don't use Boost Units. It should not affect the solution, but can be judiciously chosen
     * to avoid precision problems.
     */
    QConcentration mReferenceConcentration;

    /**
     * The reference length scale, used to scale diffusivity to the mesh size
     */
    QLength mReferenceLengthScale;

    /**
     * The reference time scale, used to scale diffusivity
     */
    QTime mReferenceTimeScale;

public:

    /**
     * Constructor
     */
    AbstractDiscreteContinuumPde();

    /**
     * Desctructor
     */
    virtual ~AbstractDiscreteContinuumPde();

    /**
     * Add a discrete source to the pde
     * @param pDiscreteSource a pointer the discrete source
     */
    void AddDiscreteSource(std::shared_ptr<DiscreteSource<SPACE_DIM> > pDiscreteSource);

    /**
     * Return the diffusion constant for isotropic diffusion
     * @return the diffusion constant
     */
    QDiffusivity ComputeIsotropicDiffusionTerm();

    /**
     * Return the collection of discrete sources
     * @return vector of pointers to the discrete sources
     */
    std::vector<std::shared_ptr<DiscreteSource<SPACE_DIM> > > GetDiscreteSources();

    /**
     * Set the isotropic diffusion constant
     * @param diffusivity the isotropic diffusion constant
     */
    void SetIsotropicDiffusionConstant(QDiffusivity diffusivity);

    /**
     * Set the reference concentration
     * @param referenceConcentration the reference concentration
     */
    void SetReferenceConcentration(QConcentration referenceConcentration);

    /**
     * Set the reference length
     * @param referenceLength the reference length
     */
    void SetReferenceLength(QLength referenceLength);

    /**
     * Set the reference time
     * @param referenceTime the reference time
     */
    void SetReferenceTime(QTime referenceTime);

    /**
     * Update the discrete source strengths
     */
    virtual void UpdateDiscreteSourceStrengths();

};

#endif /*ABSTRACTDISCRETECONTINUUMPDE_HPP_*/
