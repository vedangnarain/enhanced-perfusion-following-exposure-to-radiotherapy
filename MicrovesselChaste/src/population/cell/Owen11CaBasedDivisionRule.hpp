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

#ifndef OWEN11CABASEDDIVISIONRULE_HPP_
#define OWEN11CABASEDDIVISIONRULE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCaBasedDivisionRule.hpp"
#include "CaBasedCellPopulation.hpp"
#include "VesselNetwork.hpp"
#include "UnitCollection.hpp"
#include "GridCalculator.hpp"

template<unsigned SPACE_DIM> class CaBasedCellPopulation;
template<unsigned SPACE_DIM> class AbstractCaBasedDivisionRule;

/**
 * A division rule for CA-based cell populations to divide randomly into the available neighbours.
 * If there are no spaces then the move is aborted.
 */
template <unsigned SPACE_DIM>
class Owen11CaBasedDivisionRule : public AbstractCaBasedDivisionRule<SPACE_DIM>
{
private:

    /**
     * Vessel network
     */
    std::shared_ptr<VesselNetwork<SPACE_DIM> > mpVesselNetwork;

    /**
     * Regular Grid
     */
    std::shared_ptr<GridCalculator<SPACE_DIM> > mpGridCalculator;

    /**
     * Reference length scale
     */
    QLength mReferenceLengthScale;

    /**
     * Cancer cell carrying capacity
     */
    unsigned mCancerCellCarryingCapacity;

    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCaBasedDivisionRule<SPACE_DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     */
    Owen11CaBasedDivisionRule();

    /**
     * Empty destructor.
     */
    virtual ~Owen11CaBasedDivisionRule();

    /**
     * Return whether there is room to divide at all.
     *
     * @param pParentCell  The cell to divide
     * @param rCellPopulation  The CA-based cell population
     * @return if the site is available.
     */
    virtual bool IsRoomToDivide(CellPtr pParentCell,
                                CaBasedCellPopulation<SPACE_DIM>& rCellPopulation);

    /**
     * Set the vessel network
     *
     * @param pVesselNetwork The vessel network
     */
    void SetVesselNetwork(std::shared_ptr<VesselNetwork<SPACE_DIM> > pVesselNetwork);

    /**
     * Set the regular grid
     *
     * @param mpRegularGrid The regular grid
     */
    void SetGridCalculator(std::shared_ptr<GridCalculator<SPACE_DIM> > mpGridCalculator);

    /**
     * Set the reference length scale
     *
     * @param referenceLengthScale the reference length scale
     */
    void SetReferenceLengthScale(QLength referenceLengthScale);

    /**
     * Return the index for the Daughter node.
     * This method can be used to move cells out of the way as necessary.
     *
     * @param pNewCell  The cell to new cell
     * @param pParentCell  The parent cell
     * @param rCellPopulation  The CA-based cell population
     * @return the node index for the daughter cell.
     */
    virtual unsigned CalculateDaughterNodeIndex(CellPtr pNewCell,
                                                CellPtr pParentCell,
                                                CaBasedCellPopulation<SPACE_DIM>& rCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(Owen11CaBasedDivisionRule)

#endif // OWEN11CABASEDDIVISIONRULE_HPP_
