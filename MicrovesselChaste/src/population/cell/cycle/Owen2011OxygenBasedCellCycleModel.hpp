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

#ifndef OWEN2011OXYGENBASEDCELLCYCLEMODEL_HPP_
#define OWEN2011OXYGENBASEDCELLCYCLEMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <vector>
#include "AbstractOdeBasedPhaseBasedCellCycleModel.hpp"
#include "AbstractCellMutationState.hpp"
#include "CancerCellMutationState.hpp"
#include "Owen2011OxygenBasedCellCycleOdeSystem.hpp"
#include "UnitCollection.hpp"

/**
 * Oxygen-dependent ODE-based cell-cycle model. Published by Owen et al. 2011
 */
class Owen2011OxygenBasedCellCycleModel : public AbstractOdeBasedPhaseBasedCellCycleModel
{

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model and ODE system.
     *
     * @param archive the archive
     * @param version the archive version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeBasedPhaseBasedCellCycleModel>(*this);
        archive & mCurrentQuiescentDuration;
        archive & mCurrentQuiescenceOnsetTime;
        archive & mEnterQuiescenceOxygenConcentration;
        archive & mLeaveQuiescenceOxygenConcentration;
        archive & mCriticalQuiescentDuration;
    }

    /**
     * Adjust any ODE parameters needed before solving until currentTime.
     *
     * @param currentTime  the time up to which the system will be solved.
     */
    void AdjustOdeParameters(double currentTime);

protected:

    /**
     * ODE integration time-step
     */
    QTime mOdeIntegrationTimeStep;

    /**
     * Cell cycle fraction at which S phase begins
     */
    QDimensionless sOnset;

    /**
     * Cell cycle fraction at which G2 phase begins
     */
    QDimensionless g2Onset;

    /**
     * Cell cycle fraction at which M phase begins
     */
    QDimensionless mOnset;

    /**
     * The reference time scale
     */
    QTime mReferenceTimeScale;

    /**
     * The reference concentration scale
     */
    QConcentration mReferenceConcentrationScale;

    /**
     * The reference solubility
     */
    QSolubility mReferenceSolubility;

    /**
     * Maximum initial value allocated to phi.
     */
    QDimensionless mMaxRandInitialPhase;

    /**
     * How long the current period of quiescence has lasted.
     */
    QTime mCurrentQuiescentDuration;

    /**
     * The time when the current period of quiescence began.
     */
    QTime mCurrentQuiescenceOnsetTime;

    /**
     * Oxygen concentration below which cancerous cells enter quiescence.
     * A prolonged period of quiescence causes the cell to become apoptotic.
     */
    QPressure mEnterQuiescenceOxygenConcentration;

    /**
     * Oxygen concentration above which cancerous cells leave their state of being quiescent
     */
    QPressure mLeaveQuiescenceOxygenConcentration;

    /**
     * Critical quiescent duration.
     */
    QTime mCriticalQuiescentDuration;

    /**
     * p53 threshold above which normal cells become apoptotic in a healthy environment.
     */
    QDimensionless mp53ThresholdForApoptosisOfNormalCellsInHealthyMicroenvironment;

    /**
     * p53 threshold above which normal cells become apoptotic in a tumour environment.
     */
    QDimensionless mp53ThresholdForApoptosisOfNormalCellsInTumourMicroenvironment;

    /**
     * Threshold which defines the proportion of neighbours which must be normal in order for
     * a cell's micro-environment to be considered normal.
     */
    QDimensionless mthresholdFractionOfNormalCellNeighbours;

public:

    /**
     * Default constructor.
     *
     * @param pOdeSolver An optional pointer to a cell-cycle model ODE solver object (allows the use of different ODE solvers)
     */
    Owen2011OxygenBasedCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     * @return a cell cycle model
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Check if the oxygen concentration of the cell is below the EnterQuiescenceOxygenConcentration.
     * If it is true the label cells.
     */
    void CheckAndLabelCell();

    /**
     * @return mCurrentQuiescentDuration
     */
    QTime GetCurrentQuiescentDuration();

    /**
     * @return mCurrentQuiescenceOnsetTime
     */
    QTime GetCurrentQuiescenceOnsetTime();

    /**
     * @return mEnterQuiescenceOxygenConcentration
     */
    QPressure GetEnterQuiescenceOxygenConcentration();

    /**
     * @return mCriticalQuiescentDuration
     */
    QTime GetCriticalQuiescentDuration();

    /**
     * Get the duration of the cell's S phase.
     * @return the duration of the cell's S phase.
     */
    double GetSDuration() const;

    /**
     * Get the duration of the cell's G2 phase.
     * @return the duration of the cell's G2 phase.
     */
    double GetG2Duration() const;

    /**
     * Get the duration of the cell's M phase.
     * @return the duration of the cell's M phase.
     */
    double GetMDuration() const;

    /**
     * Get the value of phi.
     * @return the value of phi.
     */
    double GetPhi();

    /**
     * Get the value of VEGF.
     * @return the value of VEGF.
     */
    double GetVEGF();

    /**
     * Get the value of p53.
     * @return the value of p53.
     */
    double GetP53();

    /**
     * @return mLeaveQuiescenceOxygenConcentration
     */
    QPressure GetLeaveQuiescenceOxygenConcentration();

    /**
     * Initialise the cell-cycle model at the start of a simulation.
     *
     * This method will be called precisely once per cell set up in the initial
     * cell population. It is not called on cell division; use ResetForDivision(),
     * CreateCellCycleModel() and InitialiseDaughterCell() for that.
     *
     * By the time this is called, a CellPopulation will have been set up, so the model
     * can know where its cell is located in space. If relevant to the simulation,
     * any singletons will also have been initialised.
     */
    void Initialise();

    /**
     * Overridden InitialiseDaughterCell() method.
     */
    void InitialiseDaughterCell();

    /**
     * Outputs cell cycle model parameters to files.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);

    /**
     * @return whether the cell is ready to divide (enter M phase).
     *
     * The intention is that this method is called precisely once at
     * each timestep of the simulation. However this does not appear
     * to always be the case at present, and so it can cope with more
     * unusual usage patterns.
     */
    bool ReadyToDivide();

    /**
     * Resets the oxygen-based model to the start of the cell cycle
     * (this model does not cycle naturally). Cells are given a new
     * birth time and cell cycle proteins are reset. Note that the
     * oxygen concentration maintains its current value.
     *
     * Should only be called by the Cell Divide() method.
     */
    virtual void ResetForDivision();

    /**
     * Set maximum phase of cell upon initialisation.
     * @param rand_max_phase the maximum phase value at initialization
     */
    void SetMaxRandInitialPhase(QDimensionless rand_max_phase);

    /**
     * Set method for mEnterQuiescenceOxygenConcentration.
     *
     * @param enterQuiescenceOxygenConcentration the new value of mEnterQuiescenceOxygenConcentration
     */
    void SetEnterQuiescenceOxygenConcentration(QPressure enterQuiescenceOxygenConcentration);

    /**
     * Set method for mLeaveQuiescenceOxygenConcentration.
     *
     * @param leaveQuiescenceOxygenConcentration the new value of mLeaveQuiescenceOxygenConcentration
     */
    void SetLeaveQuiescenceOxygenConcentration(QPressure leaveQuiescenceOxygenConcentration);

    /**
     * Set method for mCriticalQuiescentDuration.
     *
     * @param criticalQuiescentDuration the new value of mCriticalQuiescentDuration
     */
    void SetCriticalQuiescentDuration(QTime criticalQuiescentDuration);

    /**
     * Set method for mCurrentQuiescenceOnsetTime.
     *
     * @param currentQuiescenceOnsetTime the new value of mCurrentQuiescenceOnsetTime
     */
    void SetCurrentQuiescenceOnsetTime(QTime currentQuiescenceOnsetTime);

    /**
     * Set the G2 onset time
     *
     * @param value the G2 onset time
     */
    void SetG2Onset(QDimensionless value);

    /**
     * Set the S2 onset time
     *
     * @param value the S2 onset time
     */
    void SetSOnset(QDimensionless value);

    /**
     * Set the M onset time
     *
     * @param value the M onset time
     */
    void SetMOnset(QDimensionless value);

    /**
     * Set the ODE solver time step
     *
     * @param timeStep the ODE solver time step
     */
    void SetOdeSolverTimeStep(QTime timeStep);

    /**
     * Set the reference timescale
     *
     * @param referenceTimeScale the reference timescale
     */
    void SetReferenceTimeScale(QTime referenceTimeScale);

    /**
     * Set the reference concentration scale
     *
     * @param referenceConcentrationScale the reference concentration scale
     */
    void SetReferenceConcentrationScale(QConcentration referenceConcentrationScale);

    /**
     * Set method for mthresholdFractionOfNormalCellNeighbours.
     *
     * @param value the new value of mthresholdFractionOfNormalCellNeighbours
     */
    void SetThresholdFractionOfNormalCellNeighbours(double value);

    /**
     * Update the duration for which the cell has been quiescent.
     */
    void UpdateQuiescentDuration();

    /**
     * Update cell-cycle phase.
     */
    void UpdateCellCyclePhase();
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Owen2011OxygenBasedCellCycleModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Owen2011OxygenBasedCellCycleModel)

#endif /*OWEN2011OXYGENBASEDCELLCYCLEMODEL_HPP_*/
