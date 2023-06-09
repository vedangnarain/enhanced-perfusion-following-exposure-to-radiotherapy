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

#ifndef OWEN11PARAMETERS_HPP_
#define OWEN11PARAMETERS_HPP_

#include "SmartPointers.hpp"
#include "BaseParameterInstance.hpp"
#include "ParameterInstance.hpp"

/**
 * This struct stores parameter values used in the paper Owen et al. (2011)
 */
struct Owen11Parameters
{
    /**
     * Network Inlet Pressure: P_in
     */
    static const std::shared_ptr<ParameterInstance<QPressure> > mpInletPressure;

    /**
     * Network Outlet Pressure: P_out
     */
    static const std::shared_ptr<ParameterInstance<QPressure> > mpOutletPressure;

    /**
     * Plasma viscosity \mu_{plasma}
     */
    static const std::shared_ptr<ParameterInstance<QDynamicViscosity> > mpPlasmaViscosity;

    /**
     * Minimum cell cycle period T_{min}_{normal}
     */
    static const std::shared_ptr<ParameterInstance<QTime> > mpMinimumCellCyclePeriodNormal;

    /**
     * Minimum cell cycle period T_{min}_{cancer}
     */
    static const std::shared_ptr<ParameterInstance<QTime> > mpMinimumCellCyclePeriodCancer;

    /**
     * Time for death due to sustained quiescence T_{death}
     */
    static const std::shared_ptr<ParameterInstance<QTime> > mpTimeToDeathDueToQuiescence;

    /**
     * Oxygen partial pressure at half max cell cycle rate C_{\phi}_{normal}
     */
    static const std::shared_ptr<ParameterInstance<QPressure> > mpOxygenPartialPressureAtHalfMaxCycleRateNormal;

    /**
     * Oxygen partial pressure at half max cell cycle rate C_{\phi}_{cancer}
     */
    static const std::shared_ptr<ParameterInstance<QPressure> > mpOxygenPartialPressureAtHalfMaxCycleRateCancer;

    /**
     * Oxygen partial pressure at quiescence C^{enter}_{quiesc}
     */
    static const std::shared_ptr<ParameterInstance<QPressure> > mpOxygenPartialPressureAtQuiescence;

    /**
     * Oxygen partial pressure to leave quiescence C^{leave}_{quiesc}
     */
    static const std::shared_ptr<ParameterInstance<QPressure> > mpOxygenPartialPressureLeaveQuiescence;

    /**
     * Chemotactic sensitivity \chi
     */
    static const std::shared_ptr<ParameterInstance<QDiffusivityPerConcentration> > mpChemotacticSensitivity;

    /**
     * Maximum cell motility D_{normal}
     */
    static const std::shared_ptr<ParameterInstance<QDiffusivity>  > mpCellMotilityNormal;

    /**
     * Maximum cell motility D_{cancer}
     */
    static const std::shared_ptr<ParameterInstance<QDiffusivity>  > mpCellMotilityCancer;

    /**
     * Maximum cell motility D_{endo}
     */
    static const std::shared_ptr<ParameterInstance<QDiffusivity>  > mpCellMotilityEndothelial;

//    /**
//     * p53 concentration for apoptosis p53_{THR}
//     */
//    static const std::shared_ptr<ParameterInstance<QConcentration> > mpP53ConcentrationForApoptosis;

    /**
     * Vessel permeability to oxygen \psi_{c}
     */
    static const std::shared_ptr<ParameterInstance<QMembranePermeability> > mpVesselOxygenPermeability;

    /**
     * Vessel permeability to vegf \psi_{v}
     */
    static const std::shared_ptr<ParameterInstance<QMembranePermeability> > mpVesselVegfPermeability;

    /**
     * Cell oxygen consumption rate \k_c^{cell}
     */
    static const std::shared_ptr<ParameterInstance<QRate> > mpCellOxygenConsumptionRate;

    /**
     * Cell vegf secretion rate \k_v^{cell}
     */
    static const std::shared_ptr<ParameterInstance<QConcentrationFlowRate> > mpCellVegfSecretionRate;

    /**
     * Maximum rate of sprouting P^{max}_{sprout}
     */
    static const std::shared_ptr<ParameterInstance<QRate> > mpMaximumSproutingRate;

    /**
     * VEGF concentration at half maximal vessel sprouting probability V_{sprout}
     */
    static const std::shared_ptr<ParameterInstance<QConcentration> > mpVegfConventrationAtHalfMaxProbSprouting;

    /**
     * Sprouting exclusion radius (R_{ex})
     */
    static const std::shared_ptr<ParameterInstance<QLength>  > mpSproutingExclusionRadius;

    /**
     * Maximum vessel survivial time with low wall shear stress T_{prune}
     */
    static const std::shared_ptr<ParameterInstance<QTime> > mpMaxTimeWithLowWallShearStress;

    /**
     * Critical wall shear stress for vessel pruning \tau_{wall}
     */
    static const std::shared_ptr<ParameterInstance<QPressure> > mpCriticalWallShearStress;

    /**
     * Master step time \delta t
     */
    static const std::shared_ptr<ParameterInstance<QTime> > mpMasterStepTime;

    /**
     * Simulation duration T_{final}
     */
    static const std::shared_ptr<ParameterInstance<QTime> > mpSimulationDuration;

    /**
     * Lattice spacing \delta x
     */
    static const std::shared_ptr<ParameterInstance<QLength>  > mpLatticeSpacing;

    /**
     * Oxygen diffusivity D_{c}
     */
    static const std::shared_ptr<ParameterInstance<QDiffusivity>  > mpOxygenDiffusivity;

    /**
     * Vegf diffusivity D_{v}
     */
    static const std::shared_ptr<ParameterInstance<QDiffusivity>  > mpVegfDiffusivity;

    /**
     * Vegf decay rate \delta_{v}
     */
    static const std::shared_ptr<ParameterInstance<QRate> > mpVegfDecayRate;

    /**
     * Inflow haematocrit H_{in}
     */
    static const std::shared_ptr<ParameterInstance<QDimensionless> > mpInflowHaematocrit;

    /**
     * Intracellular p53 production rate constant k_7
     */
    static const std::shared_ptr<ParameterInstance<QRate> > mpP53ProductionRateConstant;

    /**
     * Max p53 degradation rate k_dash_7
     */
    static const std::shared_ptr<ParameterInstance<QRate> > mpP53MaxDegradationRate;

    /**
     * Tissue oxygen tension for half-max p53 degradation C_{p53}
     */
    static const std::shared_ptr<ParameterInstance<QPressure> > mpOxygenTensionForHalfMaxP53Degradation;

    /**
     * Basal VEGF production rate in cell k_8
     */
    static const std::shared_ptr<ParameterInstance<QRate> > mpCellVegfProductionRate;

    /**
     * Max VEGF production rate in cell k_8_dash
     */
    static const std::shared_ptr<ParameterInstance<QRate> > mpMaxCellVegfProductionRate;

    /**
     * Effect of P53 on vegf production k_8_dash_dash
     */
    static const std::shared_ptr<ParameterInstance<QRate> > mpP53EffectOnVegfProduction;

    /**
     * Effect of VEGF on VEGF production j_5
     */
    static const std::shared_ptr<ParameterInstance<QDimensionless> > mpVegfEffectOnVegfProduction;

    /**
     * Tissue oxygen tension for half-max vegf degradation C_{VEGF}
     */
    static const std::shared_ptr<ParameterInstance<QPressure> > mpOxygenTensionForHalfMaxVegfDegradation;

    /**
     * Vessel radius update timestep \epsilon_t
     */
    static const std::shared_ptr<ParameterInstance<QTime> > mpVesselRadiusUpdateTimestep;

    /**
     * Minimum possible radius R_{MIN}
     */
    static const std::shared_ptr<ParameterInstance<QLength>  > mpMinimumRadius;

    /**
     * Maximum possible radius R_{MAX}
     */
    static const std::shared_ptr<ParameterInstance<QLength>  > mpMaximumRadius;

    /**
     * Reference flow rate for metabolic stimulus Q_{ref}
     */
    static const std::shared_ptr<ParameterInstance<QFlowRate> > mpReferenceFlowRateForMetabolicStimulus;

    /**
     * Shrinking tendency k_s
     */
    static const std::shared_ptr<ParameterInstance<QRate> > mpShrinkingTendency;

    /**
     * Shrinking to intravascaulr pressure k_p
     */
    static const std::shared_ptr<ParameterInstance<QRate> > mpSensitivityToIntravascularPressure;

    /**
     * Basal metabolic stimulus k^0_m
     */
    static const std::shared_ptr<ParameterInstance<QRate> > mpBasalMetabolicStimulus;

    /**
     * Reference partial pressure of inlet haematocrit vessels C_{Ref}
     */
    static const std::shared_ptr<ParameterInstance<QPressure> > mpReferencePartialPressure;

};

#endif /*OWEN11PARAMETERS_HPP_*/
