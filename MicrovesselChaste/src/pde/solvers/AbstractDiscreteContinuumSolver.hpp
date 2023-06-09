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

#ifndef ABSTRACTDISCRETECONTINUUMSOLVER_HPP_
#define ABSTRACTDISCRETECONTINUUMSOLVER_HPP_

#include <vector>
#include <string>
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning for now (gcc4.3)
#include <vtkSmartPointer.h>
#include "OutputFileHandler.hpp"
#include "DiscreteContinuumBoundaryCondition.hpp"
#include "RegularGrid.hpp"
#include "DiscreteContinuumMesh.hpp"
#include "UnitCollection.hpp"
#include "DensityMap.hpp"
#include "AbstractDiscreteContinuumGrid.hpp"

/**
 * Forward declare VTK members
 */
class vtkDataSet;

/**
 * Forward declaration as some PDEs can have DiscreteContinuumSolvers in their DiscreteSources
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class AbstractDiscreteContinuumPde;

/**
 * An abstract solver class for continuum-discrete field problems.
 * The class is used by the MicrovesselSolver to provide a concentration or dimensionless field for a single,
 * labelled quantity for cells and/or vessels. It contains methods for sampling on structured and unstructured grids.
 * Derived classes are responsible for updating the values of data fields in cells
 * and vessels on each call and optionally writing the solution to file.
 */
template<unsigned DIM>
class AbstractDiscreteContinuumSolver
{

protected:

    /**
     * File handler containing the output directory
     */
    std::shared_ptr<OutputFileHandler> mpOutputFileHandler;

    /**
     *  The filename for output
     */
    std::string mFilename;

    /**
     *  The label for the quantity being solved for
     */
    std::string mLabel;

    /**
     *  Has the Setup function been called.
     */
    bool IsSetupForSolve;

    /**
     *  Should the solution be written to file
     */
    bool mWriteSolution;

    /**
     * The PDE to be solved
     */
    std::shared_ptr<AbstractDiscreteContinuumPde<DIM, DIM> > mpPde;

    /**
     * The DiscreteContinuum boundary conditions, optional
     */
    std::vector<std::shared_ptr<DiscreteContinuumBoundaryCondition<DIM> > > mBoundaryConditions;

    /**
     * This is used internally to scale concentrations before and after linear system solves, reads and writes.
     * Since those functions don't use Boost Units. It should not affect the solution, but can be judiciously chosen
     * to avoid precision problems.
     */
    QConcentration mReferenceConcentration;

    /**
     * A solution field. Ordering is decided in child classes.
     */
    std::vector<double> mSolution;

    /**
     * A solution field. Ordering is decided in child classes.
     */
    std::vector<QConcentration > mConcentrations;

    /**
     * A density map
     */
    std::shared_ptr<DensityMap<DIM> > mpDensityMap;

public:

    /**
     * Constructor
     */
    AbstractDiscreteContinuumSolver();

    /**
     * Destructor
     */
    virtual ~AbstractDiscreteContinuumSolver();

    /**
     * Add a DiscreteContinuum boundary condition for the domain
     * @param pBoundaryCondition the boundary condition
     */
    void AddBoundaryCondition(std::shared_ptr<DiscreteContinuumBoundaryCondition<DIM> > pBoundaryCondition);

    /**
     * Return the value of the field with ordering determined by child classes
     * @return the value of the field with ordering determined by child classes
     */
    virtual std::vector<QConcentration > GetConcentrations();

    /**
     * Return the value of the field at all points on the supplied grid
     * @param pGrid the sampling grid
     * @return the value of the field ordered according to grid order
     */
    virtual std::vector<QConcentration > GetConcentrations(std::shared_ptr<AbstractDiscreteContinuumGrid<DIM> > pGrid);

    /**
     * Return the value of the field at the requested points
     * @param pSamplePoints vtk sample points
     * @return the value of the field ordered according to input point order
     */
    virtual std::vector<QConcentration > GetConcentrations(vtkSmartPointer<vtkPoints> pSamplePoints);

    /**
     * Return the density map
     * @return the density map
     */
    std::shared_ptr<DensityMap<DIM> > GetDensityMap();

    /**
     * Return the name of the field being solved for
     * @return a reference to the field name
     */
    const std::string& GetLabel();

    /**
     * Return the PDE
     * @return the DiscreteContinuum linear elliptic pde
     */
    std::shared_ptr<AbstractDiscreteContinuumPde<DIM, DIM> > GetPde();

    /**
     * Return the reference length scale, needs a grid first
     * @return the reference length value
     */
    QLength GetReferenceLength();

    /**
     * Return the reference concentration value.
     * @return the reference concentration value
     */
    QConcentration GetReferenceConcentration();

    /**
     * Return the value of the field with ordering determined by child classes
     * @return the value of the field with ordering determined by child classes
     */
    virtual std::vector<double> GetSolution();

    /**
     * Return the value of the field at the requested points
     * @param pSamplePoints vtk sample points
     * @return the value of the field ordered according to input point order
     */
    virtual std::vector<double> GetSolution(vtkSmartPointer<vtkPoints> pSamplePoints);

    /**
     * Return the gradients of the field at the requested points
     * @param pSamplePoints vtk sample points
     * @return the gradient of the field ordered according to input point order
     */
    virtual std::vector<c_vector<double, 3 > > GetSolutionGradients(vtkSmartPointer<vtkPoints> pSamplePoints);

    /**
     * Return the value of the field at the requested points. Used by the Python wrapper only.
     * @param pSamplePoints vtk sample points
     * @return the value of the field ordered according to input point order
     */
    virtual std::vector<double> GetSolutionP(vtkPoints* pSamplePoints);

    /**
     * Return the value of the field at all points on the supplied grid
     * @param pGrid the grid to be sampled
     * @return the value of the field ordered according to input point order
     */
    virtual std::vector<double> GetSolution(std::shared_ptr<AbstractDiscreteContinuumGrid<DIM> > pGrid);

    /**
     * Return the solution as vtk data
     * @return the vtk solution
     */
    virtual vtkSmartPointer<vtkDataSet> GetVtkSolution();

    /**
     * Set the file handler containing the working directory
     * @param pOutputFileHandler the file handler containing the working directory
     */
    void SetFileHandler(std::shared_ptr<OutputFileHandler> pOutputFileHandler);

    /**
     * Set the file name for output
     * @param rFilename the file name
     */
    void SetFileName(const std::string& rFilename);

    /**
     * Set the name of the field being solved for
     * @param rLabel a reference to the field name
     */
    void SetLabel(const std::string& rLabel);

    /**
     *  Set the PDE to be solved
     * @param pPde the pde to be solved
     */
    void SetPde(std::shared_ptr<AbstractDiscreteContinuumPde<DIM, DIM> > pPde);

    /**
     * Operations to be performed prior to the first solve
     */
    virtual void Setup() = 0;

    /**
     * Set the reference concentration
     * @param referenceConcentration the reference concentration
     */
    void SetReferenceConcentration(QConcentration referenceConcentration);

    /**
     * Set whether to write the solution to file on next solve
     * @param write write the solution
     */
    void SetWriteSolution(bool write=true);

    /**
     * Set the regular grid. Only do this if not working with discrete sinks or sources. Otherwise
     * supply a density map.
     * @param pGrid the grid
     */
    void SetGrid(std::shared_ptr<AbstractDiscreteContinuumGrid<DIM> > pGrid);

    /**
     * Set the density map
     */
    void SetDensityMap(std::shared_ptr<DensityMap<DIM> > pDensityMap);

    /**
     * Set the vessel network
     * @param pNetwork the vessel network
     */
    void SetVesselNetwork(std::shared_ptr<VesselNetwork<DIM> > pNetwork);

    /**
     * Set the cell population
     * @param rCellPopulation a reference to the cell population
     * @param cellPopulationReferenceLength the length scale for the cell population
     * @param cellPopulationReferenceConcentration the concentration scale for the cell population
     */
    void SetCellPopulation(AbstractCellPopulation<DIM>& rCellPopulation,
                           QLength cellPopulationReferenceLength,
                           QConcentration cellPopulationReferenceConcentration);

    /**
     * Do the solve
     */
    virtual void Solve() = 0;

    /**
     * Operations to be performed prior to every solve
     */
    virtual void Update() = 0;

    /**
     * Set the cell data to the values in the solution field
     */
    virtual void UpdateCellData() = 0;

    /**
     * Update the solution manually
     * @param rData solution data map
     */
    virtual void UpdateSolution(const std::vector<double>& rData);

    /**
     * Update the solution manually
     * @param rData solution data map
     */
    virtual void UpdateSolution(const std::vector<QConcentration >& rData);

    /**
     * Write the solution to file
     */
    virtual void Write() = 0;
};

#endif /* ABSTRACTDISCRETECONTINUUMSOLVER_HPP_ */
