#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <petsc/private/vecimpl.h>
#include <petsc/private/matimpl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "UnitCollection.hpp"
#include "vtkPolyData.h"
#include "AbstractFiniteDifferenceSolverBase.hpp"

#include "PythonObjectConverters.hpp"
#include "AbstractFiniteDifferenceSolverBase3.cppwg.hpp"

namespace py = pybind11;
PYBIND11_CVECTOR_TYPECASTER2();
PYBIND11_CVECTOR_TYPECASTER3();   
typedef AbstractFiniteDifferenceSolverBase<3 > AbstractFiniteDifferenceSolverBase3;
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);
PYBIND11_MAKE_OPAQUE(Vec);
PYBIND11_MAKE_OPAQUE(Mat);

class AbstractFiniteDifferenceSolverBase3_Overloads : public AbstractFiniteDifferenceSolverBase3{
    public:
    using AbstractFiniteDifferenceSolverBase3::AbstractFiniteDifferenceSolverBase;
    void AddDiscreteTermsToMatrix() override {
        PYBIND11_OVERLOAD(
            void,
            AbstractFiniteDifferenceSolverBase3,
            AddDiscreteTermsToMatrix,
            );
    }
    void AddDiscreteTermsToRhs() override {
        PYBIND11_OVERLOAD(
            void,
            AbstractFiniteDifferenceSolverBase3,
            AddDiscreteTermsToRhs,
            );
    }
    void AssembleMatrix() override {
        PYBIND11_OVERLOAD_PURE(
            void,
            AbstractFiniteDifferenceSolverBase3,
            AssembleMatrix,
            );
    }
    void AssembleVector() override {
        PYBIND11_OVERLOAD_PURE(
            void,
            AbstractFiniteDifferenceSolverBase3,
            AssembleVector,
            );
    }
    void Solve() override {
        PYBIND11_OVERLOAD(
            void,
            AbstractFiniteDifferenceSolverBase3,
            Solve,
            );
    }
    void Setup() override {
        PYBIND11_OVERLOAD(
            void,
            AbstractFiniteDifferenceSolverBase3,
            Setup,
            );
    }
    void Update() override {
        PYBIND11_OVERLOAD(
            void,
            AbstractFiniteDifferenceSolverBase3,
            Update,
            );
    }

};
void register_AbstractFiniteDifferenceSolverBase3_class(py::module &m){
py::class_<AbstractFiniteDifferenceSolverBase3 , AbstractFiniteDifferenceSolverBase3_Overloads , std::shared_ptr<AbstractFiniteDifferenceSolverBase3 >  , AbstractRegularGridDiscreteContinuumSolver<3>  >(m, "AbstractFiniteDifferenceSolverBase3")
        .def(
            "AddDiscreteTermsToMatrix", 
            (void(AbstractFiniteDifferenceSolverBase3::*)()) &AbstractFiniteDifferenceSolverBase3::AddDiscreteTermsToMatrix, 
            " "  )
        .def(
            "AddDiscreteTermsToRhs", 
            (void(AbstractFiniteDifferenceSolverBase3::*)()) &AbstractFiniteDifferenceSolverBase3::AddDiscreteTermsToRhs, 
            " "  )
        .def(
            "AssembleMatrix", 
            (void(AbstractFiniteDifferenceSolverBase3::*)()) &AbstractFiniteDifferenceSolverBase3::AssembleMatrix, 
            " "  )
        .def(
            "AssembleVector", 
            (void(AbstractFiniteDifferenceSolverBase3::*)()) &AbstractFiniteDifferenceSolverBase3::AssembleVector, 
            " "  )
        .def(
            "Solve", 
            (void(AbstractFiniteDifferenceSolverBase3::*)()) &AbstractFiniteDifferenceSolverBase3::Solve, 
            " "  )
        .def(
            "Setup", 
            (void(AbstractFiniteDifferenceSolverBase3::*)()) &AbstractFiniteDifferenceSolverBase3::Setup, 
            " "  )
        .def(
            "Update", 
            (void(AbstractFiniteDifferenceSolverBase3::*)()) &AbstractFiniteDifferenceSolverBase3::Update, 
            " "  )
        .def(
            "UpdateBoundaryConditionsEachSolve", 
            (void(AbstractFiniteDifferenceSolverBase3::*)(bool)) &AbstractFiniteDifferenceSolverBase3::UpdateBoundaryConditionsEachSolve, 
            " " , py::arg("doUpdate") )
    ;
}
