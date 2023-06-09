#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <petsc/private/vecimpl.h>
#include <petsc/private/matimpl.h>
#include <petsc/private/tsimpl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "UnitCollection.hpp"
#include "vtkPolyData.h"
#include "CoupledLumpedSystemFiniteElementSolver.hpp"

#include "PythonObjectConverters.hpp"
#include "CoupledLumpedSystemFiniteElementSolver2.cppwg.hpp"

namespace py = pybind11;
PYBIND11_CVECTOR_TYPECASTER2();
PYBIND11_CVECTOR_TYPECASTER3();   
typedef CoupledLumpedSystemFiniteElementSolver<2 > CoupledLumpedSystemFiniteElementSolver2;
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);
PYBIND11_MAKE_OPAQUE(Vec);
PYBIND11_MAKE_OPAQUE(Mat);
PYBIND11_MAKE_OPAQUE(TS);

class CoupledLumpedSystemFiniteElementSolver2_Overloads : public CoupledLumpedSystemFiniteElementSolver2{
    public:
    using CoupledLumpedSystemFiniteElementSolver2::CoupledLumpedSystemFiniteElementSolver;
    void Solve() override {
        PYBIND11_OVERLOAD(
            void,
            CoupledLumpedSystemFiniteElementSolver2,
            Solve,
            );
    }

};
void register_CoupledLumpedSystemFiniteElementSolver2_class(py::module &m){
py::class_<CoupledLumpedSystemFiniteElementSolver2 , CoupledLumpedSystemFiniteElementSolver2_Overloads , std::shared_ptr<CoupledLumpedSystemFiniteElementSolver2 >  , AbstractFiniteElementSolverBase<2>  >(m, "CoupledLumpedSystemFiniteElementSolver2")
        .def(py::init< >())
        .def_static(
            "Create", 
            (::std::shared_ptr<CoupledLumpedSystemFiniteElementSolver<2> >(*)()) &CoupledLumpedSystemFiniteElementSolver2::Create, 
            " "  )
        .def(
            "rGetIntermediateSolutions", 
            (::std::vector<std::pair<std::vector<double, std::allocator<double> >, double>, std::allocator<std::pair<std::vector<double, std::allocator<double> >, double> > > const &(CoupledLumpedSystemFiniteElementSolver2::*)()) &CoupledLumpedSystemFiniteElementSolver2::rGetIntermediateSolutions, 
            " "  , py::return_value_policy::reference_internal)
        .def(
            "SetTargetTimeIncrement", 
            (void(CoupledLumpedSystemFiniteElementSolver2::*)(double)) &CoupledLumpedSystemFiniteElementSolver2::SetTargetTimeIncrement, 
            " " , py::arg("targetIncrement") )
        .def(
            "SetStartTime", 
            (void(CoupledLumpedSystemFiniteElementSolver2::*)(double)) &CoupledLumpedSystemFiniteElementSolver2::SetStartTime, 
            " " , py::arg("startTime") )
        .def(
            "SetEndTime", 
            (void(CoupledLumpedSystemFiniteElementSolver2::*)(double)) &CoupledLumpedSystemFiniteElementSolver2::SetEndTime, 
            " " , py::arg("endTime") )
        .def(
            "SetInitialGuess", 
            (void(CoupledLumpedSystemFiniteElementSolver2::*)(::std::vector<double, std::allocator<double> > const &)) &CoupledLumpedSystemFiniteElementSolver2::SetInitialGuess, 
            " " , py::arg("rInitialGuess") )
        .def(
            "SetStoreIntermediateSolutions", 
            (void(CoupledLumpedSystemFiniteElementSolver2::*)(bool, unsigned int)) &CoupledLumpedSystemFiniteElementSolver2::SetStoreIntermediateSolutions, 
            " " , py::arg("store"), py::arg("frequency") = 1 )
        .def(
            "SetWriteIntermediateSolutions", 
            (void(CoupledLumpedSystemFiniteElementSolver2::*)(bool, unsigned int)) &CoupledLumpedSystemFiniteElementSolver2::SetWriteIntermediateSolutions, 
            " " , py::arg("write"), py::arg("frequency") = 1 )
        .def(
            "SetUseCoupling", 
            (void(CoupledLumpedSystemFiniteElementSolver2::*)(bool)) &CoupledLumpedSystemFiniteElementSolver2::SetUseCoupling, 
            " " , py::arg("useCoupling") )
        .def(
            "Solve", 
            (void(CoupledLumpedSystemFiniteElementSolver2::*)()) &CoupledLumpedSystemFiniteElementSolver2::Solve, 
            " "  )
    ;
}
