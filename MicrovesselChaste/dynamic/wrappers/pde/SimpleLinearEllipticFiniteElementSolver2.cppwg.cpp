#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "UnitCollection.hpp"
#include "vtkPolyData.h"
#include "SimpleLinearEllipticFiniteElementSolver.hpp"

#include "PythonObjectConverters.hpp"
#include "SimpleLinearEllipticFiniteElementSolver2.cppwg.hpp"

namespace py = pybind11;
PYBIND11_CVECTOR_TYPECASTER2();
PYBIND11_CVECTOR_TYPECASTER3();   
typedef SimpleLinearEllipticFiniteElementSolver<2 > SimpleLinearEllipticFiniteElementSolver2;
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

class SimpleLinearEllipticFiniteElementSolver2_Overloads : public SimpleLinearEllipticFiniteElementSolver2{
    public:
    using SimpleLinearEllipticFiniteElementSolver2::SimpleLinearEllipticFiniteElementSolver;
    void Solve() override {
        PYBIND11_OVERLOAD(
            void,
            SimpleLinearEllipticFiniteElementSolver2,
            Solve,
            );
    }

};
void register_SimpleLinearEllipticFiniteElementSolver2_class(py::module &m){
py::class_<SimpleLinearEllipticFiniteElementSolver2 , SimpleLinearEllipticFiniteElementSolver2_Overloads , std::shared_ptr<SimpleLinearEllipticFiniteElementSolver2 >  , AbstractFiniteElementSolverBase<2>  >(m, "SimpleLinearEllipticFiniteElementSolver2")
        .def(py::init< >())
        .def_static(
            "Create", 
            (::std::shared_ptr<SimpleLinearEllipticFiniteElementSolver<2> >(*)()) &SimpleLinearEllipticFiniteElementSolver2::Create, 
            " "  )
        .def(
            "Solve", 
            (void(SimpleLinearEllipticFiniteElementSolver2::*)()) &SimpleLinearEllipticFiniteElementSolver2::Solve, 
            " "  )
    ;
}
