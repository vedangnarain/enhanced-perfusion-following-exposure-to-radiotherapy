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
#include "AbstractCellMutationState.hpp"

#include "AbstractCellMutationState.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCellMutationState AbstractCellMutationState;
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

void register_AbstractCellMutationState_class(py::module &m){
py::class_<AbstractCellMutationState  , std::shared_ptr<AbstractCellMutationState >   >(m, "AbstractCellMutationState")
        .def(py::init<unsigned int >(), py::arg("colour"))
        .def(
            "GetColour", 
            (unsigned int(AbstractCellMutationState::*)() const ) &AbstractCellMutationState::GetColour, 
            " "  )
    ;
}
