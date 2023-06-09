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
#include "WallShearStressCalculator.hpp"

#include "PythonObjectConverters.hpp"
#include "WallShearStressCalculator3.cppwg.hpp"

namespace py = pybind11;
PYBIND11_CVECTOR_TYPECASTER2();
PYBIND11_CVECTOR_TYPECASTER3();   
typedef WallShearStressCalculator<3 > WallShearStressCalculator3;
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

class WallShearStressCalculator3_Overloads : public WallShearStressCalculator3{
    public:
    using WallShearStressCalculator3::WallShearStressCalculator;
    void Calculate() override {
        PYBIND11_OVERLOAD(
            void,
            WallShearStressCalculator3,
            Calculate,
            );
    }

};
void register_WallShearStressCalculator3_class(py::module &m){
py::class_<WallShearStressCalculator3 , WallShearStressCalculator3_Overloads , std::shared_ptr<WallShearStressCalculator3 >  , AbstractVesselNetworkCalculator<3>  >(m, "WallShearStressCalculator3")
        .def(py::init< >())
        .def_static(
            "Create", 
            (::std::shared_ptr<WallShearStressCalculator<3> >(*)()) &WallShearStressCalculator3::Create, 
            " "  )
        .def(
            "Calculate", 
            (void(WallShearStressCalculator3::*)()) &WallShearStressCalculator3::Calculate, 
            " "  )
    ;
}
