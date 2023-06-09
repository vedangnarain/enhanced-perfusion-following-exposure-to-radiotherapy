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
#include "ShrinkingStimulusCalculator.hpp"

#include "PythonObjectConverters.hpp"
#include "ShrinkingStimulusCalculator2.cppwg.hpp"

namespace py = pybind11;
PYBIND11_CVECTOR_TYPECASTER2();
PYBIND11_CVECTOR_TYPECASTER3();   
typedef ShrinkingStimulusCalculator<2 > ShrinkingStimulusCalculator2;
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

class ShrinkingStimulusCalculator2_Overloads : public ShrinkingStimulusCalculator2{
    public:
    using ShrinkingStimulusCalculator2::ShrinkingStimulusCalculator;
    void Calculate() override {
        PYBIND11_OVERLOAD(
            void,
            ShrinkingStimulusCalculator2,
            Calculate,
            );
    }

};
void register_ShrinkingStimulusCalculator2_class(py::module &m){
py::class_<ShrinkingStimulusCalculator2 , ShrinkingStimulusCalculator2_Overloads , std::shared_ptr<ShrinkingStimulusCalculator2 >  , AbstractVesselNetworkCalculator<2>  >(m, "ShrinkingStimulusCalculator2")
        .def(py::init< >())
        .def_static(
            "Create", 
            (::std::shared_ptr<ShrinkingStimulusCalculator<2> >(*)()) &ShrinkingStimulusCalculator2::Create, 
            " "  )
        .def(
            "GetStimulus", 
            (::QRate(ShrinkingStimulusCalculator2::*)()) &ShrinkingStimulusCalculator2::GetStimulus, 
            " "  )
        .def(
            "SetStimulus", 
            (void(ShrinkingStimulusCalculator2::*)(::QRate)) &ShrinkingStimulusCalculator2::SetStimulus, 
            " " , py::arg("stimulus") )
        .def(
            "Calculate", 
            (void(ShrinkingStimulusCalculator2::*)()) &ShrinkingStimulusCalculator2::Calculate, 
            " "  )
    ;
}
