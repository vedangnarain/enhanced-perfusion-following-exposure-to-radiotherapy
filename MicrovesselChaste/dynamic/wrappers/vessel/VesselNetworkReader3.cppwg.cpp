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
#include "VesselNetworkReader.hpp"

#include "PythonObjectConverters.hpp"
#include "VesselNetworkReader3.cppwg.hpp"

namespace py = pybind11;
PYBIND11_CVECTOR_TYPECASTER2();
PYBIND11_CVECTOR_TYPECASTER3();   
typedef VesselNetworkReader<3 > VesselNetworkReader3;
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

void register_VesselNetworkReader3_class(py::module &m){
py::class_<VesselNetworkReader3  , std::shared_ptr<VesselNetworkReader3 >   >(m, "VesselNetworkReader3")
        .def(py::init< >())
        .def_static(
            "Create", 
            (::std::shared_ptr<VesselNetworkReader<3> >(*)()) &VesselNetworkReader3::Create, 
            " "  )
        .def(
            "Read", 
            (::std::shared_ptr<VesselNetwork<3> >(VesselNetworkReader3::*)()) &VesselNetworkReader3::Read, 
            " "  )
        .def(
            "SetRadiusArrayName", 
            (void(VesselNetworkReader3::*)(::std::string const &)) &VesselNetworkReader3::SetRadiusArrayName, 
            " " , py::arg("rRadius") )
        .def(
            "SetMergeCoincidentPoints", 
            (void(VesselNetworkReader3::*)(bool)) &VesselNetworkReader3::SetMergeCoincidentPoints, 
            " " , py::arg("mergePoints") )
        .def(
            "SetTargetSegmentLength", 
            (void(VesselNetworkReader3::*)(::QLength)) &VesselNetworkReader3::SetTargetSegmentLength, 
            " " , py::arg("targetSegmentLength") )
        .def(
            "SetFileName", 
            (void(VesselNetworkReader3::*)(::std::string const &)) &VesselNetworkReader3::SetFileName, 
            " " , py::arg("rFileName") )
        .def(
            "SetReferenceLengthScale", 
            (void(VesselNetworkReader3::*)(::QLength)) &VesselNetworkReader3::SetReferenceLengthScale, 
            " " , py::arg("rReferenceLength") )
    ;
}
