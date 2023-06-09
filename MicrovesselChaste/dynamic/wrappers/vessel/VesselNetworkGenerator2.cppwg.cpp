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
#include "VesselNetworkGenerator.hpp"

#include "PythonObjectConverters.hpp"
#include "VesselNetworkGenerator2.cppwg.hpp"

namespace py = pybind11;
PYBIND11_CVECTOR_TYPECASTER2();
PYBIND11_CVECTOR_TYPECASTER3();   
typedef VesselNetworkGenerator<2 > VesselNetworkGenerator2;
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

void register_VesselNetworkGenerator2_class(py::module &m){
py::class_<VesselNetworkGenerator2  , std::shared_ptr<VesselNetworkGenerator2 >   >(m, "VesselNetworkGenerator2")
        .def(py::init< >())
        .def(
            "GenerateParallelNetwork", 
            (::std::shared_ptr<VesselNetwork<2> >(VesselNetworkGenerator2::*)(::std::shared_ptr<Part<2> >, ::QPerArea, ::VesselDistribution::Value, ::QLength, bool, ::std::vector<std::shared_ptr<Vertex<2> >, std::allocator<std::shared_ptr<Vertex<2> > > >)) &VesselNetworkGenerator2::GenerateParallelNetwork, 
            " " , py::arg("domain"), py::arg("targetDensity"), py::arg("distrbutionType"), py::arg("exclusionDistance") = 0_m, py::arg("useBbox") = false, py::arg("seeds") = std::vector<std::shared_ptr<Vertex<2> > >() )
        .def(
            "GenerateHexagonalNetwork", 
            (::std::shared_ptr<VesselNetwork<2> >(VesselNetworkGenerator2::*)(::QLength, ::QLength, ::QLength, bool)) &VesselNetworkGenerator2::GenerateHexagonalNetwork, 
            " " , py::arg("width"), py::arg("height"), py::arg("vesselLength"), py::arg("fillDomain") = false )
        .def(
            "GenerateHexagonalUnit", 
            (::std::shared_ptr<VesselNetwork<2> >(VesselNetworkGenerator2::*)(::QLength)) &VesselNetworkGenerator2::GenerateHexagonalUnit, 
            " " , py::arg("vesselLength") )
        .def(
            "GenerateBifurcationUnit", 
            (::std::shared_ptr<VesselNetwork<2> >(VesselNetworkGenerator2::*)(::QLength, ::Vertex<2>)) &VesselNetworkGenerator2::GenerateBifurcationUnit, 
            " " , py::arg("vesselLength"), py::arg("startPosition") = Vertex<2>() )
        .def(
            "GenerateSingleVessel", 
            (::std::shared_ptr<VesselNetwork<2> >(VesselNetworkGenerator2::*)(::QLength, ::Vertex<2>, unsigned int, unsigned int)) &VesselNetworkGenerator2::GenerateSingleVessel, 
            " " , py::arg("vesselLength"), py::arg("startPosition") = Vertex<2>(), py::arg("divisions") = 0, py::arg("axis") = 2 )
        .def(
            "GenerateOvalNetwork", 
            (::std::shared_ptr<VesselNetwork<2> >(VesselNetworkGenerator2::*)(::QLength, unsigned int, double, double)) &VesselNetworkGenerator2::GenerateOvalNetwork, 
            " " , py::arg("scaleFactor"), py::arg("num_increments") = 40, py::arg("a_param") = 0.5, py::arg("b_param") = 1. )
        .def(
            "GenerateFromPart", 
            (::std::shared_ptr<VesselNetwork<2> >(VesselNetworkGenerator2::*)(::std::shared_ptr<Part<2> >)) &VesselNetworkGenerator2::GenerateFromPart, 
            " " , py::arg("pPart") )
        .def(
            "PatternUnitByTranslation", 
            (void(VesselNetworkGenerator2::*)(::std::shared_ptr<VesselNetwork<2> >, ::std::array<unsigned int, 2>)) &VesselNetworkGenerator2::PatternUnitByTranslation, 
            " " , py::arg("pInputUnit"), py::arg("numberOfUnits") )
        .def(
            "MapToSphere", 
            (void(VesselNetworkGenerator2::*)(::std::shared_ptr<VesselNetwork<2> >, ::QLength, ::QLength, double, double)) &VesselNetworkGenerator2::MapToSphere, 
            " " , py::arg("pInputUnit"), py::arg("radius"), py::arg("thickess"), py::arg("azimuthExtent"), py::arg("polarExtent") )
        .def(
            "SetReferenceLengthScale", 
            (void(VesselNetworkGenerator2::*)(::QLength)) &VesselNetworkGenerator2::SetReferenceLengthScale, 
            " " , py::arg("rReferenceLength") )
    ;
}
