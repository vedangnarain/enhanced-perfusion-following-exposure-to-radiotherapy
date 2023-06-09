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
#include "DensityMap.hpp"

#include "PythonObjectConverters.hpp"
#include "DensityMap2.cppwg.hpp"

namespace py = pybind11;
PYBIND11_CVECTOR_TYPECASTER2();
PYBIND11_CVECTOR_TYPECASTER3();   
typedef DensityMap<2 > DensityMap2;
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

void register_DensityMap2_class(py::module &m){
py::class_<DensityMap2  , std::shared_ptr<DensityMap2 >   >(m, "DensityMap2")
        .def(py::init< >())
        .def_static(
            "Create", 
            (::std::shared_ptr<DensityMap<2> >(*)()) &DensityMap2::Create, 
            " "  )
        .def(
            "GetVesselNetwork", 
            (::std::shared_ptr<VesselNetwork<2> >(DensityMap2::*)()) &DensityMap2::GetVesselNetwork, 
            " "  )
        .def(
            "GetSamplingGrid", 
            (::vtkSmartPointer<vtkUnstructuredGrid>(DensityMap2::*)(::vtkSmartPointer<vtkUnstructuredGrid>)) &DensityMap2::GetSamplingGrid, 
            " " , py::arg("pGrid") )
        .def(
            "GetSamplingGrid", 
            (::vtkSmartPointer<vtkUnstructuredGrid>(DensityMap2::*)(::std::shared_ptr<RegularGrid<2> >)) &DensityMap2::GetSamplingGrid, 
            " " , py::arg("pGrid") )
        .def(
            "GetGridCalculator", 
            (::std::shared_ptr<GridCalculator<2> >(DensityMap2::*)()) &DensityMap2::GetGridCalculator, 
            " "  )
        .def(
            "rGetVesselSurfaceAreaDensity", 
            (::std::vector<double, std::allocator<double> > const &(DensityMap2::*)(bool)) &DensityMap2::rGetVesselSurfaceAreaDensity, 
            " " , py::arg("update") = true , py::return_value_policy::reference_internal)
        .def(
            "rGetVesselLineDensity", 
            (::std::vector<double, std::allocator<double> >(DensityMap2::*)(bool)) &DensityMap2::rGetVesselLineDensity, 
            " " , py::arg("update") = true )
        .def_static(
            "GetVesselLineDensity", 
            (void(*)(::vtkSmartPointer<vtkUnstructuredGrid>, ::std::shared_ptr<VesselNetwork<2> >, ::QLength)) &DensityMap2::GetVesselLineDensity, 
            " " , py::arg("pGrid"), py::arg("pNetwork"), py::arg("reference_length") )
        .def_static(
            "GetVesselTipDensity", 
            (void(*)(::vtkSmartPointer<vtkUnstructuredGrid>, ::std::shared_ptr<VesselNetwork<2> >, ::QLength)) &DensityMap2::GetVesselTipDensity, 
            " " , py::arg("pGrid"), py::arg("pNetwork"), py::arg("reference_length") )
        .def(
            "rGetPerfusedVesselSurfaceAreaDensity", 
            (::std::vector<double, std::allocator<double> > const &(DensityMap2::*)(bool)) &DensityMap2::rGetPerfusedVesselSurfaceAreaDensity, 
            " " , py::arg("update") = true , py::return_value_policy::reference_internal)
        .def(
            "rGetPerfusedVesselLineDensity", 
            (::std::vector<double, std::allocator<double> > const &(DensityMap2::*)(bool)) &DensityMap2::rGetPerfusedVesselLineDensity, 
            " " , py::arg("update") = true , py::return_value_policy::reference_internal)
        .def(
            "rGetVesselTipDensity", 
            (::std::vector<double, std::allocator<double> > const &(DensityMap2::*)(bool)) &DensityMap2::rGetVesselTipDensity, 
            " " , py::arg("update") = true , py::return_value_policy::reference_internal)
        .def(
            "rGetVesselBranchDensity", 
            (::std::vector<double, std::allocator<double> > const &(DensityMap2::*)(bool)) &DensityMap2::rGetVesselBranchDensity, 
            " " , py::arg("update") = true , py::return_value_policy::reference_internal)
        .def(
            "rGetVesselQuantityDensity", 
            (::std::vector<double, std::allocator<double> > const &(DensityMap2::*)(::std::string const &, bool)) &DensityMap2::rGetVesselQuantityDensity, 
            " " , py::arg("rQuantity"), py::arg("update") = true , py::return_value_policy::reference_internal)
        .def(
            "rGetCellDensity", 
            (::std::vector<double, std::allocator<double> > const &(DensityMap2::*)(bool)) &DensityMap2::rGetCellDensity, 
            " " , py::arg("update") = true , py::return_value_policy::reference_internal)
        .def(
            "rGetCellDensity", 
            (::std::vector<double, std::allocator<double> > const &(DensityMap2::*)(::boost::shared_ptr<AbstractCellMutationState>, bool)) &DensityMap2::rGetCellDensity, 
            " " , py::arg("pMutationState"), py::arg("update") = true , py::return_value_policy::reference_internal)
        .def_static(
            "IsPointInCell", 
            (bool(*)(::vtkSmartPointer<vtkCellLocator>, ::boost::numeric::ublas::c_vector<double, 2>, unsigned int)) &DensityMap2::IsPointInCell, 
            " " , py::arg("pCellLocator"), py::arg("loc"), py::arg("index") )
        .def_static(
            "LengthOfLineInCell", 
            (double(*)(::vtkSmartPointer<vtkUnstructuredGrid>, ::boost::numeric::ublas::c_vector<double, 2>, ::boost::numeric::ublas::c_vector<double, 2>, unsigned int, bool, bool)) &DensityMap2::LengthOfLineInCell, 
            " " , py::arg("pSamplingGrid"), py::arg("loc1"), py::arg("loc2"), py::arg("index"), py::arg("loc1InCell"), py::arg("loc2InCell") )
        .def(
            "SetVesselNetwork", 
            (void(DensityMap2::*)(::std::shared_ptr<VesselNetwork<2> >)) &DensityMap2::SetVesselNetwork, 
            " " , py::arg("pNetwork") )
        .def(
            "SetCellPopulation", 
            (void(DensityMap2::*)(::AbstractCellPopulation<2, 2> &, ::QLength, ::QConcentration)) &DensityMap2::SetCellPopulation, 
            " " , py::arg("rCellPopulation"), py::arg("cellPopulationReferenceLength"), py::arg("cellPopulationReferenceConcentration") )
        .def(
            "SetGrid", 
            (void(DensityMap2::*)(::std::shared_ptr<AbstractDiscreteContinuumGrid<2, 2> >)) &DensityMap2::SetGrid, 
            " " , py::arg("pGrid") )
        .def(
            "SetGridCalculator", 
            (void(DensityMap2::*)(::std::shared_ptr<GridCalculator<2> >)) &DensityMap2::SetGridCalculator, 
            " " , py::arg("pGridCalculator") )
    ;
}
