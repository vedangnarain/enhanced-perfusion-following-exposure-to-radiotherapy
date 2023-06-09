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
#include "AbstractVesselNetworkComponentChemicalProperties.hpp"

#include "PythonObjectConverters.hpp"
#include "AbstractVesselNetworkComponentChemicalProperties2.cppwg.hpp"

namespace py = pybind11;
PYBIND11_CVECTOR_TYPECASTER2();
PYBIND11_CVECTOR_TYPECASTER3();   
typedef AbstractVesselNetworkComponentChemicalProperties<2 > AbstractVesselNetworkComponentChemicalProperties2;
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

class AbstractVesselNetworkComponentChemicalProperties2_Overloads : public AbstractVesselNetworkComponentChemicalProperties2{
    public:
    using AbstractVesselNetworkComponentChemicalProperties2::AbstractVesselNetworkComponentChemicalProperties;
    void SetVegfConcentration(::QConcentration concentration) override {
        PYBIND11_OVERLOAD(
            void,
            AbstractVesselNetworkComponentChemicalProperties2,
            SetVegfConcentration,
            concentration);
    }
    void SetPermeability(::QMembranePermeability permeability) override {
        PYBIND11_OVERLOAD(
            void,
            AbstractVesselNetworkComponentChemicalProperties2,
            SetPermeability,
            permeability);
    }

};
void register_AbstractVesselNetworkComponentChemicalProperties2_class(py::module &m){
py::class_<AbstractVesselNetworkComponentChemicalProperties2 , AbstractVesselNetworkComponentChemicalProperties2_Overloads , std::shared_ptr<AbstractVesselNetworkComponentChemicalProperties2 >  , AbstractVesselNetworkComponentProperties<2>  >(m, "AbstractVesselNetworkComponentChemicalProperties2")
        .def(
            "GetVegfConcentration", 
            (::QConcentration(AbstractVesselNetworkComponentChemicalProperties2::*)() const ) &AbstractVesselNetworkComponentChemicalProperties2::GetVegfConcentration, 
            " "  )
        .def(
            "SetVegfConcentration", 
            (void(AbstractVesselNetworkComponentChemicalProperties2::*)(::QConcentration)) &AbstractVesselNetworkComponentChemicalProperties2::SetVegfConcentration, 
            " " , py::arg("concentration") )
        .def(
            "GetPermeability", 
            (::QMembranePermeability(AbstractVesselNetworkComponentChemicalProperties2::*)() const ) &AbstractVesselNetworkComponentChemicalProperties2::GetPermeability, 
            " "  )
        .def(
            "SetPermeability", 
            (void(AbstractVesselNetworkComponentChemicalProperties2::*)(::QMembranePermeability)) &AbstractVesselNetworkComponentChemicalProperties2::SetPermeability, 
            " " , py::arg("permeability") )
    ;
}
