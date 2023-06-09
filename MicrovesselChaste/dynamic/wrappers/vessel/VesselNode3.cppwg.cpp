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
#include "VesselNode.hpp"

#include "PythonObjectConverters.hpp"
#include "VesselNode3.cppwg.hpp"

namespace py = pybind11;
PYBIND11_CVECTOR_TYPECASTER2();
PYBIND11_CVECTOR_TYPECASTER3();   
typedef VesselNode<3 > VesselNode3;
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);
typedef ::std::map<std::basic_string<char>, double, std::less<std::basic_string<char> >, std::allocator<std::pair<const std::basic_string<char>, double> > > _std_map_lt_std_basic_string_lt_char_gt__double_std_less_lt_std_basic_string_lt_char_gt__gt__std_allocator_lt_std_pair_lt_conststd_basic_string_lt_char_gt__double_gt__gt__gt_;

class VesselNode3_Overloads : public VesselNode3{
    public:
    using VesselNode3::VesselNode;
    ::std::map<std::basic_string<char>, double, std::less<std::basic_string<char> >, std::allocator<std::pair<const std::basic_string<char>, double> > > GetOutputData() override {
        PYBIND11_OVERLOAD(
            _std_map_lt_std_basic_string_lt_char_gt__double_std_less_lt_std_basic_string_lt_char_gt__gt__std_allocator_lt_std_pair_lt_conststd_basic_string_lt_char_gt__double_gt__gt__gt_,
            VesselNode3,
            GetOutputData,
            );
    }

};
void register_VesselNode3_class(py::module &m){
py::class_<VesselNode3 , VesselNode3_Overloads , std::shared_ptr<VesselNode3 >  , AbstractVesselNetworkComponent<3>  >(m, "VesselNode3")
        .def(py::init<::QLength, ::QLength, ::QLength >(), py::arg("v1") = 0_m, py::arg("v2") = 0_m, py::arg("v3") = 0_m)
        .def(py::init<::Vertex<3> const & >(), py::arg("rLocation"))
        .def_static(
            "Create", 
            (::std::shared_ptr<VesselNode<3> >(*)(::QLength, ::QLength, ::QLength)) &VesselNode3::Create, 
            " " , py::arg("v1"), py::arg("v2") = 0_m, py::arg("v3") = 0_m )
        .def_static(
            "Create", 
            (::std::shared_ptr<VesselNode<3> >(*)(::Vertex<3> const &)) &VesselNode3::Create, 
            " " , py::arg("location") )
        .def_static(
            "Create", 
            (::std::shared_ptr<VesselNode<3> >(*)(::VesselNode<3> const &)) &VesselNode3::Create, 
            " " , py::arg("rExistingNode") )
        .def_static(
            "Create", 
            (::std::shared_ptr<VesselNode<3> >(*)(::std::shared_ptr<VesselNode<3> >)) &VesselNode3::Create, 
            " " , py::arg("pExistingNode") )
        .def(
            "GetComparisonId", 
            (unsigned int(VesselNode3::*)()) &VesselNode3::GetComparisonId, 
            " "  )
        .def(
            "GetDistance", 
            (::QLength(VesselNode3::*)(::Vertex<3> const &) const ) &VesselNode3::GetDistance, 
            " " , py::arg("rLocation") )
        .def(
            "GetFlowProperties", 
            (::std::shared_ptr<NodeFlowProperties<3> >(VesselNode3::*)() const ) &VesselNode3::GetFlowProperties, 
            " "  )
        .def(
            "GetChemicalProperties", 
            (::std::shared_ptr<NodeChemicalProperties<3> >(VesselNode3::*)() const ) &VesselNode3::GetChemicalProperties, 
            " "  )
        .def(
            "rGetLocation", 
            (::Vertex<3> const &(VesselNode3::*)() const ) &VesselNode3::rGetLocation, 
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetNumberOfSegments", 
            (unsigned int(VesselNode3::*)() const ) &VesselNode3::GetNumberOfSegments, 
            " "  )
        .def(
            "GetOutputData", 
            (::std::map<std::basic_string<char>, double, std::less<std::basic_string<char> >, std::allocator<std::pair<const std::basic_string<char>, double> > >(VesselNode3::*)()) &VesselNode3::GetOutputData, 
            " "  )
        .def(
            "GetSegment", 
            (::std::shared_ptr<VesselSegment<3> >(VesselNode3::*)(unsigned int) const ) &VesselNode3::GetSegment, 
            " " , py::arg("index") )
        .def(
            "GetSegments", 
            (::std::vector<std::shared_ptr<VesselSegment<3> >, std::allocator<std::shared_ptr<VesselSegment<3> > > >(VesselNode3::*)() const ) &VesselNode3::GetSegments, 
            " "  )
        .def(
            "GetGlobalIndex", 
            (unsigned int(VesselNode3::*)()) &VesselNode3::GetGlobalIndex, 
            " "  )
        .def(
            "GetLocalIndex", 
            (unsigned int(VesselNode3::*)()) &VesselNode3::GetLocalIndex, 
            " "  )
        .def(
            "GetOwnerRank", 
            (unsigned int(VesselNode3::*)()) &VesselNode3::GetOwnerRank, 
            " "  )
        .def(
            "IsHalo", 
            (bool(VesselNode3::*)()) &VesselNode3::IsHalo, 
            " "  )
        .def(
            "HasHalo", 
            (bool(VesselNode3::*)()) &VesselNode3::HasHalo, 
            " "  )
        .def(
            "GetOtherProcessorRank", 
            (unsigned int(VesselNode3::*)()) &VesselNode3::GetOtherProcessorRank, 
            " "  )
        .def(
            "GetOtherProcessorLocalIndex", 
            (unsigned int(VesselNode3::*)()) &VesselNode3::GetOtherProcessorLocalIndex, 
            " "  )
        .def(
            "IsAttachedTo", 
            (bool(VesselNode3::*)(::std::shared_ptr<VesselSegment<3> > const) const ) &VesselNode3::IsAttachedTo, 
            " " , py::arg("pSegment") )
        .def(
            "IsCoincident", 
            (bool(VesselNode3::*)(::Vertex<3> const &) const ) &VesselNode3::IsCoincident, 
            " " , py::arg("rLocation") )
        .def(
            "IsMigrating", 
            (bool(VesselNode3::*)() const ) &VesselNode3::IsMigrating, 
            " "  )
        .def(
            "SetComparisonId", 
            (void(VesselNode3::*)(unsigned int)) &VesselNode3::SetComparisonId, 
            " " , py::arg("id") )
        .def(
            "SetFlowProperties", 
            (void(VesselNode3::*)(::NodeFlowProperties<3> const &)) &VesselNode3::SetFlowProperties, 
            " " , py::arg("rFlowProperties") )
        .def(
            "SetChemicalProperties", 
            (void(VesselNode3::*)(::NodeChemicalProperties<3> const &)) &VesselNode3::SetChemicalProperties, 
            " " , py::arg("rChemicalProperties") )
        .def(
            "SetIsMigrating", 
            (void(VesselNode3::*)(bool)) &VesselNode3::SetIsMigrating, 
            " " , py::arg("isMigrating") )
        .def(
            "SetLocation", 
            (void(VesselNode3::*)(::Vertex<3> const &)) &VesselNode3::SetLocation, 
            " " , py::arg("rLocation") )
        .def(
            "SetLocation", 
            (void(VesselNode3::*)(::QLength, ::QLength, ::QLength)) &VesselNode3::SetLocation, 
            " " , py::arg("v1"), py::arg("v2") = 0_m, py::arg("v3") = 0_m )
        .def(
            "SetGlobalIndex", 
            (void(VesselNode3::*)(unsigned int)) &VesselNode3::SetGlobalIndex, 
            " " , py::arg("index") )
        .def(
            "SetLocalIndex", 
            (void(VesselNode3::*)(unsigned int)) &VesselNode3::SetLocalIndex, 
            " " , py::arg("index") )
        .def(
            "SetOwnerRank", 
            (void(VesselNode3::*)(unsigned int)) &VesselNode3::SetOwnerRank, 
            " " , py::arg("rank") )
        .def(
            "SetIsHalo", 
            (void(VesselNode3::*)(bool)) &VesselNode3::SetIsHalo, 
            " " , py::arg("isHalo") )
        .def(
            "SetHasHalo", 
            (void(VesselNode3::*)(bool)) &VesselNode3::SetHasHalo, 
            " " , py::arg("hasHalo") )
        .def(
            "SetOtherProcessorRank", 
            (void(VesselNode3::*)(unsigned int)) &VesselNode3::SetOtherProcessorRank, 
            " " , py::arg("otherRank") )
        .def(
            "SetOtherProcessorLocalIndex", 
            (void(VesselNode3::*)(unsigned int)) &VesselNode3::SetOtherProcessorLocalIndex, 
            " " , py::arg("otherIndex") )
    ;
}
