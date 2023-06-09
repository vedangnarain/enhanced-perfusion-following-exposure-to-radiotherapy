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
#include "LatticeBasedSproutingRule.hpp"

#include "PythonObjectConverters.hpp"
#include "LatticeBasedSproutingRule3.cppwg.hpp"

namespace py = pybind11;
PYBIND11_CVECTOR_TYPECASTER2();
PYBIND11_CVECTOR_TYPECASTER3();   
typedef LatticeBasedSproutingRule<3 > LatticeBasedSproutingRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);
typedef ::std::vector<std::shared_ptr<VesselNode<3> >, std::allocator<std::shared_ptr<VesselNode<3> > > > _std_vector_lt_std_shared_ptr_lt_VesselNode_lt_3_gt__gt__std_allocator_lt_std_shared_ptr_lt_VesselNode_lt_3_gt__gt__gt__gt_;

class LatticeBasedSproutingRule3_Overloads : public LatticeBasedSproutingRule3{
    public:
    using LatticeBasedSproutingRule3::LatticeBasedSproutingRule;
    ::std::vector<std::shared_ptr<VesselNode<3> >, std::allocator<std::shared_ptr<VesselNode<3> > > > GetSprouts(::std::vector<std::shared_ptr<VesselNode<3> >, std::allocator<std::shared_ptr<VesselNode<3> > > > const & rNodes) override {
        PYBIND11_OVERLOAD(
            _std_vector_lt_std_shared_ptr_lt_VesselNode_lt_3_gt__gt__std_allocator_lt_std_shared_ptr_lt_VesselNode_lt_3_gt__gt__gt__gt_,
            LatticeBasedSproutingRule3,
            GetSprouts,
            rNodes);
    }
    void SetGridCalculator(::std::shared_ptr<GridCalculator<3> > pGrid) override {
        PYBIND11_OVERLOAD(
            void,
            LatticeBasedSproutingRule3,
            SetGridCalculator,
            pGrid);
    }

};
void register_LatticeBasedSproutingRule3_class(py::module &m){
py::class_<LatticeBasedSproutingRule3 , LatticeBasedSproutingRule3_Overloads , std::shared_ptr<LatticeBasedSproutingRule3 >  , AbstractSproutingRule<3>  >(m, "LatticeBasedSproutingRule3")
        .def(py::init< >())
        .def_static(
            "Create", 
            (::std::shared_ptr<LatticeBasedSproutingRule<3> >(*)()) &LatticeBasedSproutingRule3::Create, 
            " "  )
        .def(
            "GetSprouts", 
            (::std::vector<std::shared_ptr<VesselNode<3> >, std::allocator<std::shared_ptr<VesselNode<3> > > >(LatticeBasedSproutingRule3::*)(::std::vector<std::shared_ptr<VesselNode<3> >, std::allocator<std::shared_ptr<VesselNode<3> > > > const &)) &LatticeBasedSproutingRule3::GetSprouts, 
            " " , py::arg("rNodes") )
        .def(
            "SetGridCalculator", 
            (void(LatticeBasedSproutingRule3::*)(::std::shared_ptr<GridCalculator<3> >)) &LatticeBasedSproutingRule3::SetGridCalculator, 
            " " , py::arg("pGrid") )
    ;
}
