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
#include "CellStateDependentDiscreteSource.hpp"

#include "PythonObjectConverters.hpp"
#include "CellStateDependentDiscreteSource2.cppwg.hpp"

namespace py = pybind11;
PYBIND11_CVECTOR_TYPECASTER2();
PYBIND11_CVECTOR_TYPECASTER3();   
typedef CellStateDependentDiscreteSource<2 > CellStateDependentDiscreteSource2;
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);
typedef ::std::vector<RQuantity<std::ratio<0, 1>, std::ratio<-3, 1>, std::ratio<-1, 1>, std::ratio<1, 1>, std::ratio<0, 1> >, std::allocator<RQuantity<std::ratio<0, 1>, std::ratio<-3, 1>, std::ratio<-1, 1>, std::ratio<1, 1>, std::ratio<0, 1> > > > _std_vector_lt_RQuantity_lt_std_ratio_lt_0_1_gt__std_ratio_lt_neg3_1_gt__std_ratio_lt_neg1_1_gt__std_ratio_lt_1_1_gt__std_ratio_lt_0_1_gt__gt__std_allocator_lt_RQuantity_lt_std_ratio_lt_0_1_gt__std_ratio_lt_neg3_1_gt__std_ratio_lt_neg1_1_gt__std_ratio_lt_1_1_gt__std_ratio_lt_0_1_gt__gt__gt__gt_;
typedef ::std::vector<RQuantity<std::ratio<0, 1>, std::ratio<0, 1>, std::ratio<-1, 1>, std::ratio<0, 1>, std::ratio<0, 1> >, std::allocator<RQuantity<std::ratio<0, 1>, std::ratio<0, 1>, std::ratio<-1, 1>, std::ratio<0, 1>, std::ratio<0, 1> > > > _std_vector_lt_RQuantity_lt_std_ratio_lt_0_1_gt__std_ratio_lt_0_1_gt__std_ratio_lt_neg1_1_gt__std_ratio_lt_0_1_gt__std_ratio_lt_0_1_gt__gt__std_allocator_lt_RQuantity_lt_std_ratio_lt_0_1_gt__std_ratio_lt_0_1_gt__std_ratio_lt_neg1_1_gt__std_ratio_lt_0_1_gt__std_ratio_lt_0_1_gt__gt__gt__gt_;

class CellStateDependentDiscreteSource2_Overloads : public CellStateDependentDiscreteSource2{
    public:
    using CellStateDependentDiscreteSource2::CellStateDependentDiscreteSource;
    ::std::vector<RQuantity<std::ratio<0, 1>, std::ratio<-3, 1>, std::ratio<-1, 1>, std::ratio<1, 1>, std::ratio<0, 1> >, std::allocator<RQuantity<std::ratio<0, 1>, std::ratio<-3, 1>, std::ratio<-1, 1>, std::ratio<1, 1>, std::ratio<0, 1> > > > GetConstantInUValues() override {
        PYBIND11_OVERLOAD(
            _std_vector_lt_RQuantity_lt_std_ratio_lt_0_1_gt__std_ratio_lt_neg3_1_gt__std_ratio_lt_neg1_1_gt__std_ratio_lt_1_1_gt__std_ratio_lt_0_1_gt__gt__std_allocator_lt_RQuantity_lt_std_ratio_lt_0_1_gt__std_ratio_lt_neg3_1_gt__std_ratio_lt_neg1_1_gt__std_ratio_lt_1_1_gt__std_ratio_lt_0_1_gt__gt__gt__gt_,
            CellStateDependentDiscreteSource2,
            GetConstantInUValues,
            );
    }
    ::std::vector<RQuantity<std::ratio<0, 1>, std::ratio<0, 1>, std::ratio<-1, 1>, std::ratio<0, 1>, std::ratio<0, 1> >, std::allocator<RQuantity<std::ratio<0, 1>, std::ratio<0, 1>, std::ratio<-1, 1>, std::ratio<0, 1>, std::ratio<0, 1> > > > GetLinearInUValues() override {
        PYBIND11_OVERLOAD(
            _std_vector_lt_RQuantity_lt_std_ratio_lt_0_1_gt__std_ratio_lt_0_1_gt__std_ratio_lt_neg1_1_gt__std_ratio_lt_0_1_gt__std_ratio_lt_0_1_gt__gt__std_allocator_lt_RQuantity_lt_std_ratio_lt_0_1_gt__std_ratio_lt_0_1_gt__std_ratio_lt_neg1_1_gt__std_ratio_lt_0_1_gt__std_ratio_lt_0_1_gt__gt__gt__gt_,
            CellStateDependentDiscreteSource2,
            GetLinearInUValues,
            );
    }

};
void register_CellStateDependentDiscreteSource2_class(py::module &m){
py::class_<CellStateDependentDiscreteSource2 , CellStateDependentDiscreteSource2_Overloads , std::shared_ptr<CellStateDependentDiscreteSource2 >  , DiscreteSource<2>  >(m, "CellStateDependentDiscreteSource2")
        .def(py::init< >())
        .def_static(
            "Create", 
            (::std::shared_ptr<CellStateDependentDiscreteSource<2> >(*)()) &CellStateDependentDiscreteSource2::Create, 
            " "  )
        .def(
            "GetConstantInUValues", 
            (::std::vector<RQuantity<std::ratio<0, 1>, std::ratio<-3, 1>, std::ratio<-1, 1>, std::ratio<1, 1>, std::ratio<0, 1> >, std::allocator<RQuantity<std::ratio<0, 1>, std::ratio<-3, 1>, std::ratio<-1, 1>, std::ratio<1, 1>, std::ratio<0, 1> > > >(CellStateDependentDiscreteSource2::*)()) &CellStateDependentDiscreteSource2::GetConstantInUValues, 
            " "  )
        .def(
            "GetLinearInUValues", 
            (::std::vector<RQuantity<std::ratio<0, 1>, std::ratio<0, 1>, std::ratio<-1, 1>, std::ratio<0, 1>, std::ratio<0, 1> >, std::allocator<RQuantity<std::ratio<0, 1>, std::ratio<0, 1>, std::ratio<-1, 1>, std::ratio<0, 1>, std::ratio<0, 1> > > >(CellStateDependentDiscreteSource2::*)()) &CellStateDependentDiscreteSource2::GetLinearInUValues, 
            " "  )
        .def(
            "SetStateRateMap", 
            (void(CellStateDependentDiscreteSource2::*)(::std::map<unsigned int, RQuantity<std::ratio<0, 1>, std::ratio<-3, 1>, std::ratio<-1, 1>, std::ratio<1, 1>, std::ratio<0, 1> >, std::less<unsigned int>, std::allocator<std::pair<const unsigned int, RQuantity<std::ratio<0, 1>, std::ratio<-3, 1>, std::ratio<-1, 1>, std::ratio<1, 1>, std::ratio<0, 1> > > > >)) &CellStateDependentDiscreteSource2::SetStateRateMap, 
            " " , py::arg("stateRateMap") )
        .def(
            "SetStateRateThresholdMap", 
            (void(CellStateDependentDiscreteSource2::*)(::std::map<unsigned int, RQuantity<std::ratio<0, 1>, std::ratio<-3, 1>, std::ratio<0, 1>, std::ratio<1, 1>, std::ratio<0, 1> >, std::less<unsigned int>, std::allocator<std::pair<const unsigned int, RQuantity<std::ratio<0, 1>, std::ratio<-3, 1>, std::ratio<0, 1>, std::ratio<1, 1>, std::ratio<0, 1> > > > >)) &CellStateDependentDiscreteSource2::SetStateRateThresholdMap, 
            " " , py::arg("stateThresholdMap") )
    ;
}
