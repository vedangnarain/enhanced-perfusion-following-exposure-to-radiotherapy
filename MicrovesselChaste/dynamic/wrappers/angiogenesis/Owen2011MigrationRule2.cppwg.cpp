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
#include "Owen2011MigrationRule.hpp"

#include "PythonObjectConverters.hpp"
#include "Owen2011MigrationRule2.cppwg.hpp"

namespace py = pybind11;
PYBIND11_CVECTOR_TYPECASTER2();
PYBIND11_CVECTOR_TYPECASTER3();   
typedef Owen2011MigrationRule<2 > Owen2011MigrationRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);
typedef ::std::vector<int, std::allocator<int> > _std_vector_lt_int_std_allocator_lt_int_gt__gt_;
typedef ::std::vector<double, std::allocator<double> > _std_vector_lt_double_std_allocator_lt_double_gt__gt_;

class Owen2011MigrationRule2_Overloads : public Owen2011MigrationRule2{
    public:
    using Owen2011MigrationRule2::Owen2011MigrationRule;
    ::std::vector<int, std::allocator<int> > GetIndices(::std::vector<std::shared_ptr<VesselNode<2> >, std::allocator<std::shared_ptr<VesselNode<2> > > > const & rNodes) override {
        PYBIND11_OVERLOAD(
            _std_vector_lt_int_std_allocator_lt_int_gt__gt_,
            Owen2011MigrationRule2,
            GetIndices,
            rNodes);
    }

};
void register_Owen2011MigrationRule2_class(py::module &m){
py::class_<Owen2011MigrationRule2 , Owen2011MigrationRule2_Overloads , std::shared_ptr<Owen2011MigrationRule2 >  , LatticeBasedMigrationRule<2>  >(m, "Owen2011MigrationRule2")
        .def(py::init< >())
        .def_static(
            "Create", 
            (::std::shared_ptr<Owen2011MigrationRule<2> >(*)()) &Owen2011MigrationRule2::Create, 
            " "  )
        .def(
            "GetIndices", 
            (::std::vector<int, std::allocator<int> >(Owen2011MigrationRule2::*)(::std::vector<std::shared_ptr<VesselNode<2> >, std::allocator<std::shared_ptr<VesselNode<2> > > > const &)) &Owen2011MigrationRule2::GetIndices, 
            " " , py::arg("rNodes") )
        .def(
            "SetCellChemotacticParameter", 
            (void(Owen2011MigrationRule2::*)(::QDiffusivityPerConcentration)) &Owen2011MigrationRule2::SetCellChemotacticParameter, 
            " " , py::arg("cellChemotacticParameter") )
        .def(
            "SetCellMotilityParameter", 
            (void(Owen2011MigrationRule2::*)(::QDiffusivity)) &Owen2011MigrationRule2::SetCellMotilityParameter, 
            " " , py::arg("cellMotility") )
    ;
}
