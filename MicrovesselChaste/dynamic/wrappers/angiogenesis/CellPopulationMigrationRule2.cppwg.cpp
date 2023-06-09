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
#include "CellPopulationMigrationRule.hpp"

#include "PythonObjectConverters.hpp"
#include "CellPopulationMigrationRule2.cppwg.hpp"

namespace py = pybind11;
PYBIND11_CVECTOR_TYPECASTER2();
PYBIND11_CVECTOR_TYPECASTER3();   
typedef CellPopulationMigrationRule<2 > CellPopulationMigrationRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);
typedef ::std::vector<int, std::allocator<int> > _std_vector_lt_int_std_allocator_lt_int_gt__gt_;
typedef ::std::vector<double, std::allocator<double> > _std_vector_lt_double_std_allocator_lt_double_gt__gt_;

class CellPopulationMigrationRule2_Overloads : public CellPopulationMigrationRule2{
    public:
    using CellPopulationMigrationRule2::CellPopulationMigrationRule;
    ::std::vector<int, std::allocator<int> > GetIndices(::std::vector<std::shared_ptr<VesselNode<2> >, std::allocator<std::shared_ptr<VesselNode<2> > > > const & rNodes) override {
        PYBIND11_OVERLOAD(
            _std_vector_lt_int_std_allocator_lt_int_gt__gt_,
            CellPopulationMigrationRule2,
            GetIndices,
            rNodes);
    }
    ::std::vector<double, std::allocator<double> > GetNeighbourMovementProbabilities(::std::shared_ptr<VesselNode<2> > pNode, ::std::vector<unsigned int, std::allocator<unsigned int> > neighbourIndices, unsigned int gridIndex) override {
        PYBIND11_OVERLOAD(
            _std_vector_lt_double_std_allocator_lt_double_gt__gt_,
            CellPopulationMigrationRule2,
            GetNeighbourMovementProbabilities,
            pNode, 
neighbourIndices, 
gridIndex);
    }

};
void register_CellPopulationMigrationRule2_class(py::module &m){
py::class_<CellPopulationMigrationRule2 , CellPopulationMigrationRule2_Overloads , std::shared_ptr<CellPopulationMigrationRule2 >  , LatticeBasedMigrationRule<2>  >(m, "CellPopulationMigrationRule2")
        .def(py::init< >())
        .def_static(
            "Create", 
            (::std::shared_ptr<CellPopulationMigrationRule<2> >(*)()) &CellPopulationMigrationRule2::Create, 
            " "  )
        .def(
            "GetIndices", 
            (::std::vector<int, std::allocator<int> >(CellPopulationMigrationRule2::*)(::std::vector<std::shared_ptr<VesselNode<2> >, std::allocator<std::shared_ptr<VesselNode<2> > > > const &)) &CellPopulationMigrationRule2::GetIndices, 
            " " , py::arg("rNodes") )
        .def(
            "SetVolumeFraction", 
            (void(CellPopulationMigrationRule2::*)(::boost::shared_ptr<AbstractCellMutationState>, double)) &CellPopulationMigrationRule2::SetVolumeFraction, 
            " " , py::arg("mutation_state"), py::arg("volume_fraction") )
        .def(
            "GetOccupyingVolumeFraction", 
            (double(CellPopulationMigrationRule2::*)(::boost::shared_ptr<AbstractCellMutationState>)) &CellPopulationMigrationRule2::GetOccupyingVolumeFraction, 
            " " , py::arg("mutation_state") )
    ;
}
