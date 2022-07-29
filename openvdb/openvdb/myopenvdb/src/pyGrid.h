// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0
//
/// @file pyGrid.h
/// @author Peter Cucka
/// @brief Boost.Python wrapper for openvdb::Grid

#ifndef OPENVDB_PYGRID_HAS_BEEN_INCLUDED
#define OPENVDB_PYGRID_HAS_BEEN_INCLUDED

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
//#include <pybind11/embed.h>
//#include <arrayobject.h> // for PyArray_Descr (see pyGrid::arrayTypeId())
#include "openvdb/tools/MeshToVolume.h"
#include "openvdb/tools/VolumeToMesh.h" // for tools::volumeToMesh()
#include "openvdb/openvdb.h"
#include "openvdb/io/Stream.h"
#include "openvdb/math/Math.h" // for math::isExactlyEqual()
#include "openvdb/points/PointDataGrid.h"
#include "openvdb/tools/LevelSetSphere.h"
#include "openvdb/tools/Dense.h"
#include "openvdb/tools/ChangeBackground.h"
#include "openvdb/tools/Prune.h"
#include "openvdb/tools/SignedFloodFill.h"
#include "pyutil.h"
#include "pyAccessor.h" // for pyAccessor::AccessorWrap
#include "pyopenvdb.h"
#include <algorithm> // for std::max()
#include <cstring> // for memcpy()
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace py = pybind11;
//using namespace pybind11::literals;

#ifdef __clang__
// This is a private header, so it's OK to include a "using namespace" directive.
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wheader-hygiene"
#endif

using namespace openvdb::OPENVDB_VERSION_NAME;

#ifdef __clang__
#pragma clang diagnostic pop
#endif


namespace pyopenvdb {

inline py::object
getPyObjectFromGrid(const GridBase::Ptr& grid)
{
    if (!grid) return py::none();

#define CONVERT_BASE_TO_GRID(GridType, grid) \
    if (grid->isType<GridType>()) { \
        return py::cast(grid); \
    }

    CONVERT_BASE_TO_GRID(FloatGrid, grid);
    CONVERT_BASE_TO_GRID(Vec3SGrid, grid);
    CONVERT_BASE_TO_GRID(BoolGrid, grid);
#ifdef PY_OPENVDB_WRAP_ALL_GRID_TYPES
    CONVERT_BASE_TO_GRID(DoubleGrid, grid);
    CONVERT_BASE_TO_GRID(Int32Grid, grid);
    CONVERT_BASE_TO_GRID(Int64Grid, grid);
    CONVERT_BASE_TO_GRID(Vec3IGrid, grid);
    CONVERT_BASE_TO_GRID(Vec3DGrid, grid);
    CONVERT_BASE_TO_GRID(points::PointDataGrid, grid);
#endif
#undef CONVERT_BASE_TO_GRID

    OPENVDB_THROW(TypeError, grid->type() + " is not a supported OpenVDB grid type");
}


inline openvdb::GridBase::Ptr
getGridFromPyObject(const py::object& gridObj)
{
    if (!gridObj) return GridBase::Ptr();

#define CONVERT_GRID_TO_BASE(GridPtrType) \
    { \
        try { \
            return py::cast<GridBase::Ptr>(gridObj); \
        } catch (py::cast_error) { \
        } \
    }

    // Extract a grid pointer of one of the supported types
    // from the input object, then cast it to a base pointer.
    CONVERT_GRID_TO_BASE(FloatGrid::Ptr);
    CONVERT_GRID_TO_BASE(Vec3SGrid::Ptr);
    CONVERT_GRID_TO_BASE(BoolGrid::Ptr);
#ifdef PY_OPENVDB_WRAP_ALL_GRID_TYPES
    CONVERT_GRID_TO_BASE(DoubleGrid::Ptr);
    CONVERT_GRID_TO_BASE(Int32Grid::Ptr);
    CONVERT_GRID_TO_BASE(Int64Grid::Ptr);
    CONVERT_GRID_TO_BASE(Vec3IGrid::Ptr);
    CONVERT_GRID_TO_BASE(Vec3DGrid::Ptr);
    CONVERT_GRID_TO_BASE(points::PointDataGrid::Ptr);
#endif
#undef CONVERT_GRID_TO_BASE

    OPENVDB_THROW(TypeError,
        pyutil::className(gridObj) + " is not a supported OpenVDB grid type");
}


inline openvdb::GridBase::Ptr
getGridFromPyObject(PyObject* gridObj)
{
    return getGridFromPyObject(pyutil::pyBorrow(gridObj));
}

} // namespace pyopenvdb


////////////////////////////////////////


namespace pyGrid {

inline py::object
getGridFromGridBase(GridBase::Ptr grid)
{
    py::object obj;
    try {
        obj = pyopenvdb::getPyObjectFromGrid(grid);
    } catch (openvdb::TypeError& e) {
        PyErr_SetString(PyExc_TypeError, e.what());
        py::error_already_set();
        return pybind11::none();
    }
    return obj;
}


/// GridBase is not exposed in Python because it isn't really needed
/// (and because exposing it would be complicated, requiring wrapping
/// pure virtual functions like GridBase::baseTree()), but there are
/// a few cases where, internally, we need to extract a GridBase::Ptr
/// from a py::object.  Hence this converter.
inline GridBase::Ptr
getGridBaseFromGrid(py::object gridObj)
{
    GridBase::Ptr grid;
    try {
        grid = pyopenvdb::getGridFromPyObject(gridObj);
    } catch (openvdb::TypeError& e) {
        PyErr_SetString(PyExc_TypeError, e.what());
        py::error_already_set();
        return GridBase::Ptr();
    }
    return grid;
}


////////////////////////////////////////


/// Variant of pyutil::castArg() that uses the class name of a given grid type
template<typename GridType, typename T>
inline T
castValueArg(
    py::object obj,
    const char* functionName,
    int argIdx = 0, // args are numbered starting from 1
    const char* expectedType = nullptr)
{
    return pyutil::castArg<T>(obj,
        functionName, pyutil::GridTraits<GridType>::name(), argIdx, expectedType);
}


/// @brief Variant of pyutil::castArg() that uses the class name
/// and @c ValueType of a given grid type
template<typename GridType>
inline typename GridType::ValueType
castValueArg(
    py::object obj,
    const char* functionName,
    int argIdx = 0, // args are numbered starting from 1
    const char* expectedType = nullptr)
{
    return castValueArg<GridType, typename GridType::ValueType>(
        obj, functionName, argIdx, expectedType);
}


////////////////////////////////////////


template<typename GridType>
inline typename GridType::Ptr
copyGrid(GridType& grid)
{
    return grid.copy();
}


template<typename GridType>
inline bool
sharesWith(const GridType& grid, py::object other)
{
    typename GridType::ConstPtr otherGrid;
    try {
        otherGrid = py::cast<typename GridType::Ptr>(other);
    } catch (py::cast_error) {
        return false;
    }
    return (&otherGrid->tree() == &grid.tree());
}


////////////////////////////////////////


template<typename GridType>
inline std::string
getValueType()
{
    return pyutil::GridTraits<GridType>::valueTypeName();
}


template<typename GridType>
inline typename GridType::ValueType
getZeroValue()
{
    return openvdb::zeroVal<typename GridType::ValueType>();
}


template<typename GridType>
inline typename GridType::ValueType
getOneValue()
{
    using ValueT = typename GridType::ValueType;
    return ValueT(openvdb::zeroVal<ValueT>() + 1);
}


template<typename GridType>
inline bool
notEmpty(const GridType& grid)
{
    return !grid.empty();
}


template<typename GridType>
inline typename GridType::ValueType
getGridBackground(const GridType& grid)
{
    return grid.background();
}


template<typename GridType>
inline void
setGridBackground(GridType& grid, py::object obj)
{
    tools::changeBackground(grid.tree(), castValueArg<GridType>(obj, "setBackground"));
}


inline void
setGridName(GridBase::Ptr grid, py::object strObj)
{
    if (grid) {
        if (!strObj) { // if name is None
            grid->removeMeta(GridBase::META_GRID_NAME);
        } else {
            const std::string name = pyutil::castArg<std::string>(
                strObj, "setName", /*className=*/nullptr, /*argIdx=*/1, "str");
            grid->setName(name);
        }
    }
}


inline void
setGridCreator(GridBase::Ptr grid, py::object strObj)
{
    if (grid) {
        if (!strObj) { // if name is None
            grid->removeMeta(GridBase::META_GRID_CREATOR);
        } else {
            const std::string name = pyutil::castArg<std::string>(
                strObj, "setCreator", /*className=*/nullptr, /*argIdx=*/1, "str");
            grid->setCreator(name);
        }
    }
}


inline std::string
getGridClass(GridBase::ConstPtr grid)
{
    return GridBase::gridClassToString(grid->getGridClass());
}


inline void
setGridClass(GridBase::Ptr grid, py::object strObj)
{
    if (!strObj) {
        grid->clearGridClass();
    } else {
        const std::string name = pyutil::castArg<std::string>(
            strObj, "setGridClass", /*className=*/nullptr, /*argIdx=*/1, "str");
        grid->setGridClass(GridBase::stringToGridClass(name));
    }
}


inline std::string
getVecType(GridBase::ConstPtr grid)
{
    return GridBase::vecTypeToString(grid->getVectorType());
}


inline void
setVecType(GridBase::Ptr grid, py::object strObj)
{
    if (!strObj) {
        grid->clearVectorType();
    } else {
        const std::string name = pyutil::castArg<std::string>(
            strObj, "setVectorType", /*className=*/nullptr, /*argIdx=*/1, "str");
        grid->setVectorType(GridBase::stringToVecType(name));
    }
}


inline std::string
gridInfo(GridBase::ConstPtr grid, int verbosity)
{
    std::ostringstream ostr;
    grid->print(ostr, std::max<int>(1, verbosity));
    return ostr.str();
}


////////////////////////////////////////


inline void
setGridTransform(GridBase::Ptr grid, py::object xformObj)
{
    if (grid) {
        if (math::Transform::Ptr xform = pyutil::castArg<math::Transform::Ptr>(
            xformObj, "setTransform", /*className=*/nullptr, /*argIdx=*/1, "Transform"))
        {
            grid->setTransform(xform);
        } else {
            PyErr_SetString(PyExc_ValueError, "null transform");
            py::error_already_set();
        }
    }
}


////////////////////////////////////////


// Helper class to construct a pyAccessor::AccessorWrap for a given grid,
// permitting partial specialization for const vs. non-const grids
template<typename GridType>
struct AccessorHelper
{
    using Wrapper = typename pyAccessor::AccessorWrap<GridType>;
    static Wrapper wrap(typename GridType::Ptr grid)
    {
        if (!grid) {
            PyErr_SetString(PyExc_ValueError, "null grid");
            py::error_already_set();
        }
        return Wrapper(grid);
    }
};

// Specialization for const grids
template<typename GridType>
struct AccessorHelper<const GridType>
{
    using Wrapper = typename pyAccessor::AccessorWrap<const GridType>;
    static Wrapper wrap(typename GridType::ConstPtr grid)
    {
        if (!grid) {
            PyErr_SetString(PyExc_ValueError, "null grid");
            py::error_already_set();
        }
        return Wrapper(grid);
    }
};


/// Return a non-const accessor (wrapped in a pyAccessor::AccessorWrap) for the given grid.
template<typename GridType>
inline typename AccessorHelper<GridType>::Wrapper
getAccessor(typename GridType::Ptr grid)
{
    return AccessorHelper<GridType>::wrap(grid);
}

/// @brief Return a const accessor (wrapped in a pyAccessor::AccessorWrap) for the given grid.
/// @internal Note that the grid pointer is non-const, even though the grid is
/// treated as const.  This is because we don't expose a const grid type in Python.
template<typename GridType>
inline typename AccessorHelper<const GridType>::Wrapper
getConstAccessor(typename GridType::Ptr grid)
{
    return AccessorHelper<const GridType>::wrap(grid);
}


////////////////////////////////////////


template<typename GridType>
inline py::tuple
evalLeafBoundingBox(const GridType& grid)
{
    CoordBBox bbox;
    grid.tree().evalLeafBoundingBox(bbox);
    return py::make_tuple(bbox.min(), bbox.max());
}


template<typename GridType>
inline Coord
evalLeafDim(const GridType& grid)
{
    Coord dim;
    grid.tree().evalLeafDim(dim);
    return dim;
}


template<typename GridType>
inline py::tuple
evalActiveVoxelBoundingBox(const GridType& grid)
{
    CoordBBox bbox = grid.evalActiveVoxelBoundingBox();
    return py::make_tuple(bbox.min(), bbox.max());
}


template<typename GridType>
inline py::tuple
getNodeLog2Dims(const GridType& grid)
{
    std::vector<Index> dims;
    grid.tree().getNodeLog2Dims(dims);
    py::list lst;
    for (size_t i = 0, N = dims.size(); i < N; ++i) {
        lst.append(dims[i]);
    }
    return py::tuple(lst);
}


template<typename GridType>
inline Index
treeDepth(const GridType& grid)
{
    return grid.tree().treeDepth();
}


template<typename GridType>
inline Index32
leafCount(const GridType& grid)
{
    return grid.tree().leafCount();
}


template<typename GridType>
inline Index32
nonLeafCount(const GridType& grid)
{
    return grid.tree().nonLeafCount();
}


template<typename GridType>
inline Index64
activeLeafVoxelCount(const GridType& grid)
{
    return grid.tree().activeLeafVoxelCount();
}


template<typename GridType>
inline py::tuple
evalMinMax(const GridType& grid)
{
    typename GridType::ValueType vmin, vmax;
    grid.tree().evalMinMax(vmin, vmax);
    return py::make_tuple(vmin, vmax);
}


template<typename GridType>
inline py::tuple
getIndexRange(const GridType& grid)
{
    CoordBBox bbox;
    grid.tree().getIndexRange(bbox);
    return py::make_tuple(bbox.min(), bbox.max());
}


////////////////////////////////////////


inline py::dict
getAllMetadata(GridBase::ConstPtr grid)
{
    if (grid) return py::cast(static_cast<const MetaMap&>(*grid));
    return py::dict();
}


inline void
replaceAllMetadata(GridBase::Ptr grid, const MetaMap& metadata)
{
    if (grid) {
        grid->clearMetadata();
        for (MetaMap::ConstMetaIterator it = metadata.beginMeta();
            it != metadata.endMeta(); ++it)
        {
            if (it->second) grid->insertMeta(it->first, *it->second);
        }
    }
}


inline void
updateMetadata(GridBase::Ptr grid, const MetaMap& metadata)
{
    if (grid) {
        for (MetaMap::ConstMetaIterator it = metadata.beginMeta();
            it != metadata.endMeta(); ++it)
        {
            if (it->second) {
                grid->removeMeta(it->first);
                grid->insertMeta(it->first, *it->second);
            }
        }
    }
}


inline py::dict
getStatsMetadata(GridBase::ConstPtr grid)
{
    MetaMap::ConstPtr metadata;
    if (grid) metadata = grid->getStatsMetadata();
    return py::cast(*metadata);
}


inline py::object
getMetadataKeys(GridBase::ConstPtr grid)
{
    if (grid) {
        // Return an iterator over the "keys" view of a dict.
        return py::make_iterator(grid->beginMeta(), grid->endMeta());
    }
    return py::none();
}


inline py::object
getMetadata(GridBase::ConstPtr grid, const std::string& name)
{
    if (!grid) return py::none();

    Metadata::ConstPtr metadata = (*grid)[name];
    if (!metadata) {
        throw py::key_error(name);
    }
    return py::cast(metadata);;
}


inline void
setMetadata(GridBase::Ptr grid, const std::string& name, const Metadata& metadata)
{
    if (!grid) return;

    grid->removeMeta(name);
    grid->insertMeta(name, metadata);
}


inline void
removeMetadata(GridBase::Ptr grid, const std::string& name)
{
    if (grid) {
        Metadata::Ptr metadata = (*grid)[name];
        if (!metadata) {
            throw py::key_error(name);
        }
        grid->removeMeta(name);
    }
}


inline bool
hasMetadata(GridBase::ConstPtr grid, const std::string& name)
{
    if (grid) return ((*grid)[name].get() != nullptr);
    return false;
}


////////////////////////////////////////


template<typename GridType>
inline void
prune(GridType& grid, py::object tolerance)
{
    tools::prune(grid.tree(), castValueArg<GridType>(tolerance, "prune"));
}


template<typename GridType>
inline void
pruneInactive(GridType& grid, py::object valObj)
{
    if (valObj.is_none()) {
        tools::pruneInactive(grid.tree());
    } else {
        tools::pruneInactiveWithValue(
            grid.tree(), castValueArg<GridType>(valObj, "pruneInactive"));
    }
}


template<typename GridType>
inline void
fill(GridType& grid, py::object minObj, py::object maxObj,
    py::object valObj, bool active)
{
    const Coord
        bmin = castValueArg<GridType, Coord>(minObj, "fill", 1, "tuple(int, int, int)"),
        bmax = castValueArg<GridType, Coord>(maxObj, "fill", 2, "tuple(int, int, int)");
    grid.fill(CoordBBox(bmin, bmax), castValueArg<GridType>(valObj, "fill", 3), active);
}


template<typename GridType>
inline void
signedFloodFill(GridType& grid)
{
    tools::signedFloodFill(grid.tree());
}


////////////////////////////////////////


template<typename GridType>
inline void
copyToArray(GridType& grid, py::buffer buffer, const Coord& coord)
{
    using ValueT = typename GridType::ValueType;
    using DenseT = typename tools::Dense<ValueT>;

    py::buffer_info info = buffer.request();
    if (info.ndim != 3) {
        std::ostringstream os;
        os << "expected 3-dimensional array, found " << info.ndim << "-dimensional array";
        throw py::value_error(os.str());
    }
    Coord maxCoord(coord[0] + info.shape[0], coord[1] + info.shape[1], coord[2] + info.shape[2]);
    CoordBBox bbox(coord, maxCoord);
    DenseT valArray(bbox, static_cast<ValueT*>(info.ptr));
    tools::copyToDense(grid, valArray);
}


template<typename GridType>
py::array_t<typename GridType::ValueType>
toArray(GridType& grid, const Coord& coord)
{
    using ValueT = typename GridType::ValueType;
    using DenseT = typename tools::Dense<ValueT>;

    const std::size_t itemsize = sizeof(ValueT);
    const std::string format = py::format_descriptor<ValueT>::value;
    const std::size_t ndim = 3;
    const std::vector<std::size_t> shape = {
        8,  // 8
        8,  // 8
        8   // 8
    };
    const std::vector<std::size_t> strides = {
        8 * 8 * itemsize,
        8 * itemsize,
        itemsize,
    };

    Coord maxCoord(coord[0] + shape[0] - 1, coord[1] + shape[1] - 1, coord[2] + shape[2] - 1);
    CoordBBox bbox(coord, maxCoord);
    DenseT valArray(bbox);
    tools::copyToDense(grid, valArray);
    return py::array_t<ValueT>(py::buffer_info(
        valArray.data(),  // Pointer to the underlying storage
        itemsize,         // Size of individual items in bytes
        format,           // set to format_descriptor<T>::format()
        ndim,             // Number of dimensions
        shape,            // Shape of the tensor (1 entry per dimension)
        strides           // Number of bytes between adjacent entries
    ));
}


template<typename GridType>
inline void
copyFromArray(GridType& grid, py::buffer buffer, const Coord& coord, const typename GridType::ValueType& tolerance)
{
    using ValueT = typename GridType::ValueType;
    using DenseT = typename tools::Dense<ValueT>;

    py::buffer_info info = buffer.request();
    if (info.ndim != 3) {
        std::ostringstream os;
        os << "expected 3-dimensional array, found " << info.ndim << "-dimensional array";
        throw py::value_error(os.str());
    }
    Coord maxCoord(coord[0] + info.shape[0], coord[1] + info.shape[1], coord[2] + info.shape[2]);
    CoordBBox bbox(coord, maxCoord);
    DenseT valArray(bbox, static_cast<ValueT*>(info.ptr));
    // tools::copyFromDense(grid, valArray, tolerance);
}


////////////////////////////////////////


template<typename GridType, typename IterType>
inline void
applyMap(const char* methodName, GridType& grid, py::object funcObj)
{
    using ValueT = typename GridType::ValueType;

    for (IterType it = grid.tree().template begin<IterType>(); it; ++it) {
        // Evaluate the functor.
        py::object result = funcObj(*it);

        // Verify that the result is of type GridType::ValueType.
        ValueT val;
        try {
            val = py::cast<ValueT>(result);
        } catch (py::cast_error) {
            PyErr_Format(PyExc_TypeError,
                "expected callable argument to %s.%s() to return %s, found %s",
                pyutil::GridTraits<GridType>::name(),
                methodName,
                openvdb::typeNameAsString<ValueT>(),
                pyutil::className(result).c_str());
            py::error_already_set();
        }

        it.setValue(val);
    }
}


template<typename GridType>
inline void
mapOn(GridType& grid, py::object funcObj)
{
    applyMap<GridType, typename GridType::ValueOnIter>("mapOn", grid, funcObj);
}


template<typename GridType>
inline void
mapOff(GridType& grid, py::object funcObj)
{
    applyMap<GridType, typename GridType::ValueOffIter>("mapOff", grid, funcObj);
}


template<typename GridType>
inline void
mapAll(GridType& grid, py::object funcObj)
{
    applyMap<GridType, typename GridType::ValueAllIter>("mapAll", grid, funcObj);
}


////////////////////////////////////////


template<typename GridType>
struct TreeCombineOp
{
    using TreeT = typename GridType::TreeType;
    using ValueT = typename GridType::ValueType;

    TreeCombineOp(py::object _op): op(_op) {}
    void operator()(const ValueT& a, const ValueT& b, ValueT& result)
    {
        py::object resultObj = op(a, b);

        ValueT val;
        try {
            val = py::cast<ValueT>(resultObj);
        } catch (py::cast_error) {
            PyErr_Format(PyExc_TypeError,
                "expected callable argument to %s.combine() to return %s, found %s",
                pyutil::GridTraits<GridType>::name(),
                openvdb::typeNameAsString<ValueT>(),
                pyutil::className(resultObj).c_str());
            py::error_already_set();
        }

        result = val;
    }
    py::object op;
};


template<typename GridType>
inline void
combine(GridType& grid, py::object otherGridObj, py::object funcObj)
{
    using GridPtr = typename GridType::Ptr;
    GridPtr otherGrid = castValueArg<GridType, GridPtr>(otherGridObj,
        "combine", 1, pyutil::GridTraits<GridType>::name());
    TreeCombineOp<GridType> op(funcObj);
    grid.tree().combine(otherGrid->tree(), op, /*prune=*/true);
}


////////////////////////////////////////


template<typename GridType>
inline typename GridType::Ptr
createLevelSetSphere(float radius, const openvdb::Vec3f& center, float voxelSize, float halfWidth)
{
    return tools::createLevelSetSphere<GridType>(radius, center, voxelSize, halfWidth);
}


////////////////////////////////////////


template<typename GridT, typename IterT> class IterWrap; // forward declaration

//
// Type traits for various iterators
//
template<typename GridT, typename IterT> struct IterTraits
{
    // IterT    the type of the iterator
    // name()   function returning the base name of the iterator type (e.g., "ValueOffIter")
    // descr()  function returning a string describing the iterator
    // begin()  function returning a begin iterator for a given grid
};

template<typename GridT> struct IterTraits<GridT, typename GridT::ValueOnCIter>
{
    using IterT = typename GridT::ValueOnCIter;
    static std::string name() { return "ValueOnCIter"; }
    static std::string descr()
    {
        return std::string("Read-only iterator over the active values (tile and voxel)\nof a ")
            + pyutil::GridTraits<typename std::remove_const<GridT>::type>::name();
    }
    static IterWrap<const GridT, IterT> begin(typename GridT::Ptr g)
    {
        return IterWrap<const GridT, IterT>(g, g->cbeginValueOn());
    }
}; // IterTraits<ValueOnCIter>

template<typename GridT> struct IterTraits<GridT, typename GridT::ValueOffCIter>
{
    using IterT = typename GridT::ValueOffCIter;
    static std::string name() { return "ValueOffCIter"; }
    static std::string descr()
    {
        return std::string("Read-only iterator over the inactive values (tile and voxel)\nof a ")
            + pyutil::GridTraits<typename std::remove_const<GridT>::type>::name();
    }
    static IterWrap<const GridT, IterT> begin(typename GridT::Ptr g)
    {
        return IterWrap<const GridT, IterT>(g, g->cbeginValueOff());
    }
}; // IterTraits<ValueOffCIter>

template<typename GridT> struct IterTraits<GridT, typename GridT::ValueAllCIter>
{
    using IterT = typename GridT::ValueAllCIter;
    static std::string name() { return "ValueAllCIter"; }
    static std::string descr()
    {
        return std::string("Read-only iterator over all tile and voxel values of a ")
            + pyutil::GridTraits<typename std::remove_const<GridT>::type>::name();
    }
    static IterWrap<const GridT, IterT> begin(typename GridT::Ptr g)
    {
        return IterWrap<const GridT, IterT>(g, g->cbeginValueAll());
    }
}; // IterTraits<ValueAllCIter>

template<typename GridT> struct IterTraits<GridT, typename GridT::ValueOnIter>
{
    using IterT = typename GridT::ValueOnIter;
    static std::string name() { return "ValueOnIter"; }
    static std::string descr()
    {
        return std::string("Read/write iterator over the active values (tile and voxel)\nof a ")
            + pyutil::GridTraits<typename std::remove_const<GridT>::type>::name();
    }
    static IterWrap<GridT, IterT> begin(typename GridT::Ptr g)
    {
        return IterWrap<GridT, IterT>(g, g->beginValueOn());
    }
}; // IterTraits<ValueOnIter>

template<typename GridT> struct IterTraits<GridT, typename GridT::ValueOffIter>
{
    using IterT = typename GridT::ValueOffIter;
    static std::string name() { return "ValueOffIter"; }
    static std::string descr()
    {
        return std::string("Read/write iterator over the inactive values (tile and voxel)\nof a ")
            + pyutil::GridTraits<typename std::remove_const<GridT>::type>::name();
    }
    static IterWrap<GridT, IterT> begin(typename GridT::Ptr g)
    {
        return IterWrap<GridT, IterT>(g, g->beginValueOff());
    }
}; // IterTraits<ValueOffIter>

template<typename GridT> struct IterTraits<GridT, typename GridT::ValueAllIter>
{
    using IterT = typename GridT::ValueAllIter;
    static std::string name() { return "ValueAllIter"; }
    static std::string descr()
    {
        return std::string("Read/write iterator over all tile and voxel values of a ")
            + pyutil::GridTraits<typename std::remove_const<GridT>::type>::name();
    }
    static IterWrap<GridT, IterT> begin(typename GridT::Ptr g)
    {
        return IterWrap<GridT, IterT>(g, g->beginValueAll());
    }
}; // IterTraits<ValueAllIter>


////////////////////////////////////////


// Helper class to modify a grid through a non-const iterator
template<typename GridT, typename IterT>
struct IterItemSetter
{
    using ValueT = typename GridT::ValueType;
    static void setValue(const IterT& iter, const ValueT& val) { iter.setValue(val); }
    static void setActive(const IterT& iter, bool on) { iter.setActiveState(on); }
};

// Partial specialization for const iterators
template<typename GridT, typename IterT>
struct IterItemSetter<const GridT, IterT>
{
    using ValueT = typename GridT::ValueType;
    static void setValue(const IterT&, const ValueT&)
    {
        PyErr_SetString(PyExc_AttributeError, "can't set attribute 'value'");
        py::error_already_set();
    }
    static void setActive(const IterT&, bool /*on*/)
    {
        PyErr_SetString(PyExc_AttributeError, "can't set attribute 'active'");
        py::error_already_set();
    }
};


/// @brief Value returned by the next() method of a grid's value iterator
/// @details This class allows both dictionary-style (e.g., items["depth"]) and
/// attribute access (e.g., items.depth) to the items returned by an iterator.
/// @todo Create a reusable base class for "named dicts" like this?
template<typename _GridT, typename _IterT>
class IterValueProxy
{
public:
    using GridT = _GridT;
    using IterT = _IterT;
    using ValueT = typename GridT::ValueType;
    using SetterT = IterItemSetter<GridT, IterT>;

    IterValueProxy(typename GridT::ConstPtr grid, const IterT& iter): mGrid(grid), mIter(iter) {}

    IterValueProxy copy() const { return *this; }

    typename GridT::ConstPtr parent() const { return mGrid; }

    ValueT getValue() const { return *mIter; }
    bool getActive() const { return mIter.isValueOn(); }
    Index getDepth() const { return mIter.getDepth(); }
    Coord getBBoxMin() const { return mIter.getBoundingBox().min(); }
    Coord getBBoxMax() const { return mIter.getBoundingBox().max(); }
    Index64 getVoxelCount() const { return mIter.getVoxelCount(); }

    void setValue(const ValueT& val) { SetterT::setValue(mIter, val); }
    void setActive(bool on) { SetterT::setActive(mIter, on); }

    /// Return this dictionary's keys as a list of C strings.
    static const char* const * keys()
    {
        static const char* const sKeys[] = {
            "value", "active", "depth", "min", "max", "count", nullptr
        };
        return sKeys;
    }

    /// Return @c true if the given string is a valid key.
    static bool hasKey(const std::string& key)
    {
        for (int i = 0; keys()[i] != nullptr; ++i) {
            if (key == keys()[i]) return true;
        }
        return false;
    }

    /// Return this dictionary's keys as a Python list of Python strings.
    static py::list getKeys()
    {
        py::list keyList;
        for (int i = 0; keys()[i] != nullptr; ++i) keyList.append(keys()[i]);
        return keyList;
    }

    /// @brief Return the value for the given key.
    /// @throw KeyError if the key is invalid
    py::object getItem(py::object keyObj) const
    {
        std::string key;
        try {
            key = py::cast<std::string>(keyObj);
        } catch (py::cast_error) {
            PyErr_SetObject(PyExc_KeyError, PyUnicode_FromFormat("%s", keyObj.attr("__repr__")()));
            py::error_already_set();
            return py::none();
        }
        if (key == "value") return py::cast(this->getValue());
        else if (key == "active") return py::cast(this->getActive());
        else if (key == "depth") return py::cast(this->getDepth());
        else if (key == "min") return py::cast(this->getBBoxMin());
        else if (key == "max") return py::cast(this->getBBoxMax());
        else if (key == "count") return py::cast(this->getVoxelCount());
        return py::none();
    }

    /// @brief Set the value for the given key.
    /// @throw KeyError if the key is invalid
    /// @throw AttributeError if the key refers to a read-only item
    void setItem(py::object keyObj, py::object valObj)
    {
        std::string key;
        try {
            key = py::cast<std::string>(keyObj);
        } catch (py::cast_error) {
            PyErr_SetObject(PyExc_KeyError,
                PyUnicode_FromFormat("'%s'", keyObj.attr("__repr__")()));
            py::error_already_set();
        }
        if (key == "value") {
            this->setValue(py::cast<ValueT>(valObj)); return;
        } else if (key == "active") {
            this->setActive(py::cast<bool>(valObj)); return;
        } else if (this->hasKey(key)) {
            PyErr_SetObject(PyExc_AttributeError,
                PyUnicode_FromFormat("can't set attribute '%s'", keyObj.attr("__repr__")()));
            py::error_already_set();
        }
    }

    bool operator==(const IterValueProxy& other) const
    {
        return (other.getActive() == this->getActive()
            && other.getDepth() == this->getDepth()
            && math::isExactlyEqual(other.getValue(), this->getValue())
            && other.getBBoxMin() == this->getBBoxMin()
            && other.getBBoxMax() == this->getBBoxMax()
            && other.getVoxelCount() == this->getVoxelCount());
    }
    bool operator!=(const IterValueProxy& other) const { return !(*this == other); }

    /// Print this dictionary to a stream.
    std::ostream& put(std::ostream& os) const
    {
        // valuesAsStrings = ["%s: %s" % key, repr(this[key]) for key in this.keys()]
        py::list valuesAsStrings;
        for (int i = 0; this->keys()[i] != nullptr; ++i) {
            py::str
                key(this->keys()[i]),
                val(this->getItem(key).attr("__repr__")());
            valuesAsStrings.append(PyUnicode_FromFormat("'%s': %s", key, val));
        }
        // print ", ".join(valuesAsStrings)
        py::object joined = py::str(", ").attr("join")(valuesAsStrings);
        std::string s = py::cast<std::string>(joined);
        os << "{" << s << "}";
        return os;
    }
    /// Return a string describing this dictionary.
    std::string info() const { std::ostringstream os; os << *this; return os.str(); }

private:
    // To keep the iterator's grid from being deleted (leaving the iterator dangling),
    // store a shared pointer to the grid.
    const typename GridT::ConstPtr mGrid;
    const IterT mIter; // the iterator may not be incremented
}; // class IterValueProxy


template<typename GridT, typename IterT>
inline std::ostream&
operator<<(std::ostream& os, const IterValueProxy<GridT, IterT>& iv) { return iv.put(os); }


////////////////////////////////////////


/// Wrapper for a grid's value iterator classes
template<typename _GridT, typename _IterT>
class IterWrap
{
public:
    using GridT = _GridT;
    using IterT = _IterT;
    using ValueT = typename GridT::ValueType;
    using IterValueProxyT = IterValueProxy<GridT, IterT>;
    using Traits = IterTraits<GridT, IterT>;

    IterWrap(typename GridT::ConstPtr grid, const IterT& iter): mGrid(grid), mIter(iter) {}

    typename GridT::ConstPtr parent() const { return mGrid; }

    /// Return an IterValueProxy for the current iterator position.
    IterValueProxyT next()
    {
        if (!mIter) {
            PyErr_SetString(PyExc_StopIteration, "no more values");
            py::error_already_set();
        }
        IterValueProxyT result(mGrid, mIter);
        ++mIter;
        return result;
    }

    static py::object returnSelf(const py::object& obj) { return obj; }

    /// @brief Define a Python wrapper class for this C++ class and another for
    /// the IterValueProxy class returned by iterators of this type.
    static void wrap(py::module_ &m)
    {
        const std::string
            gridClassName = pyutil::GridTraits<typename std::remove_const<GridT>::type>::name(),
            iterClassName = gridClassName + Traits::name(),
            valueClassName = gridClassName + Traits::name() + "Value";

        py::class_<IterWrap>(m,
            iterClassName.c_str(),
            /*docstring=*/Traits::descr().c_str())

            .def_property_readonly("parent", &IterWrap::parent,
                ("the " + gridClassName + " over which to iterate").c_str())

            .def("next", &IterWrap::next, ("next() -> " + valueClassName).c_str())
            .def("__next__", &IterWrap::next, ("__next__() -> " + valueClassName).c_str())
            .def("__iter__", &returnSelf);

        py::class_<IterValueProxyT>(m,
            valueClassName.c_str(),
            /*docstring=*/("Proxy for a tile or voxel value in a " + gridClassName).c_str())

            .def("copy", &IterValueProxyT::copy,
                ("copy() -> " + valueClassName + "\n\n"
                "Return a shallow copy of this value, i.e., one that shares\n"
                "its data with the original.").c_str())

            .def_property_readonly("parent", &IterValueProxyT::parent,
                ("the " + gridClassName + " to which this value belongs").c_str())

            .def("__str__", &IterValueProxyT::info)
            .def("__repr__", &IterValueProxyT::info)

            .def("__eq__", &IterValueProxyT::operator==)
            .def("__ne__", &IterValueProxyT::operator!=)

            .def_property("value", &IterValueProxyT::getValue, &IterValueProxyT::setValue,
                "value of this tile or voxel")
            .def_property("active", &IterValueProxyT::getActive, &IterValueProxyT::setActive,
                "active state of this tile or voxel")
            .def_property_readonly("depth", &IterValueProxyT::getDepth,
                "tree depth at which this value is stored")
            .def_property_readonly("min", &IterValueProxyT::getBBoxMin,
                "lower bound of the axis-aligned bounding box of this tile or voxel")
            .def_property_readonly("max", &IterValueProxyT::getBBoxMax,
                "upper bound of the axis-aligned bounding box of this tile or voxel")
            .def_property_readonly("count", &IterValueProxyT::getVoxelCount,
                "number of voxels spanned by this value")

            .def_static("keys", &IterValueProxyT::getKeys,
                "keys() -> list\n\n"
                "Return a list of keys for this tile or voxel.")
            .def_static("__contains__", &IterValueProxyT::hasKey,
                "__contains__(key) -> bool\n\n"
                "Return True if the given key exists.")
            .def("__getitem__", &IterValueProxyT::getItem,
                "__getitem__(key) -> value\n\n"
                "Return the value of the item with the given key.")
            .def("__setitem__", &IterValueProxyT::getItem,
                "__setitem__(key, value)\n\n"
                "Set the value of the item with the given key.");
    }

private:
    // To keep this iterator's grid from being deleted, leaving the iterator dangling,
    // store a shared pointer to the grid.
    const typename GridT::ConstPtr mGrid;
    IterT mIter;
}; // class IterWrap


////////////////////////////////////////


/// Return a tuple representing the state of the given Grid.
template<typename GridT>
py::tuple getstate(const py::object& self)
{
    using GridPtrT = typename GridT::Ptr;

    // Extract a Grid from the Python object.
    GridPtrT grid;
    try {
        grid = py::cast<GridPtrT>(self);
    } catch (py::cast_error) {
        throw std::runtime_error("Invalid cast!");
    }

    // Serialize the Grid to a string.
    std::ostringstream ostr(std::ios_base::binary);
    {
        openvdb::io::Stream strm(ostr);
        strm.setGridStatsMetadataEnabled(false);
        strm.write(openvdb::GridPtrVec(1, grid));
    }
    // Construct a state tuple comprising the Python object's __dict__
    // and the serialized Grid.=
    py::bytes bytesObj(ostr.str());
    return py::make_tuple(self.attr("__dict__"), bytesObj);
}

/// Restore the given Grid to a saved state.
template<typename GridT>
std::pair<GridT, py::dict> setstate(const py::tuple &state)
{
    using GridPtrT = typename GridT::Ptr;

    GridT grid;

    bool badState = (py::len(state) != 2);

    std::string serialized;
    if (!badState) {
        // Extract the sequence containing the serialized Grid.
        py::object bytesObj = state[1];
        badState = true;
        if (PyBytes_Check(bytesObj.ptr())) {
            // Convert the "bytes" sequence to a byte string.
            char* buf = nullptr;
            Py_ssize_t length = 0;
            if (-1 != PyBytes_AsStringAndSize(bytesObj.ptr(), &buf, &length)) {
                if (buf != nullptr && length > 0) {
                    serialized.assign(buf, buf + length);
                    badState = false;
                }
            }
        }
    }

    if (badState) {
        throw py::value_error("expected (dict, bytes) tuple in call to __setstate__;");
    }

    // Restore the internal state of the C++ object.
    GridPtrVecPtr grids;
    {
        std::istringstream istr(serialized, std::ios_base::binary);
        io::Stream strm(istr);
        grids = strm.getGrids(); // (note: file-level metadata is ignored)
    }
    if (grids && !grids->empty()) {
        if (GridPtrT savedGrid = gridPtrCast<GridT>((*grids)[0])) {
            grid.MetaMap::operator=(*savedGrid); ///< @todo add a Grid::setMetadata() method?
            grid.setTransform(savedGrid->transformPtr());
            grid.setTree(savedGrid->treePtr());
        }
    }
    return std::make_pair(grid, py::cast<py::dict>(state[0]));
}


////////////////////////////////////////


/// Create a Python wrapper for a particular template instantiation of Grid.
template<typename GridType>
inline void
exportGrid(py::module_ &m)
{
    using ValueT = typename GridType::ValueType;
    using GridPtr = typename GridType::Ptr;
    using Traits = pyutil::GridTraits<GridType>;
    using DenseT = typename tools::Dense<ValueT>;

    using ValueOnCIterT = typename GridType::ValueOnCIter;
    using ValueOffCIterT = typename GridType::ValueOffCIter;
    using ValueAllCIterT = typename GridType::ValueAllCIter;
    using ValueOnIterT = typename GridType::ValueOnIter;
    using ValueOffIterT = typename GridType::ValueOffIter;
    using ValueAllIterT = typename GridType::ValueAllIter;

    math::Transform::Ptr (GridType::*getTransform)() = &GridType::transformPtr;

    const std::string pyGridTypeName = Traits::name();

    // Define the Grid wrapper class and make it the current scope.
    {
        py::class_<GridType, /*HeldType=*/GridPtr> clss(m,
            /*classname=*/pyGridTypeName.c_str(),
            /*docstring=*/(Traits::descr()).c_str()
        );

        clss.def(py::init<const ValueT&>(), py::arg("background")=pyGrid::getZeroValue<GridType>(),
                "Initialize with the given background value.")

            .def("copy", &pyGrid::copyGrid<GridType>,
                ("copy() -> " + pyGridTypeName + "\n\n"
                "Return a shallow copy of this grid, i.e., a grid\n"
                "that shares its voxel data with this grid.").c_str())
            .def("deepCopy", &GridType::deepCopy,
                ("deepCopy() -> " + pyGridTypeName + "\n\n"
                "Return a deep copy of this grid.\n").c_str())

            .def(py::pickle(&pyGrid::getstate<GridType>, &pyGrid::setstate<GridType>))

            .def("sharesWith", &pyGrid::sharesWith<GridType>,
                ("sharesWith(" + pyGridTypeName + ") -> bool\n\n"
                "Return True if this grid shares its voxel data with the given grid.").c_str())

            /// @todo Any way to set a docstring for a class property?
            .def_property_readonly("valueTypeName", &pyGrid::getValueType<GridType>)
                /// @todo docstring = "name of this grid's value type"
            .def_property_readonly("zeroValue", &pyGrid::getZeroValue<GridType>)
                /// @todo docstring = "zero, as expressed in this grid's value type"
            .def_property_readonly("oneValue", &pyGrid::getOneValue<GridType>)
                /// @todo docstring = "one, as expressed in this grid's value type"
            /// @todo Is Grid.typeName ever needed?
            //.add_static_property("typeName", &GridType::gridType)
                /// @todo docstring = to "name of this grid's type"

            .def_property("background",
                &pyGrid::getGridBackground<GridType>, &pyGrid::setGridBackground<GridType>,
                "value of this grid's background voxels")
            .def_property("name", &GridType::getName, &pyGrid::setGridName,
                "this grid's name")
            .def_property("creator", &GridType::getCreator, &pyGrid::setGridCreator,
                "description of this grid's creator")

            .def_property("transform", getTransform, &pyGrid::setGridTransform,
                "transform associated with this grid")

            .def_property("gridClass", &pyGrid::getGridClass, &pyGrid::setGridClass,
                "the class of volumetric data (level set, fog volume, etc.)\nstored in this grid")

            .def_property("vectorType", &pyGrid::getVecType, &pyGrid::setVecType,
                "how transforms are applied to values stored in this grid")

            .def("getAccessor", &pyGrid::getAccessor<GridType>,
                ("getAccessor() -> " + pyGridTypeName + "Accessor\n\n"
                "Return an accessor that provides random read and write access\n"
                "to this grid's voxels.").c_str())
            .def("getConstAccessor", &pyGrid::getConstAccessor<GridType>,
                ("getConstAccessor() -> " + pyGridTypeName + "Accessor\n\n"
                "Return an accessor that provides random read-only access\n"
                "to this grid's voxels.").c_str())

            .def_property("metadata", &pyGrid::getAllMetadata, &pyGrid::replaceAllMetadata,
                "dict of this grid's metadata\n\n"
                "Setting this attribute replaces all of this grid's metadata,\n"
                "but mutating it in place has no effect on the grid, since\n"
                "the value of this attribute is a only a copy of the metadata.\n"
                "Use either indexing or updateMetadata() to mutate metadata in place.")
            .def("updateMetadata", &pyGrid::updateMetadata,
                "updateMetadata(dict)\n\n"
                "Add metadata to this grid, replacing any existing items\n"
                "having the same names as the new items.")

            .def("addStatsMetadata", &GridType::addStatsMetadata,
                "addStatsMetadata()\n\n"
                "Add metadata to this grid comprising the current values\n"
                "of statistics like the active voxel count and bounding box.\n"
                "(This metadata is not automatically kept up-to-date with\n"
                "changes to this grid.)")
            .def("getStatsMetadata", &pyGrid::getStatsMetadata,
                "getStatsMetadata() -> dict\n\n"
                "Return a (possibly empty) dict containing just the metadata\n"
                "that was added to this grid with addStatsMetadata().")

            .def("__getitem__", &pyGrid::getMetadata,
                "__getitem__(name) -> value\n\n"
                "Return the metadata value associated with the given name.")
            .def("__setitem__", &pyGrid::setMetadata,
                "__setitem__(name, value)\n\n"
                "Add metadata to this grid, replacing any existing item having\n"
                "the same name as the new item.")
            .def("__delitem__", &pyGrid::removeMetadata,
                "__delitem__(name)\n\n"
                "Remove the metadata with the given name.")
            .def("__contains__", &pyGrid::hasMetadata,
                "__contains__(name) -> bool\n\n"
                "Return True if this grid contains metadata with the given name.")
            .def("__iter__", &pyGrid::getMetadataKeys,
                "__iter__() -> iterator\n\n"
                "Return an iterator over this grid's metadata keys.")
            .def("iterkeys", &pyGrid::getMetadataKeys,
                "iterkeys() -> iterator\n\n"
                "Return an iterator over this grid's metadata keys.")

            .def_property("saveFloatAsHalf",
                &GridType::saveFloatAsHalf, &GridType::setSaveFloatAsHalf,
                "if True, write floating-point voxel values as 16-bit half floats")

            //
            // Statistics
            //
            .def("memUsage", &GridType::memUsage,
                "memUsage() -> int\n\n"
                "Return the memory usage of this grid in bytes.")

            .def("evalLeafBoundingBox", &pyGrid::evalLeafBoundingBox<GridType>,
                "evalLeafBoundingBox() -> xyzMin, xyzMax\n\n"
                "Return the coordinates of opposite corners of the axis-aligned\n"
                "bounding box of all leaf nodes.")
            .def("evalLeafDim", &pyGrid::evalLeafDim<GridType>,
                "evalLeafDim() -> x, y, z\n\n"
                "Return the dimensions of the axis-aligned bounding box\n"
                "of all leaf nodes.")

            .def("evalActiveVoxelBoundingBox", &pyGrid::evalActiveVoxelBoundingBox<GridType>,
                "evalActiveVoxelBoundingBox() -> xyzMin, xyzMax\n\n"
                "Return the coordinates of opposite corners of the axis-aligned\n"
                "bounding box of all active voxels.")
            .def("evalActiveVoxelDim", &GridType::evalActiveVoxelDim,
                "evalActiveVoxelDim() -> x, y, z\n\n"
                "Return the dimensions of the axis-aligned bounding box of all\n"
                "active voxels.")

            .def_property_readonly("treeDepth", &pyGrid::treeDepth<GridType>,
                "depth of this grid's tree from root node to leaf node")
            .def("nodeLog2Dims", &pyGrid::getNodeLog2Dims<GridType>,
                "list of Log2Dims of the nodes of this grid's tree\n"
                "in order from root to leaf")

            .def("leafCount", &pyGrid::leafCount<GridType>,
                "leafCount() -> int\n\n"
                "Return the number of leaf nodes in this grid's tree.")
            .def("nonLeafCount", &pyGrid::nonLeafCount<GridType>,
                "nonLeafCount() -> int\n\n"
                "Return the number of non-leaf nodes in this grid's tree.")

            .def("activeVoxelCount", &GridType::activeVoxelCount,
                "activeVoxelCount() -> int\n\n"
                "Return the number of active voxels in this grid.")
            .def("activeLeafVoxelCount", &pyGrid::activeLeafVoxelCount<GridType>,
                "activeLeafVoxelCount() -> int\n\n"
                "Return the number of active voxels that are stored\n"
                "in the leaf nodes of this grid's tree.")

            .def("evalMinMax", &pyGrid::evalMinMax<GridType>,
                "evalMinMax() -> min, max\n\n"
                "Return the minimum and maximum active values in this grid.")

            .def("getIndexRange", &pyGrid::getIndexRange<GridType>,
                "getIndexRange() -> min, max\n\n"
                "Return the minimum and maximum coordinates that are represented\n"
                "in this grid.  These might include background voxels.")

            .def("info", &pyGrid::gridInfo,
                py::arg("verbosity")=1,
                "info(verbosity=1) -> str\n\n"
                "Return a string containing information about this grid\n"
                "with a specified level of verbosity.\n")

            //
            // Tools
            //
            .def("fill", &pyGrid::fill<GridType>,
                py::arg("min"), py::arg("max"), py::arg("value"), py::arg("active")=true,
                "fill(min, max, value, active=True)\n\n"
                "Set all voxels within a given axis-aligned box to\n"
                "a constant value (either active or inactive).")
            .def("signedFloodFill", &pyGrid::signedFloodFill<GridType>,
                "signedFloodFill()\n\n"
                "Propagate the sign from a narrow-band level set into inactive\n"
                "voxels and tiles.")

            .def("copyFromArray", &pyGrid::copyFromArray<GridType>,
                py::arg("array"), py::arg("ijk")=py::make_tuple(0, 0, 0),
                py::arg("tolerance")=pyGrid::getZeroValue<GridType>(),
                ("copyFromArray(array, ijk=(0, 0, 0), tolerance=0)\n\n"
                "Populate this grid, starting at voxel (i, j, k), with values\nfrom a "
                + std::string(openvdb::VecTraits<ValueT>::IsVec ? "four" : "three")
                + "-dimensional array.  Mark voxels as inactive\n"
                "if and only if their values are equal to this grid's\n"
                "background value within the given tolerance.").c_str())
            .def("copyToArray", &pyGrid::copyToArray<GridType>,
                py::arg("array"), py::arg("ijk")=py::make_tuple(0, 0, 0),
                ("copyToArray(array, ijk=(0, 0, 0))\n\nPopulate a "
                + std::string(openvdb::VecTraits<ValueT>::IsVec ? "four" : "three")
                + "-dimensional array with values\n"
                "from this grid, starting at voxel (i, j, k).").c_str())
            .def("toArray", &pyGrid::toArray<GridType>,
                py::arg("ijk")=py::make_tuple(0, 0, 0),
                ("toArray(ijk=(0, 0, 0))\n\nPopulate a "
                + std::string(openvdb::VecTraits<ValueT>::IsVec ? "four" : "three")
                + "-dimensional array with values\n"
                "from this grid, starting at voxel (i, j, k).").c_str())

            // .def("convertToQuads",
            //     &pyGrid::volumeToQuadMesh<GridType>,
            //     py::arg("isovalue")=0,
            //     "convertToQuads(isovalue=0) -> points, quads\n\n"
            //     "Uniformly mesh a scalar grid that has a continuous isosurface\n"
            //     "at the given isovalue.  Return a NumPy array of world-space\n"
            //     "points and a NumPy array of 4-tuples of point indices, which\n"
            //     "specify the vertices of the quadrilaterals that form the mesh.")
            // .def("convertToPolygons",
            //     &pyGrid::volumeToMesh<GridType>,
            //     py::arg("isovalue")=0, py::arg("adaptivity")=0,
            //     "convertToPolygons(isovalue=0, adaptivity=0) -> points, triangles, quads\n\n"
            //     "Adaptively mesh a scalar grid that has a continuous isosurface\n"
            //     "at the given isovalue.  Return a NumPy array of world-space\n"
            //     "points and NumPy arrays of 3- and 4-tuples of point indices,\n"
            //     "which specify the vertices of the triangles and quadrilaterals\n"
            //     "that form the mesh.  Adaptivity can vary from 0 to 1, where 0\n"
            //     "produces a high-polygon-count mesh that closely approximates\n"
            //     "the isosurface, and 1 produces a lower-polygon-count mesh\n"
            //     "with some loss of surface detail.")
            // .def_static("createLevelSetFromPolygons",
            //     &pyGrid::meshToLevelSet<GridType>,
            //     py::arg("points"),
            //         py::arg("triangles")=py::none(),
            //         py::arg("quads")=py::none(),
            //         py::arg("transform")=py::none(),
            //         py::arg("halfWidth")=openvdb::LEVEL_SET_HALF_WIDTH,
            //     ("createLevelSetFromPolygons(points, triangles=None, quads=None,\n"
            //      "    transform=None, halfWidth="
            //      + std::to_string(openvdb::LEVEL_SET_HALF_WIDTH) + ") -> "
            //      + pyGridTypeName + "\n\n"
            //     "Convert a triangle and/or quad mesh to a narrow-band level set volume.\n"
            //     "The mesh must form a closed surface, but the surface need not be\n"
            //     "manifold and may have self intersections and degenerate faces.\n"
            //     "The mesh is described by a NumPy array of world-space points\n"
            //     "and NumPy arrays of 3- and 4-tuples of point indices that specify\n"
            //     "the vertices of the triangles and quadrilaterals that form the mesh.\n"
            //     "Either the triangle or the quad array may be empty or None.\n"
            //     "The resulting volume will have the given transform (or the identity\n"
            //     "transform if no transform is given) and a narrow band width of\n"
            //     "2 x halfWidth voxels.").c_str())

            .def("prune", &pyGrid::prune<GridType>,
                py::arg("tolerance")=0,
                "prune(tolerance=0)\n\n"
                "Remove nodes whose values all have the same active state\n"
                "and are equal to within a given tolerance.")
            .def("pruneInactive", &pyGrid::pruneInactive<GridType>,
                py::arg("value")=py::none(),
                "pruneInactive(value=None)\n\n"
                "Remove nodes whose values are all inactive and replace them\n"
                "with either background tiles or tiles of the given value\n"
                "(if the value is not None).")

            .def("empty", &GridType::empty,
                "empty() -> bool\n\n"
                "Return True if this grid contains only background voxels.")
            .def("__nonzero__", &pyGrid::notEmpty<GridType>)

            .def("clear", &GridType::clear,
                "clear()\n\n"
                "Remove all tiles from this grid and all nodes other than the root node.")

            .def("merge", &GridType::merge,
                ("merge(" + pyGridTypeName + ")\n\n"
                "Move child nodes from the other grid into this grid wherever\n"
                "those nodes correspond to constant-value tiles in this grid,\n"
                "and replace leaf-level inactive voxels in this grid with\n"
                "corresponding voxels in the other grid that are active.\n\n"
                "Note: this operation always empties the other grid.").c_str())

            .def("mapOn", &pyGrid::mapOn<GridType>,
                py::arg("function"),
                "mapOn(function)\n\n"
                "Iterate over all the active (\"on\") values (tile and voxel)\n"
                "of this grid and replace each value with function(value).\n\n"
                "Example: grid.mapOn(lambda x: x * 2 if x < 0.5 else x)")

            .def("mapOff", &pyGrid::mapOff<GridType>,
                py::arg("function"),
                "mapOff(function)\n\n"
                "Iterate over all the inactive (\"off\") values (tile and voxel)\n"
                "of this grid and replace each value with function(value).\n\n"
                "Example: grid.mapOff(lambda x: x * 2 if x < 0.5 else x)")

            .def("mapAll", &pyGrid::mapAll<GridType>,
                py::arg("function"),
                "mapAll(function)\n\n"
                "Iterate over all values (tile and voxel) of this grid\n"
                "and replace each value with function(value).\n\n"
                "Example: grid.mapAll(lambda x: x * 2 if x < 0.5 else x)")

            .def("combine", &pyGrid::combine<GridType>,
                py::arg("grid"), py::arg("function"),
                "combine(grid, function)\n\n"
                "Compute function(self, other) over all corresponding pairs\n"
                "of values (tile or voxel) of this grid and the other grid\n"
                "and store the result in this grid.\n\n"
                "Note: this operation always empties the other grid.\n\n"
                "Example: grid.combine(otherGrid, lambda a, b: min(a, b))")

            //
            // Iterators
            //
            .def("citerOnValues", &pyGrid::IterTraits<GridType, ValueOnCIterT>::begin,
                "citerOnValues() -> iterator\n\n"
                "Return a read-only iterator over this grid's active\ntile and voxel values.")
            .def("citerOffValues", &pyGrid::IterTraits<GridType, ValueOffCIterT>::begin,
                "iterOffValues() -> iterator\n\n"
                "Return a read-only iterator over this grid's inactive\ntile and voxel values.")
            .def("citerAllValues", &pyGrid::IterTraits<GridType, ValueAllCIterT>::begin,
                "iterAllValues() -> iterator\n\n"
                "Return a read-only iterator over all of this grid's\ntile and voxel values.")

            .def("iterOnValues", &pyGrid::IterTraits<GridType, ValueOnIterT>::begin,
                "iterOnValues() -> iterator\n\n"
                "Return a read/write iterator over this grid's active\ntile and voxel values.")
            .def("iterOffValues", &pyGrid::IterTraits<GridType, ValueOffIterT>::begin,
                "iterOffValues() -> iterator\n\n"
                "Return a read/write iterator over this grid's inactive\ntile and voxel values.")
            .def("iterAllValues", &pyGrid::IterTraits<GridType, ValueAllIterT>::begin,
                "iterAllValues() -> iterator\n\n"
                "Return a read/write iterator over all of this grid's\ntile and voxel values.")

            ; // py::class_<Grid>

        // Wrap const and non-const value accessors and expose them
        // as nested classes of the Grid class.
        pyAccessor::AccessorWrap<const GridType>::wrap(m);
        pyAccessor::AccessorWrap<GridType>::wrap(m);

        // Wrap tree value iterators and expose them as nested classes of the Grid class.
        IterWrap<const GridType, ValueOnCIterT>::wrap(m);
        IterWrap<const GridType, ValueOffCIterT>::wrap(m);
        IterWrap<const GridType, ValueAllCIterT>::wrap(m);
        IterWrap<GridType, ValueOnIterT>::wrap(m);
        IterWrap<GridType, ValueOffIterT>::wrap(m);
        IterWrap<GridType, ValueAllIterT>::wrap(m);

    } // gridClassScope

    // Add the Python type object for this grid type to the module-level list.
    py::list gridTypes = m.attr("GridTypes");
    gridTypes.append(pyGridTypeName);

    // Python buffer protocol for dense arrays
    py::class_<DenseT>(m, (std::string("Dense") + pyGridTypeName).c_str(), py::buffer_protocol())
    .def_buffer([](DenseT &d) -> py::buffer_info {
        return py::buffer_info(
            d.data(),                                // Pointer to buffer
            sizeof(ValueT),                          // Size of one scalar
            py::format_descriptor<ValueT>::format(), // Python struct-style format descriptor
            3,                                       // Number of dimensions
            { d.bbox().dim()[0],                     // Buffer dimensions
              d.bbox().dim()[1],
              d.bbox().dim()[2] },
            { d.xStride(),                           // Strides (in bytes) for each index
              d.yStride(),
              d.zStride() }
        );
    });
}

} // namespace pyGrid

#endif // OPENVDB_PYGRID_HAS_BEEN_INCLUDED
