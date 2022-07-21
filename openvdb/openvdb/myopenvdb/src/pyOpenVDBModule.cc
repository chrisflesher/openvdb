// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0

#include <iostream> // must be included before python on macos
// #include <cstring> // for strncmp(), strrchr(), etc.
// #include <limits>
// #include <string>
// #include <utility> // for std::make_pair()
#include <pybind11/pybind11.h>
#include "openvdb/openvdb.h"
// #include "pyopenvdb.h"
// #include "pyGrid.h"
#include "pyutil.h"

namespace py = pybind11;


// Forward declarations
void exportTransform(py::module_ &m);
// void exportMetadata(py::module_ &m);
void exportFloatGrid(py::module_ &m);
// void exportIntGrid(py::module_ &m);
// void exportVec3Grid(py::module_ &m);
// void exportPointGrid(py::module_ &m);


void print(openvdb::Coord xyz) {
    std::cout << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;
}

// namespace _openvdbmodule {

using namespace openvdb;

// ////////////////////////////////////////

namespace pybind11 { namespace detail {

    /// Helper class to convert between a Python numeric sequence
    /// (tuple, list, etc.) and an openvdb::Coord
    template <> struct type_caster<openvdb::Coord> {
    public:
        PYBIND11_TYPE_CASTER(openvdb::Coord, const_name("Coord"));

        bool load(handle src, bool) {
            PyObject *obj = src.ptr();
            if (!PySequence_Check(obj)) {
                return false;
            }
            switch (PySequence_Length(obj)) {
            case 1:
                value.reset(pyutil::getSequenceItem<openvdb::Int32>(obj, 0));
                break;
            case 3:
                value.reset(
                    pyutil::getSequenceItem<openvdb::Int32>(obj, 0),
                    pyutil::getSequenceItem<openvdb::Int32>(obj, 1),
                    pyutil::getSequenceItem<openvdb::Int32>(obj, 2));
                break;
            default:
                throw py::value_error("expected a sequence of three integers");
                break;
            }
            return !(PyErr_Occurred());
        }

        static handle cast(openvdb::Coord src, return_value_policy, handle) {
            py::object obj = py::make_tuple(src[0], src[1], src[2]);
            Py_INCREF(obj.ptr());
                ///< @todo is this the right way to ensure that the object
                ///< doesn't get freed on exit?
            return obj.ptr();
        }
    };


    /// Helper class to convert between a Python numeric sequence
    /// (tuple, list, etc.) and an openvdb::Vec
    template <> struct type_caster<openvdb::Vec2i> {
    public:
        PYBIND11_TYPE_CASTER(openvdb::Vec2i, const_name("VecT"));

        bool load(handle src, bool) {
            PyObject *obj = src.ptr();
            if (!PySequence_Check(obj)) {
                return false;
            }
            if (PySequence_Length(obj) != openvdb::Vec2i::size) {
                return false;
            }
            for (int n = 0; n < openvdb::Vec2i::size; ++n) {
                value[n] = pyutil::getSequenceItem<openvdb::Vec2i::value_type>(obj, n);
            }
            return !(PyErr_Occurred());
        }

        static handle cast(openvdb::Vec2i src, return_value_policy, handle) {
            OPENVDB_NO_UNREACHABLE_CODE_WARNING_BEGIN
            py::object obj;
            switch (openvdb::Vec2i::size) { // compile-time constant
                case 2: obj = py::make_tuple(src[0], src[1]); break;
                case 3: obj = py::make_tuple(src[0], src[1], src[2]); break;
                case 4: obj = py::make_tuple(src[0], src[1], src[2], src[3]); break;
                default:
                {
                    py::list lst;
                    for (int n = 0; n < openvdb::Vec2i::size; ++n) lst.append(src[n]);
                    obj = lst;
                }
            }
            OPENVDB_NO_UNREACHABLE_CODE_WARNING_END
            Py_INCREF(obj.ptr());
            return obj.ptr();
        }
    };


    /// Helper class to convert between a 2D Python numeric sequence
    /// (tuple, list, etc.) and an openvdb::Mat
    template <> struct type_caster<openvdb::Mat4s> {
    public:
        PYBIND11_TYPE_CASTER(openvdb::Mat4s, const_name("MatT"));

        bool load(handle src, bool) {
            PyObject *obj = src.ptr();
            if (!PySequence_Check(obj)) {
                return false;
            }
            // value = openvdb::Mat4s::zero();
            if (py::len(obj) == openvdb::Mat4s::size) {
                for (int i = 0; i < openvdb::Mat4s::size; ++i) {
                    PyObject *rowObj = PySequence_GetItem(obj, i);
                    if (!PySequence_Check(rowObj)) {
                        return false;
                    }
                    if (py::len(rowObj) != openvdb::Mat4s::size) {
                        return false;
                    }
                    for (int j = 0; j < openvdb::Mat4s::size; ++j) {
                        value(i, j) = py::cast<typename openvdb::Mat4s::value_type>(PySequence_GetItem(rowObj, j));
                    }
                }
            }
            return !(PyErr_Occurred());
        }

        /// Return the given matrix as a Python list of lists.
        static handle cast(openvdb::Mat4s src, return_value_policy, handle) {
            py::list obj;
            for (int i = 0; i < openvdb::Mat4s::size; ++i) {
                py::list rowObj;
                for (int j = 0; j < openvdb::Mat4s::size; ++j) {
                    rowObj.append(src(i, j));
                }
                obj.append(rowObj);
            }
            return std::move(obj);
        }
    };

    /// Helper class to convert between a Python integer and a openvdb::PointIndex
    template <> struct type_caster<PointDataIndex32> {
    public:
        PYBIND11_TYPE_CASTER(PointDataIndex32, const_name("PointIndexT"));

        /// Convert from a Python object to a PointIndex.
        bool load(handle src, bool) {
            PyObject *obj = src.ptr();
            PyObject *tmp = PyNumber_Long(obj);
            if (!tmp) {
                return false;
            }
            value = static_cast<PointDataIndex32>(PyLong_AsLong(tmp));
            Py_DECREF(tmp);
            return !(PyErr_Occurred());
        }

        /// @return a Python integer object equivalent to the given PointIndex.
        static handle cast(PointDataIndex32 src, return_value_policy, handle) {
            return PyLong_FromLong(src);
        }
    };

    // /// Helper class to convert between a Python dict and an openvdb::MetaMap
    // /// @todo Consider implementing a separate, templated converter for
    // /// the various Metadata types.
    // template <> struct type_caster<MetaMap> {
    // public:
    //     PYBIND11_TYPE_CASTER(MetaMap, const_name("MetaMapT"));

    //     /// Convert from a Python object to a PointIndex.
    //     bool load(handle src, bool) {
    //         PyObject *obj = src.ptr();
    //         if (!PyMapping_Check(obj)) {
    //             return false;
    //         }

    //         // Populate the map.
    //         py::dict pyDict(pyutil::pyBorrow(obj));
    //         py::list keys = pyDict.keys();
    //         for (size_t i = 0, N = py::len(keys); i < N; ++i) {
    //             std::string name;
    //             py::object key = keys[i];
    //             if (py::cast<std::string>(key).check()) {
    //                 name = py::cast<std::string>(key);
    //             } else {
    //                 const std::string
    //                     keyAsStr = py::cast<std::string>(key.attr("__str__")()),
    //                     keyType = pyutil::className(key);
    //                 PyErr_Format(PyExc_TypeError,
    //                     "expected string as metadata name, found object"
    //                     " \"%s\" of type %s", keyAsStr.c_str(), keyType.c_str());
    //                 py::throw_error_already_set();
    //             }

    //             // Note: the order of the following tests is significant, as it
    //             // avoids unnecessary type promotion (e.g., of ints to floats).
    //             py::object val = pyDict[keys[i]];
    //             Metadata::Ptr meta;
    //             if (py::cast<std::string>(val).check()) {
    //                 meta.reset(new StringMetadata(py::cast<std::string>(val)));
    //             } else if (bool(PyBool_Check(val.ptr()))) {
    //                 meta.reset(new BoolMetadata(py::cast<bool>(val)));
    //             } else if (py::cast<Int64>(val).check()) {
    //                 const Int64 n = py::cast<Int64>(val);
    //                 if (n <= std::numeric_limits<Int32>::max()
    //                     && n >= std::numeric_limits<Int32>::min())
    //                 {
    //                     meta.reset(new Int32Metadata(static_cast<Int32>(n)));
    //                 } else {
    //                     meta.reset(new Int64Metadata(n));
    //                 }
    //             //} else if (py::cast<float>(val).check()) {
    //             //    meta.reset(new FloatMetadata(py::cast<float>(val)));
    //             } else if (py::cast<double>(val).check()) {
    //                 meta.reset(new DoubleMetadata(py::cast<double>(val)));
    //             } else if (py::cast<Vec2i>(val).check()) {
    //                 meta.reset(new Vec2IMetadata(py::cast<Vec2i>(val)));
    //             } else if (py::cast<Vec2d>(val).check()) {
    //                 meta.reset(new Vec2DMetadata(py::cast<Vec2d>(val)));
    //             } else if (py::cast<Vec2s>(val).check()) {
    //                 meta.reset(new Vec2SMetadata(py::cast<Vec2s>(val)));
    //             } else if (py::cast<Vec3i>(val).check()) {
    //                 meta.reset(new Vec3IMetadata(py::cast<Vec3i>(val)));
    //             } else if (py::cast<Vec3d>(val).check()) {
    //                 meta.reset(new Vec3DMetadata(py::cast<Vec3d>(val)));
    //             } else if (py::cast<Vec3s>(val).check()) {
    //                 meta.reset(new Vec3SMetadata(py::cast<Vec3s>(val)));
    //             } else if (py::cast<Vec4i>(val).check()) {
    //                 meta.reset(new Vec4IMetadata(py::cast<Vec4i>(val)));
    //             } else if (py::cast<Vec4d>(val).check()) {
    //                 meta.reset(new Vec4DMetadata(py::cast<Vec4d>(val)));
    //             } else if (py::cast<Vec4s>(val).check()) {
    //                 meta.reset(new Vec4SMetadata(py::cast<Vec4s>(val)));
    //             } else if (py::cast<Mat4d>(val).check()) {
    //                 meta.reset(new Mat4DMetadata(py::cast<Mat4d>(val)));
    //             } else if (py::cast<Mat4s>(val).check()) {
    //                 meta.reset(new Mat4SMetadata(py::cast<Mat4s>(val)));
    //             } else if (py::cast<Metadata::Ptr>(val).check()) {
    //                 meta = py::cast<Metadata::Ptr>(val);
    //             } else {
    //                 const std::string
    //                     valAsStr = py::cast<std::string>(val.attr("__str__")()),
    //                     valType = pyutil::className(val);
    //                 PyErr_Format(PyExc_TypeError,
    //                     "metadata value \"%s\" of type %s is not allowed",
    //                     valAsStr.c_str(), valType.c_str());
    //                 py::throw_error_already_set();
    //             }
    //             if (meta) {
    //                 value.insertMeta(name, *meta);
    //             }
    //         }
    //     }

    //     /// @return a Python dict object equivalent to the given MetaMap.
    //     static handle cast(MetaMap src, return_value_policy, handle) {
    //         py::dict ret;
    //         for (MetaMap::ConstMetaIterator it = src.beginMeta(); it != src.endMeta(); ++it) {
    //             if (Metadata::Ptr meta = it->second) {
    //                 py::object obj(meta);
    //                 const std::string typeName = meta->typeName();
    //                 if (typeName == StringMetadata::staticTypeName()) {
    //                     obj = py::str(static_cast<StringMetadata&>(*meta).value());
    //                 } else if (typeName == DoubleMetadata::staticTypeName()) {
    //                     obj = py::object(static_cast<DoubleMetadata&>(*meta).value());
    //                 } else if (typeName == FloatMetadata::staticTypeName()) {
    //                     obj = py::object(static_cast<FloatMetadata&>(*meta).value());
    //                 } else if (typeName == Int32Metadata::staticTypeName()) {
    //                     obj = py::object(static_cast<Int32Metadata&>(*meta).value());
    //                 } else if (typeName == Int64Metadata::staticTypeName()) {
    //                     obj = py::object(static_cast<Int64Metadata&>(*meta).value());
    //                 } else if (typeName == BoolMetadata::staticTypeName()) {
    //                     obj = py::object(static_cast<BoolMetadata&>(*meta).value());
    //                 } else if (typeName == Vec2DMetadata::staticTypeName()) {
    //                     const Vec2d v = static_cast<Vec2DMetadata&>(*meta).value();
    //                     obj = py::make_tuple(v[0], v[1]);
    //                 } else if (typeName == Vec2IMetadata::staticTypeName()) {
    //                     const Vec2i v = static_cast<Vec2IMetadata&>(*meta).value();
    //                     obj = py::make_tuple(v[0], v[1]);
    //                 } else if (typeName == Vec2SMetadata::staticTypeName()) {
    //                     const Vec2s v = static_cast<Vec2SMetadata&>(*meta).value();
    //                     obj = py::make_tuple(v[0], v[1]);
    //                 } else if (typeName == Vec3DMetadata::staticTypeName()) {
    //                     const Vec3d v = static_cast<Vec3DMetadata&>(*meta).value();
    //                     obj = py::make_tuple(v[0], v[1], v[2]);
    //                 } else if (typeName == Vec3IMetadata::staticTypeName()) {
    //                     const Vec3i v = static_cast<Vec3IMetadata&>(*meta).value();
    //                     obj = py::make_tuple(v[0], v[1], v[2]);
    //                 } else if (typeName == Vec3SMetadata::staticTypeName()) {
    //                     const Vec3s v = static_cast<Vec3SMetadata&>(*meta).value();
    //                     obj = py::make_tuple(v[0], v[1], v[2]);
    //                 } else if (typeName == Vec4DMetadata::staticTypeName()) {
    //                     const Vec4d v = static_cast<Vec4DMetadata&>(*meta).value();
    //                     obj = py::make_tuple(v[0], v[1], v[2], v[3]);
    //                 } else if (typeName == Vec4IMetadata::staticTypeName()) {
    //                     const Vec4i v = static_cast<Vec4IMetadata&>(*meta).value();
    //                     obj = py::make_tuple(v[0], v[1], v[2], v[3]);
    //                 } else if (typeName == Vec4SMetadata::staticTypeName()) {
    //                     const Vec4s v = static_cast<Vec4SMetadata&>(*meta).value();
    //                     obj = py::make_tuple(v[0], v[1], v[2], v[3]);
    //                 } else if (typeName == Mat4SMetadata::staticTypeName()) {
    //                     const Mat4s m = static_cast<Mat4SMetadata&>(*meta).value();
    //                     obj = MatConverter<Mat4s>::toList(m);
    //                 } else if (typeName == Mat4DMetadata::staticTypeName()) {
    //                     const Mat4d m = static_cast<Mat4DMetadata&>(*meta).value();
    //                     obj = MatConverter<Mat4d>::toList(m);
    //                 }
    //                 ret[it->first] = obj;
    //             }
    //         }
    //         Py_INCREF(ret.ptr());
    //         return ret.ptr();
    //     }
    // };

}} // namespace pybind11::detail


// ////////////////////////////////////////


// template<typename T> void translateException(const T&) {}

// /// @brief Define a function that translates an OpenVDB exception into
// /// the equivalent Python exception.
// /// @details openvdb::Exception::what() typically returns a string of the form
// /// "<exception>: <description>".  To avoid duplication of the exception name in Python
// /// stack traces, the function strips off the "<exception>: " prefix.  To do that,
// /// it needs the class name in the form of a string, hence the preprocessor macro.
// #define PYOPENVDB_CATCH(_openvdbname, _pyname)                      \
//     template<>                                                      \
//     void translateException<_openvdbname>(const _openvdbname& e)    \
//     {                                                               \
//         const char* name = #_openvdbname;                           \
//         if (const char* c = std::strrchr(name, ':')) name = c + 1;  \
//         const int namelen = int(std::strlen(name));                 \
//         const char* msg = e.what();                                 \
//         if (0 == std::strncmp(msg, name, namelen)) msg += namelen;  \
//         if (0 == std::strncmp(msg, ": ", 2)) msg += 2;              \
//         PyErr_SetString(_pyname, msg);                              \
//     }


// /// Define an overloaded function that translate all OpenVDB exceptions into
// /// their Python equivalents.
// /// @todo LookupError is redundant and should someday be removed.
// PYOPENVDB_CATCH(openvdb::ArithmeticError,       PyExc_ArithmeticError)
// PYOPENVDB_CATCH(openvdb::IndexError,            PyExc_IndexError)
// PYOPENVDB_CATCH(openvdb::IoError,               PyExc_IOError)
// PYOPENVDB_CATCH(openvdb::KeyError,              PyExc_KeyError)
// PYOPENVDB_CATCH(openvdb::LookupError,           PyExc_LookupError)
// PYOPENVDB_CATCH(openvdb::NotImplementedError,   PyExc_NotImplementedError)
// PYOPENVDB_CATCH(openvdb::ReferenceError,        PyExc_ReferenceError)
// PYOPENVDB_CATCH(openvdb::RuntimeError,          PyExc_RuntimeError)
// PYOPENVDB_CATCH(openvdb::TypeError,             PyExc_TypeError)
// PYOPENVDB_CATCH(openvdb::ValueError,            PyExc_ValueError)

// #undef PYOPENVDB_CATCH


// ////////////////////////////////////////


// py::object readFromFile(const std::string&, const std::string&);
// py::tuple readAllFromFile(const std::string&);
// py::dict readFileMetadata(const std::string&);
// py::object readGridMetadataFromFile(const std::string&, const std::string&);
// py::list readAllGridMetadataFromFile(const std::string&);
// void writeToFile(const std::string&, py::object, py::object);


// py::object
// readFromFile(const std::string& filename, const std::string& gridName)
// {
//     io::File vdbFile(filename);
//     vdbFile.open();

//     if (!vdbFile.hasGrid(gridName)) {
//         PyErr_Format(PyExc_KeyError,
//             "file %s has no grid named \"%s\"",
//             filename.c_str(), gridName.c_str());
//         py::throw_error_already_set();
//     }

//     return pyGrid::getGridFromGridBase(vdbFile.readGrid(gridName));
// }


// py::tuple
// readAllFromFile(const std::string& filename)
// {
//     io::File vdbFile(filename);
//     vdbFile.open();

//     GridPtrVecPtr grids = vdbFile.getGrids();
//     MetaMap::Ptr metadata = vdbFile.getMetadata();
//     vdbFile.close();

//     py::list gridList;
//     for (GridPtrVec::const_iterator it = grids->begin(); it != grids->end(); ++it) {
//         gridList.append(pyGrid::getGridFromGridBase(*it));
//     }

//     return py::make_tuple(gridList, py::dict(*metadata));
// }


// py::dict
// readFileMetadata(const std::string& filename)
// {
//     io::File vdbFile(filename);
//     vdbFile.open();

//     MetaMap::Ptr metadata = vdbFile.getMetadata();
//     vdbFile.close();

//     return py::dict(*metadata);
// }


// py::object
// readGridMetadataFromFile(const std::string& filename, const std::string& gridName)
// {
//     io::File vdbFile(filename);
//     vdbFile.open();

//     if (!vdbFile.hasGrid(gridName)) {
//         PyErr_Format(PyExc_KeyError,
//             "file %s has no grid named \"%s\"",
//             filename.c_str(), gridName.c_str());
//         py::throw_error_already_set();
//     }

//     return pyGrid::getGridFromGridBase(vdbFile.readGridMetadata(gridName));
// }


// py::list
// readAllGridMetadataFromFile(const std::string& filename)
// {
//     io::File vdbFile(filename);
//     vdbFile.open();
//     GridPtrVecPtr grids = vdbFile.readAllGridMetadata();
//     vdbFile.close();

//     py::list gridList;
//     for (GridPtrVec::const_iterator it = grids->begin(); it != grids->end(); ++it) {
//         gridList.append(pyGrid::getGridFromGridBase(*it));
//     }
//     return gridList;
// }


// void
// writeToFile(const std::string& filename, py::object gridOrSeqObj, py::object dictObj)
// {
//     GridPtrVec gridVec;
//     try {
//         GridBase::Ptr base = pyopenvdb::getGridFromPyObject(gridOrSeqObj);
//         gridVec.push_back(base);
//     } catch (openvdb::TypeError&) {
//         for (py::stl_input_iterator<py::object> it(gridOrSeqObj), end; it != end; ++it) {
//             if (GridBase::Ptr base = pyGrid::getGridBaseFromGrid(*it)) {
//                 gridVec.push_back(base);
//             }
//         }
//     }

//     io::File vdbFile(filename);
//     if (dictObj.is_none()) {
//         vdbFile.write(gridVec);
//     } else {
//         MetaMap metadata = py::cast<MetaMap>(dictObj);
//         vdbFile.write(gridVec, metadata);
//     }
//     vdbFile.close();
// }


// ////////////////////////////////////////


// std::string getLoggingLevel();
// void setLoggingLevel(py::object);
// void setProgramName(py::object, bool);


std::string
getLoggingLevel()
{
    switch (logging::getLevel()) {
        case logging::Level::Debug: return "debug";
        case logging::Level::Info:  return "info";
        case logging::Level::Warn:  return "warn";
        case logging::Level::Error: return "error";
        case logging::Level::Fatal: break;
    }
    return "fatal";
}


void
setLoggingLevel(std::string level)
{
    if (level == "debug") { logging::setLevel(logging::Level::Debug); }
    else if (level == "info") { logging::setLevel(logging::Level::Info); }
    else if (level == "warn") { logging::setLevel(logging::Level::Warn); }
    else if (level == "error") { logging::setLevel(logging::Level::Error); }
    else if (level == "fatal") { logging::setLevel(logging::Level::Fatal); }
    else {
        throw py::value_error(
            "expected logging level \"debug\", \"info\", \"warn\", \"error\", or \"fatal\"");
    }
}


void
setProgramName(std::string name, bool color)
{
    logging::setProgramName(name, color);
}


// ////////////////////////////////////////


// // Descriptor for the openvdb::GridClass enum (for use with pyutil::StringEnum)
// struct GridClassDescr
// {
//     static const char* name() { return "GridClass"; }
//     static const char* doc()
//     {
//         return "Classes of volumetric data (level set, fog volume, etc.)";
//     }
//     static pyutil::CStringPair item(int i)
//     {
//         static const int sCount = 4;
//         static const char* const sStrings[sCount][2] = {
//             { "UNKNOWN",    strdup(GridBase::gridClassToString(GRID_UNKNOWN).c_str()) },
//             { "LEVEL_SET",  strdup(GridBase::gridClassToString(GRID_LEVEL_SET).c_str()) },
//             { "FOG_VOLUME", strdup(GridBase::gridClassToString(GRID_FOG_VOLUME).c_str()) },
//             { "STAGGERED",  strdup(GridBase::gridClassToString(GRID_STAGGERED).c_str()) }
//         };
//         if (i >= 0 && i < sCount) return pyutil::CStringPair(&sStrings[i][0], &sStrings[i][1]);
//         return pyutil::CStringPair(static_cast<char**>(nullptr), static_cast<char**>(nullptr));
//     }
// };


// // Descriptor for the openvdb::VecType enum (for use with pyutil::StringEnum)
// struct VecTypeDescr
// {
//     static const char* name() { return "VectorType"; }
//     static const char* doc()
//     {
//         return
//             "The type of a vector determines how transforms are applied to it.\n"
//             "  - INVARIANT:\n"
//             "      does not transform (e.g., tuple, uvw, color)\n"
//             "  - COVARIANT:\n"
//             "      apply inverse-transpose transformation with w = 0\n"
//             "      and ignore translation (e.g., gradient/normal)\n"
//             "  - COVARIANT_NORMALIZE:\n"
//             "      apply inverse-transpose transformation with w = 0\n"
//             "      and ignore translation, vectors are renormalized\n"
//             "      (e.g., unit normal)\n"
//             "  - CONTRAVARIANT_RELATIVE:\n"
//             "      apply \"regular\" transformation with w = 0 and ignore\n"
//             "      translation (e.g., displacement, velocity, acceleration)\n"
//             "  - CONTRAVARIANT_ABSOLUTE:\n"
//             "      apply \"regular\" transformation with w = 1 so that\n"
//             "      vector translates (e.g., position)\n";
//     }
//     static pyutil::CStringPair item(int i)
//     {
//         static const int sCount = 5;
//         static const char* const sStrings[sCount][2] = {
//             { "INVARIANT", strdup(GridBase::vecTypeToString(openvdb::VEC_INVARIANT).c_str()) },
//             { "COVARIANT", strdup(GridBase::vecTypeToString(openvdb::VEC_COVARIANT).c_str()) },
//             { "COVARIANT_NORMALIZE",
//                 strdup(GridBase::vecTypeToString(openvdb::VEC_COVARIANT_NORMALIZE).c_str()) },
//             { "CONTRAVARIANT_RELATIVE",
//                 strdup(GridBase::vecTypeToString(openvdb::VEC_CONTRAVARIANT_RELATIVE).c_str()) },
//             { "CONTRAVARIANT_ABSOLUTE",
//                 strdup(GridBase::vecTypeToString(openvdb::VEC_CONTRAVARIANT_ABSOLUTE).c_str()) }
//         };
//         if (i >= 0 && i < sCount) return std::make_pair(&sStrings[i][0], &sStrings[i][1]);
//         return pyutil::CStringPair(static_cast<char**>(nullptr), static_cast<char**>(nullptr));
//     }
// };

// } // namespace _openvdbmodule


////////////////////////////////////////


// #ifdef DWA_OPENVDB
// #define PY_OPENVDB_MODULE_NAME  _openvdb
// extern "C" { void init_openvdb(); }
// #else
// #define PY_OPENVDB_MODULE_NAME  pyopenvdb
// extern "C" { void initpyopenvdb(); }
// #endif

PYBIND11_MODULE(_core, m) // PY_OPENVDB_MODULE_NAME
{
    // // Don't auto-generate ugly, C++-style function signatures.
    // py::docstring_options docOptions;
    // docOptions.disable_signatures();
    // docOptions.enable_user_defined();

// #ifdef PY_OPENVDB_USE_NUMPY
//     // Initialize NumPy.
// #ifdef PY_OPENVDB_USE_BOOST_PYTHON_NUMPY
//     boost::python::numpy::initialize();
// #else
// #if PY_MAJOR_VERSION >= 3
//     if (_import_array()) {}
// #else
//     import_array();
// #endif
// #endif
// #endif

    // using namespace openvdb::OPENVDB_VERSION_NAME;

    // Initialize OpenVDB.
    initialize();

    // _openvdbmodule::VecConverter<Vec2i>::registerConverter();
    // _openvdbmodule::VecConverter<Vec2I>::registerConverter();
    // _openvdbmodule::VecConverter<Vec2s>::registerConverter();
    // _openvdbmodule::VecConverter<Vec2d>::registerConverter();

    // _openvdbmodule::VecConverter<Vec3i>::registerConverter();
    // _openvdbmodule::VecConverter<Vec3I>::registerConverter();
    // _openvdbmodule::VecConverter<Vec3s>::registerConverter();
    // _openvdbmodule::VecConverter<Vec3d>::registerConverter();

    // _openvdbmodule::VecConverter<Vec4i>::registerConverter();
    // _openvdbmodule::VecConverter<Vec4I>::registerConverter();
    // _openvdbmodule::VecConverter<Vec4s>::registerConverter();
    // _openvdbmodule::VecConverter<Vec4d>::registerConverter();

    // _openvdbmodule::MatConverter<Mat4s>::registerConverter();
    // _openvdbmodule::MatConverter<Mat4d>::registerConverter();

    // _openvdbmodule::PointIndexConverter<PointDataIndex32>::registerConverter();

    // _openvdbmodule::MetaMapConverter::registerConverter();

// #define PYOPENVDB_TRANSLATE_EXCEPTION(_classname) \
//     py::register_exception_translator<_classname>(&_openvdbmodule::translateException<_classname>)

//     PYOPENVDB_TRANSLATE_EXCEPTION(ArithmeticError);
//     PYOPENVDB_TRANSLATE_EXCEPTION(IndexError);
//     PYOPENVDB_TRANSLATE_EXCEPTION(IoError);
//     PYOPENVDB_TRANSLATE_EXCEPTION(KeyError);
//     PYOPENVDB_TRANSLATE_EXCEPTION(LookupError);
//     PYOPENVDB_TRANSLATE_EXCEPTION(NotImplementedError);
//     PYOPENVDB_TRANSLATE_EXCEPTION(ReferenceError);
//     PYOPENVDB_TRANSLATE_EXCEPTION(RuntimeError);
//     PYOPENVDB_TRANSLATE_EXCEPTION(TypeError);
//     PYOPENVDB_TRANSLATE_EXCEPTION(ValueError);

// #undef PYOPENVDB_TRANSLATE_EXCEPTION


    // Export the python bindings.
    exportTransform(m);
    // exportMetadata(m);
    exportFloatGrid(m);
    // exportIntGrid(m);
    // exportVec3Grid(m);
    // exportPointGrid(m);

    m.def("print",
          &print,
          py::arg("xyz"),
          "print(xyz) -> None\n\n"
          "Print yo int yo.");

    // m.def("read",
    //     &_openvdbmodule::readFromFile,
    //     (py::arg("filename"), py::arg("gridname")),
    //     "read(filename, gridname) -> Grid\n\n"
    //     "Read a single grid from a .vdb file.");

   //  m.def("readAll",
   //      &_openvdbmodule::readAllFromFile,
   //      py::arg("filename"),
   //      "readAll(filename) -> list, dict\n\n"
   //      "Read a .vdb file and return a list of grids and\n"
   //      "a dict of file-level metadata.");

   //  m.def("readMetadata",
   //      &_openvdbmodule::readFileMetadata,
   //      py::arg("filename"),
   //      "readMetadata(filename) -> dict\n\n"
   //      "Read file-level metadata from a .vdb file.");

   //  m.def("readGridMetadata",
   //      &_openvdbmodule::readGridMetadataFromFile,
   //      (m.arg("filename"), py::arg("gridname")),
   //      "readGridMetadata(filename, gridname) -> Grid\n\n"
   //      "Read a single grid's metadata and transform (but not its tree)\n"
   //      "from a .vdb file.");

   //  m.def("readAllGridMetadata",
   //      &_openvdbmodule::readAllGridMetadataFromFile,
   //      py::arg("filename"),
   //      "readAllGridMetadata(filename) -> list\n\n"
   //      "Read a .vdb file and return a list of grids populated with\n"
   //      "their metadata and transforms, but not their trees.");

   //  m.def("write",
   //      &_openvdbmodule::writeToFile,
   //      (py::arg("filename"), py::arg("grids"), py::arg("metadata") = py::object()),
   //      "write(filename, grids, metadata=None)\n\n"
   //      "Write a grid or a sequence of grids and, optionally, a dict\n"
   //      "of (name, value) metadata pairs to a .vdb file.");

   m.def("getLoggingLevel", &getLoggingLevel,
        "getLoggingLevel() -> str\n\n"
        "Return the severity threshold (\"debug\", \"info\", \"warn\", \"error\",\n"
        "or \"fatal\") for error messages.");
    m.def("setLoggingLevel", &setLoggingLevel,
        (py::arg("level")),
        "setLoggingLevel(level)\n\n"
        "Specify the severity threshold (\"debug\", \"info\", \"warn\", \"error\",\n"
        "or \"fatal\") for error messages.  Messages of lower severity\n"
        "will be suppressed.");
    m.def("setProgramName", &setProgramName,
        py::arg("name"), py::arg("color") = true,
        "setProgramName(name, color=True)\n\n"
        "Specify the program name to be displayed in error messages,\n"
        "and optionally specify whether to print error messages in color.");

   //  // Add some useful module-level constants.
   //  py::scope().attr("LIBRARY_VERSION") = py::make_tuple(
   //      openvdb::OPENVDB_LIBRARY_MAJOR_VERSION,
   //      openvdb::OPENVDB_LIBRARY_MINOR_VERSION,
   //      openvdb::OPENVDB_LIBRARY_PATCH_VERSION);
   //  py::scope().attr("FILE_FORMAT_VERSION") = openvdb::OPENVDB_FILE_VERSION;
   //  py::scope().attr("COORD_MIN") = openvdb::Coord::min();
   //  py::scope().attr("COORD_MAX") = openvdb::Coord::max();
   //  py::scope().attr("LEVEL_SET_HALF_WIDTH") = openvdb::LEVEL_SET_HALF_WIDTH;

   //  pyutil::StringEnum<_openvdbmodule::GridClassDescr>::wrap();
   //  pyutil::StringEnum<_openvdbmodule::VecTypeDescr>::wrap();

} // PYBIND11_MODULE
