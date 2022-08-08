#include <iostream>
#include <pybind11/pybind11.h>
#include "openvdb/openvdb.h"

namespace py = pybind11;

namespace pybind11 { namespace detail {

    /// Helper class to convert between a Python numeric sequence
    /// (tuple, list, etc.) and an openvdb::Coord
    template <> struct type_caster<openvdb::Coord> {
    public:
        PYBIND11_TYPE_CASTER(openvdb::Coord, const_name("Coord"));

        bool load(handle src, bool) {
            py::sequence obj;
            try {
                obj = py::cast<py::sequence>(src);  // HACK: this performs a copy... maybe use reinterpret_steal?
            } catch (py::cast_error) {
                return false;
            }
            switch (py::len(obj)) {
            case 1:
                value.reset(py::cast<openvdb::Int32>(obj[0]));
                break;
            case 3:
                value.reset(py::cast<openvdb::Int32>(obj[0]),
                            py::cast<openvdb::Int32>(obj[1]),
                            py::cast<openvdb::Int32>(obj[2]));
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

}} // namespace pybind11::detail


void hello() {
    std::cout << "Hello, World!" << std::endl;
}

int return_two() {
    return 2;
}

openvdb::Coord pass_coord(openvdb::Coord coord) {
    return coord;
}


PYBIND11_MODULE(_hello, m) {
    m.doc() = "_hello";
    m.def("hello", &hello, "Prints \"Hello, World!\"");
    m.def("return_two", &return_two, "Returns 2");
    m.def("pass_coord", &pass_coord, py::arg("coord")=openvdb::Coord(0), "Passthrough coord");
}
