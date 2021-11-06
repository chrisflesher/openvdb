// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0
//
/// @file pyIntGrid.cc
/// @brief Boost.Python wrappers for scalar, integer-valued openvdb::Grid types

#include "pyGrid.h"


void exportIntGrid(py::module_ &m);


void
exportIntGrid(py::module_ &m)
{
    pyGrid::exportGrid<BoolGrid>(m);
#ifdef PY_OPENVDB_WRAP_ALL_GRID_TYPES
    pyGrid::exportGrid<Int32Grid>(m);
    pyGrid::exportGrid<Int64Grid>(m);
#endif
}
