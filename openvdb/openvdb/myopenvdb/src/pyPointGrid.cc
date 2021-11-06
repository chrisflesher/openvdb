// Copyright Contributors to the OpenVDB Project
// SPDX-License-Identifier: MPL-2.0
//
/// @file pyPointGrid.cc
/// @brief Boost.Python wrappers for point openvdb::Grid types

#include "pyGrid.h"


void exportPointGrid(py::module_ &m);


void
exportPointGrid(py::module_ &m)
{
#ifdef PY_OPENVDB_WRAP_ALL_GRID_TYPES
    pyGrid::exportGrid<points::PointDataGrid>(m);
#endif
}
