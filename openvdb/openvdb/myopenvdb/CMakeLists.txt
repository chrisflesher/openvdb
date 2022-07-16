cmake_minimum_required(VERSION 3.15...3.19)

project(scikit_build_example VERSION "0.0.1")

if(SKBUILD)
  # Scikit-Build does not add your site-packages to the search path automatically,
  # so we need to add it _or_ the pybind11 specific directory here.
  execute_process(
    COMMAND
      "${PYTHON_EXECUTABLE}" -c
      "import pybind11; print(pybind11.get_cmake_dir())"
    OUTPUT_VARIABLE _tmp_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
  list(APPEND CMAKE_PREFIX_PATH "${_tmp_dir}")
endif()

# Now we can find pybind11
find_package(pybind11 CONFIG REQUIRED)

set(OPENVDB_PYTHON_MODULE_SOURCE_FILES
  src/pyOpenVDBModule.cc
  src/pyFloatGrid.cc
#  src/pyIntGrid.cc
#  src/pyPointGrid.cc
#  src/pyPointGrid.cc
  src/pyTransform.cc
#  src/pyVec3Grid.cc
)

pybind11_add_module(_core MODULE ${OPENVDB_PYTHON_MODULE_SOURCE_FILES})

target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})

target_link_libraries(_core PUBLIC
  openvdb
)

install(TARGETS _core DESTINATION .)
