### setup project ###
# https://numpy.org/doc/stable/f2py/buildtools/skbuild.html
cmake_minimum_required(VERSION 3.18)

project(${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  DESCRIPTION "Utilities for reading WRF output"
  LANGUAGES C Fortran
  )

# Safety net
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n"
  )
endif()

find_package(Python 3.8 REQUIRED COMPONENTS Development.Module NumPy)

# Ensure scikit-build modules
if (NOT SKBUILD)
  # Kanged --> https://github.com/Kitware/torch_liberator/blob/master/CMakeLists.txt
  # If skbuild is not the driver; include its utilities in CMAKE_MODULE_PATH
  execute_process(
    COMMAND "${Python_EXECUTABLE}"
    -c "import os, skbuild; print(os.path.dirname(skbuild.__file__))"
    OUTPUT_VARIABLE SKBLD_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  list(APPEND CMAKE_MODULE_PATH "${SKBLD_DIR}/resources/cmake")
  message(STATUS "Looking in ${SKBLD_DIR}/resources/cmake for CMake modules")
endif()

# Grab the variables from a local Python installation
# NumPy headers
# F2PY headers
execute_process(
  COMMAND "${Python_EXECUTABLE}"
  -c "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

find_package(OpenMP COMPONENTS Fortran)
set_source_files_properties(fortran/ompgen.F90 PROPERTIES Fortran_PREPROCESS ON)
# TODO: Figure out the conditionals for running the C Preprocessor on Fortran files
# I think the main thing to be changed is -E -cpp
# Intel is -fpp -save-temps or /fpp on Windows
# or call fpp instead of the fortran compiler to get it to stop after preprocessing
if (${OpenMP_Fortran_FOUND})
  # This would probably be cleaner if I shoved it in the subdirectory
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran")
  add_executable(sizes "${CMAKE_SOURCE_DIR}/fortran/build_help/omp_sizes.f90")
  target_link_libraries(sizes PUBLIC OpenMP::OpenMP_Fortran)
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/fortran/ompgen.F90"
    DEPENDS "${CMAKE_SOURCE_DIR}/fortran/ompgen.F90.template"
      ${CMAKE_SOURCE_DIR}/fortran/build_help/sub_sizes.py
      sizes
    COMMAND
      ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/fortran/build_help/sub_sizes.py
        ${CMAKE_SOURCE_DIR}/fortran/ompgen.F90.template
	${CMAKE_CURRENT_BINARY_DIR}/fortran/ompgen.F90
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/fortran/omp.f90"
    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/fortran/ompgen.F90"
    COMMAND ${CMAKE_Fortran_COMPILER} -E "${CMAKE_CURRENT_BINARY_DIR}/fortran/ompgen.F90"
      -o "${CMAKE_CURRENT_BINARY_DIR}/fortran/omp.f90" ${OpenMP_Fortran_FLAGS} -cpp
    )
else()
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/fortran/omp.f90"
    DEPENDS "${CMAKE_SOURCE_DIR}/fortran/ompgen.F90"
    COMMAND ${CMAKE_Fortran_COMPILER} -E fortran/ompgen.F90 -o fortran/omp.f90 -cpp
    )
endif()

# Prepping the module
set(f2py_module_name "_wrffortran")
set(fortran_src_files
    "${CMAKE_SOURCE_DIR}/fortran/wrf_constants.f90"
    "${CMAKE_SOURCE_DIR}/fortran/wrf_testfunc.f90"
    "${CMAKE_SOURCE_DIR}/fortran/wrf_user.f90"
    "${CMAKE_SOURCE_DIR}/fortran/rip_cape.f90"
    "${CMAKE_SOURCE_DIR}/fortran/wrf_cloud_fracf.f90"
    "${CMAKE_SOURCE_DIR}/fortran/wrf_fctt.f90"
    "${CMAKE_SOURCE_DIR}/fortran/wrf_user_dbz.f90"
    "${CMAKE_SOURCE_DIR}/fortran/wrf_relhl.f90"
    "${CMAKE_SOURCE_DIR}/fortran/calc_uh.f90"
    "${CMAKE_SOURCE_DIR}/fortran/wrf_user_latlon_routines.f90"
    "${CMAKE_SOURCE_DIR}/fortran/wrf_pvo.f90"
    "${CMAKE_SOURCE_DIR}/fortran/eqthecalc.f90"
    "${CMAKE_SOURCE_DIR}/fortran/wrf_rip_phys_routines.f90"
    "${CMAKE_SOURCE_DIR}/fortran/wrf_pw.f90"
    "${CMAKE_SOURCE_DIR}/fortran/wrf_vinterp.f90"
    "${CMAKE_SOURCE_DIR}/fortran/wrf_wind.f90"
    "${CMAKE_CURRENT_BINARY_DIR}/fortran/omp.f90")
set(python_src_files
    "${CMAKE_SOURCE_DIR}/src/wrf/__init__.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/api.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/cache.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/computation.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/config.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/constants.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/contrib.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/coordpair.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/decorators.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/destag.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/extension.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/g_cape.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/g_cloudfrac.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/g_ctt.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/g_dbz.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/g_dewpoint.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/g_geoht.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/g_helicity.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/g_latlon.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/g_omega.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/g_precip.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/g_pressure.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/g_pw.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/g_rh.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/g_slp.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/g_temp.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/g_terrain.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/g_times.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/g_uvmet.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/g_vorticity.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/g_wind.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/geobnds.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/interp.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/interputils.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/latlonutils.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/metadecorators.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/projection.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/projutils.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/py3compat.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/routines.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/specialdec.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/units.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/util.py"
    "${CMAKE_SOURCE_DIR}/src/wrf/version.py"
)
set(f2py_module_c "${f2py_module_name}module.c")

# Target for enforcing dependencies
add_custom_target(genpyf
  DEPENDS "${fortran_src_files}"
)
add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}"
	    "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}-f2pywrappers.f"
     "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}-f2pywrappers2.f90"
  COMMAND ${Python_EXECUTABLE}  -m "numpy.f2py"
                   -m "${f2py_module_name}"
                   --lower # Important
                   ${fortran_src_files}
  DEPENDS "${fortran_src_files}" # Fortran source
)

Python_add_library(${f2py_module_name} MODULE
            "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}"
	    "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}-f2pywrappers.f"
            "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}-f2pywrappers2.f90"
            "${F2PY_INCLUDE_DIR}/fortranobject.c"
            "${fortran_src_files}")

target_include_directories(${f2py_module_name} PUBLIC
                           ${F2PY_INCLUDE_DIR}
                           ${Python_NumPy_INCLUDE_DIRS}
                           ${Python_INCLUDE_DIRS})
set_target_properties(${f2py_module_name} PROPERTIES SUFFIX ".${Python_SOABI}${CMAKE_SHARED_MODULE_SUFFIX}")
set_target_properties(${f2py_module_name} PROPERTIES PREFIX "")

# https://scikit-build-core.readthedocs.io/en/latest/getting_started.html
target_link_libraries(${f2py_module_name} PRIVATE Python::NumPy)
if (${OpenMP_Fortran_FOUND})
  target_link_libraries(${f2py_module_name} PUBLIC OpenMP::OpenMP_Fortran)
endif()

# Linker fixes
if (UNIX)
  if (APPLE)
    set_target_properties(${f2py_module_name} PROPERTIES
    LINK_FLAGS  '-Wl,-dylib,-undefined,dynamic_lookup')
  else()
    set_target_properties(${f2py_module_name} PROPERTIES
  LINK_FLAGS  '-Wl,--allow-shlib-undefined')
  endif()
endif()

add_dependencies(${f2py_module_name} genpyf)

if (NOT SKBUILD)
  string(REGEX REPLACE "^/(usr/(local/)?)?" "" Python_SITEARCH_INSTALL ${Python_SITEARCH})
  string(REGEX REPLACE "^/(usr/(local/)?)?" "" Python_SITELIB_INSTALL ${Python_SITELIB})
  # string(SUBSTRING ${Python_SITEARCH} 1 -1 Python_SITEARCH_INSTALL)
  # string(SUBSTRING ${Python_SITELIB} 1 -1 Python_SITELIB_INSTALL)

  install(TARGETS ${f2py_module_name} DESTINATION "${Python_SITEARCH_INSTALL}/wrf/")
  install(FILES ${python_src_files} DESTINATION "${Python_SITELIB_INSTALL}/wrf/")
  install(FILES src/wrf/data/psadilookup.dat DESTINATION "${Python_SITELIB_INSTALL}/wrf")
else()
  # https://scikit-build-core.readthedocs.io/en/latest/cmakelists.html#install-directories
  install(TARGETS ${f2py_module_name} DESTINATION "${SKBUILD_PLATLIB_DIR}/wrf/")
endif()
